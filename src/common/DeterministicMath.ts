/*
 * Planck.js — deterministic math fork
 *
 * Copyright (c) Erin Catto, Ali Shakiba
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 *
 * ---------------------------------------------------------------------------
 *
 * WHY THIS MODULE EXISTS
 *
 * Upstream planck calls the host `Math.sin`, `Math.cos`, `Math.atan2` and
 * `Math.random`. IEEE-754 guarantees `+ - * /` and `sqrt` are correctly
 * rounded, so those are identical on every platform — but ECMAScript
 * explicitly permits an implementation-dependent approximation for the
 * transcendental functions, and every engine forwards them to the host libm.
 * Those results differ between ARM64/JavaScriptCore and x86-64/V8.
 *
 * `Rot.setAngle` and `Sweep.getTransform` call sin/cos for every rotating body
 * on every step, so a simulation shared between a mobile build and a desktop
 * build diverges as soon as anything rotates. That breaks replay portability,
 * server-side replay verification, and any leaderboard built on them.
 *
 * Everything here is therefore built from `+ - * /` and comparisons only, and
 * is bit-identical on every platform.
 *
 * Polynomial coefficients and the Cody-Waite range reduction follow fdlibm
 * (Sun Microsystems, freely distributable), the reference implementation most
 * libms derive from. Accuracy is within a few ulp of the host libm — far
 * tighter than the tolerance of any physical quantity a simulation carries.
 */

/** @internal */ const PI = 3.141592653589793;
/** @internal */ const TWO_OVER_PI = 0.6366197723675814;
/** @internal */ const PIO2_1 = 1.5707963267341256;
/** @internal */ const PIO2_1T = 6.077100506506192e-11;
/** @internal */ const LN2_HI = 0.6931471803691238;
/** @internal */ const LN2_LO = 1.9082149292705877e-10;
/** @internal */ const INV_LN2 = 1.4426950408889634;
/** @internal */ const HALF_PI = 1.5707963267948966;

// fdlibm __kernel_sin coefficients (minimax over [-PI/4, PI/4]).
/** @internal */ const S1 = -1.6666666666666632e-1;
/** @internal */ const S2 = 8.333333333332249e-3;
/** @internal */ const S3 = -1.984126982985795e-4;
/** @internal */ const S4 = 2.755731370707079e-6;
/** @internal */ const S5 = -2.5050760253406863e-8;
/** @internal */ const S6 = 1.58969099521155e-10;

// fdlibm __kernel_cos coefficients (minimax over [-PI/4, PI/4]).
/** @internal */ const C1 = 4.166666666666602e-2;
/** @internal */ const C2 = -1.3888888888741097e-3;
/** @internal */ const C3 = 2.480158728947673e-5;
/** @internal */ const C4 = -2.7557314351390663e-7;
/** @internal */ const C5 = 2.0875723212981748e-9;
/** @internal */ const C6 = -1.1359647557788195e-11;

// fdlibm atan coefficients.
/** @internal */ const AT0 = 3.333333333332932e-1;
/** @internal */ const AT1 = -1.9999999999876483e-1;
/** @internal */ const AT2 = 1.4285714272503466e-1;
/** @internal */ const AT3 = -1.1111110405462356e-1;
/** @internal */ const AT4 = 9.090887133436507e-2;
/** @internal */ const AT5 = -7.691876205044829e-2;
/** @internal */ const AT6 = 6.661073137387531e-2;
/** @internal */ const AT7 = -5.833570133790573e-2;
/** @internal */ const AT8 = 4.976877994615932e-2;
/** @internal */ const AT9 = -3.6531572744216916e-2;
/** @internal */ const AT10 = 1.6285820115365782e-2;

/** @internal */ const ATAN_HI = [
  4.636476090008061e-1, 7.853981633974483e-1, 9.827937232473291e-1, 1.5707963267948966,
];
/** @internal */ const ATAN_LO = [
  2.2698777452961687e-17, 3.061616997868383e-17, 1.3903311031230998e-17,
  6.123233995736766e-17,
];

/** @internal */
function isFiniteNumber(x: number): boolean {
  return x === x && x !== Infinity && x !== -Infinity;
}

/**
 * Round to the nearest integer using only arithmetic.
 *
 * Adding then subtracting 2^52 forces integer rounding under IEEE-754
 * round-to-nearest-even without any library call.
 * @internal
 */
function roundToInt(x: number): number {
  const MAGIC = 4503599627370496; // 2^52
  if (x >= 0) {
    if (x >= MAGIC) return x;
    return x + MAGIC - MAGIC;
  }
  if (-x >= MAGIC) return x;
  return x - MAGIC + MAGIC;
}

/** sin on the reduced interval [-PI/4, PI/4]. @internal */
function kernelSin(x: number): number {
  const z = x * x;
  const r = S2 + z * (S3 + z * (S4 + z * (S5 + z * S6)));
  return x + x * z * (S1 + z * r);
}

/** cos on the reduced interval [-PI/4, PI/4]. @internal */
function kernelCos(x: number): number {
  const z = x * x;
  const r = z * z * (C1 + z * (C2 + z * (C3 + z * (C4 + z * (C5 + z * C6)))));
  return 1 - 0.5 * z + r;
}

/**
 * Reduce x to r in [-PI/4, PI/4] plus a quadrant index, so x = n*(PI/2) + r.
 *
 * Two-term Cody-Waite: `n * PIO2_1` is exact because PIO2_1's low 33 bits are
 * zero, and PIO2_1 + PIO2_1T reproduces PI/2 to about 1e-27. That keeps the
 * reduction accurate even for the large angles a long-running body accumulates.
 * @internal
 */
function reduceQuadrant(x: number): { n: number; r: number } {
  const n = roundToInt(x * TWO_OVER_PI);
  const r = x - n * PIO2_1 - n * PIO2_1T;
  let q = n - roundToInt(n / 4) * 4;
  if (q < 0) q += 4;
  return { n: q, r };
}

/** Cross-platform deterministic sine. */
export function sin(x: number): number {
  if (!isFiniteNumber(x)) return NaN;
  const { n, r } = reduceQuadrant(x);
  switch (n) {
    case 0:
      return kernelSin(r);
    case 1:
      return kernelCos(r);
    case 2:
      return -kernelSin(r);
    default:
      return -kernelCos(r);
  }
}

/** Cross-platform deterministic cosine. */
export function cos(x: number): number {
  if (!isFiniteNumber(x)) return NaN;
  const { n, r } = reduceQuadrant(x);
  switch (n) {
    case 0:
      return kernelCos(r);
    case 1:
      return -kernelSin(r);
    case 2:
      return -kernelCos(r);
    default:
      return kernelSin(r);
  }
}

/** 2^k for |k| <= 500, by exact repeated multiplication. @internal */
function twoTo(k: number): number {
  let result = 1;
  const factor = k >= 0 ? 2 : 0.5;
  const count = k >= 0 ? k : -k;
  for (let i = 0; i < count; i++) result = result * factor;
  return result;
}

/** Exact scaling by a power of two — no rounding occurs. @internal */
function scaleByPowerOfTwo(x: number, k: number): number {
  let result = x;
  let n = k;
  while (n > 0) {
    const step = n > 500 ? 500 : n;
    result = result * twoTo(step);
    n -= step;
  }
  while (n < 0) {
    const step = n < -500 ? -500 : n;
    result = result * twoTo(step);
    n -= step;
  }
  return result;
}

/**
 * Cross-platform deterministic exponential.
 *
 * Range-reduce to x = k*ln2 + r with |r| <= ln2/2, evaluate e^r by Taylor
 * series, then scale by 2^k exactly.
 */
export function exp(x: number): number {
  if (!isFiniteNumber(x)) return x > 0 ? Infinity : NaN;
  if (x > 709.78) return Infinity;
  if (x < -745.2) return 0;

  const k = roundToInt(x * INV_LN2);
  const hi = x - k * LN2_HI;
  const lo = k * LN2_LO;
  const r = hi - lo;

  let term = 1;
  let sum = 1;
  for (let i = 1; i <= 15; i++) {
    term = (term * r) / i;
    sum = sum + term;
  }

  return scaleByPowerOfTwo(sum, k);
}

/** Cross-platform deterministic arctangent. */
export function atan(x: number): number {
  if (!isFiniteNumber(x)) {
    if (x === Infinity) return HALF_PI;
    if (x === -Infinity) return -HALF_PI;
    return NaN;
  }
  const negative = x < 0;
  let v = negative ? -x : x;
  let id: number;

  if (v < 0.4375) {
    id = -1;
  } else if (v < 0.6875) {
    id = 0;
    v = (2 * v - 1) / (2 + v);
  } else if (v < 1.1875) {
    id = 1;
    v = (v - 1) / (v + 1);
  } else if (v < 2.4375) {
    id = 2;
    v = (v - 1.5) / (1 + 1.5 * v);
  } else {
    id = 3;
    v = -1 / v;
  }

  const z = v * v;
  const w = z * z;
  const s1 = z * (AT0 + w * (AT2 + w * (AT4 + w * (AT6 + w * (AT8 + w * AT10)))));
  const s2 = w * (AT1 + w * (AT3 + w * (AT5 + w * (AT7 + w * AT9))));

  let result: number;
  if (id < 0) {
    result = v - v * (s1 + s2);
  } else {
    const hi = ATAN_HI[id] as number;
    const lo = ATAN_LO[id] as number;
    result = hi - (v * (s1 + s2) - lo - v);
  }
  return negative ? -result : result;
}

/** Cross-platform deterministic two-argument arctangent. */
export function atan2(y: number, x: number): number {
  if (x === 0 && y === 0) return 0;
  if (x === 0) return y > 0 ? HALF_PI : -HALF_PI;
  if (x > 0) return atan(y / x);
  return y >= 0 ? atan(y / x) + PI : atan(y / x) - PI;
}

/**
 * Square root.
 *
 * The one host call retained. IEEE-754 §5.4.1 requires `sqrt` to be correctly
 * rounded and every engine compiles it to the hardware instruction, so it is
 * bit-identical across platforms in a way `sin`/`exp` are not.
 */
export function sqrt(x: number): number {
  return Math.sqrt(x);
}

/** @internal */ let randomState = 0x9e3779b9;

/**
 * Seed the deterministic PRNG.
 *
 * Call before a run and before replaying it. Upstream planck uses
 * `Math.random`, which cannot be reproduced and so cannot appear anywhere in a
 * simulation whose replay must be verifiable.
 */
export function seedRandom(seed: number): void {
  randomState = seed | 0;
}

/** Deterministic replacement for `Math.random` (mulberry32). */
export function random(): number {
  randomState = (randomState + 0x6d2b79f5) | 0;
  let t = randomState;
  t = Math.imul(t ^ (t >>> 15), t | 1);
  t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
  return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
}
