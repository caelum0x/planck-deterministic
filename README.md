# planck-deterministic

A fork of [planck.js](https://github.com/piqnt/planck.js) (MIT, Erin Catto / Ali Shakiba)
with a **cross-platform deterministic math kernel**.

## Why

Upstream planck calls the host `Math.sin`, `Math.cos`, `Math.atan2` and `Math.random`.
IEEE-754 guarantees `+ - * /` and `sqrt` are correctly rounded, so those are identical
everywhere — but ECMAScript explicitly permits an *implementation-dependent approximation*
for the transcendental functions, and every engine forwards them to the host libm. Those
results differ between ARM64/JavaScriptCore (iOS) and x86-64/V8 (desktop).

`Rot.setAngle` and `Sweep.getTransform` call sin/cos for **every rotating body on every
step**. So a simulation shared between a mobile build and a desktop build diverges as soon
as anything rotates. That breaks replay portability, server-side replay verification, and
any leaderboard built on re-simulating a replay.

## What changed

- Added `src/common/DeterministicMath.ts` — `sin`, `cos`, `exp`, `atan`, `atan2`, `sqrt` and
  a seeded PRNG, built from `+ - * /` and comparisons only. Coefficients and the Cody-Waite
  range reduction follow fdlibm (Sun, freely distributable), the reference implementation
  most libms derive from. Accuracy is within a few ulp of the host libm.
- Rewired the nine hoisted aliases in `Math.ts`, `Matrix.ts`, `Rot.ts`, `Sweep.ts` and
  `dynamics/Position.ts` to use it. **No host libm call remains anywhere in `src/`.**
- Exported the kernel (`DeterministicMath`, `seedRandom`) so consumers can share the exact
  same functions the solver uses and seed the PRNG for reproducible replays.

Every upstream test still passes, plus new determinism tests including a check that a
simulated free fall matches the closed-form analytic solution.

```js
import { World, seedRandom, DeterministicMath } from "planck-deterministic";

seedRandom(12345); // before a run, and before replaying it
```

## Upstream

This tracks planck.js v1.5.0. Upstream docs and API below.

---

# Planck.js

Planck.js is JavaScript/TypeScript rewrite of Box2D physics engine for cross-platform HTML5 game development.

#### Motivations

- Taking advantage of Box2D's efforts and achievements
- Developing readable and editable code in JavaScript/TypeScript
- Providing idiomatic JavaScript/TypeScript API
- Optimizing the library for web and mobile platforms

#### [Documentation](https://piqnt.com/planck.js/docs/)

#### [Examples](https://piqnt.com/planck.js/)

#### [Discord](https://discord.com/invite/znjh6J7)

#### [Made with Planck.js](https://github.com/piqnt/planck.js/wiki/)

#### [Report Issues](https://github.com/piqnt/planck.js/issues)
To speed up resolving issues, please provide [testbed](https://piqnt.com/planck.js/docs/testbed) code to reproduce the issue.
