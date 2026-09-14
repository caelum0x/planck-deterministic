import { describe, it, expect } from "vitest";

import * as DM from "../DeterministicMath";
import { World, Box, Edge } from "../../index";

/**
 * The reason this fork exists.
 *
 * Upstream planck calls the host libm for sin/cos/atan2. ECMAScript permits
 * those to be implementation-dependent approximations, and they genuinely
 * differ between ARM64/JavaScriptCore and x86-64/V8. Since `Rot.setAngle` and
 * `Sweep.getTransform` call them for every rotating body on every step, a
 * simulation shared between a mobile and a desktop build drifts apart as soon
 * as anything rotates — which breaks replay portability and any leaderboard
 * that depends on re-simulating a replay.
 */
describe("deterministic math kernel", () => {
  it("tracks the host libm closely", () => {
    for (let i = -720; i <= 720; i++) {
      const x = (i * Math.PI) / 180;
      expect(Math.abs(DM.sin(x) - Math.sin(x))).toBeLessThan(1e-14);
      expect(Math.abs(DM.cos(x) - Math.cos(x))).toBeLessThan(1e-14);
    }
    for (let i = 0; i <= 200; i++) {
      const x = -i * 0.25;
      expect(Math.abs(DM.exp(x) - Math.exp(x))).toBeLessThan(1e-14);
    }
    for (const x of [0, 0.4375, 0.6875, 1, 1.1875, 2.4375, 100, -5]) {
      expect(Math.abs(DM.atan(x) - Math.atan(x))).toBeLessThan(1e-14);
    }
  });

  it("stays accurate for large accumulated angles", () => {
    for (const x of [1e3, 1e4, 12345.6789, -98765.4321]) {
      expect(Math.abs(DM.sin(x) - Math.sin(x))).toBeLessThan(1e-10);
      expect(Math.abs(DM.cos(x) - Math.cos(x))).toBeLessThan(1e-10);
    }
  });

  it("holds the Pythagorean identity", () => {
    for (let i = 0; i < 500; i++) {
      const x = -50 + i * 0.2;
      expect(Math.abs(DM.sin(x) ** 2 + DM.cos(x) ** 2 - 1)).toBeLessThan(1e-14);
    }
  });

  it("does NOT merely delegate to the host libm", () => {
    // If this ever starts failing it means someone replaced the kernel with
    // Math.* and the fork no longer does anything.
    let ours = 0;
    let host = 0;
    for (let i = 1; i <= 2000; i++) {
      const x = i * 0.037;
      ours += DM.sin(x) + DM.cos(x) + DM.atan(x);
      host += Math.sin(x) + Math.cos(x) + Math.atan(x);
    }
    expect(Math.abs(ours - host)).toBeLessThan(1e-9); // still accurate
    expect(ours).not.toBe(host); // but independently computed
  });

  it("seeds a reproducible PRNG", () => {
    DM.seedRandom(99);
    const first = [DM.random(), DM.random(), DM.random()];
    DM.seedRandom(99);
    expect([DM.random(), DM.random(), DM.random()]).toEqual(first);
    DM.seedRandom(100);
    expect(DM.random()).not.toBe(first[0]);
  });
});

describe("deterministic simulation", () => {
  /** A falling, rotating, bouncing box — exercises the patched call sites. */
  function trajectory(steps: number): number[] {
    DM.seedRandom(12345);
    const world = new World({ x: 0, y: -1.625 }); // lunar gravity

    const ground = world.createBody();
    ground.createFixture(new Edge({ x: -100, y: 0 }, { x: 100, y: 0 }), {
      friction: 0.8,
    });

    const body = world.createDynamicBody({
      position: { x: 0, y: 40 },
      angle: 0.3,
      angularVelocity: 1.7,
      linearVelocity: { x: 2.5, y: 0 },
    });
    body.createFixture(new Box(1.5, 1.5), {
      density: 1.2,
      friction: 0.7,
      restitution: 0.1,
    });

    const out: number[] = [];
    for (let i = 0; i < steps; i++) {
      world.step(1 / 120);
      const p = body.getPosition();
      out.push(p.x, p.y, body.getAngle());
    }
    return out;
  }

  it("reproduces an identical trajectory across runs", () => {
    expect(trajectory(600)).toEqual(trajectory(600));
  });

  it("matches the analytic free-fall solution", () => {
    // 5 s of lunar gravity from 40 m. Semi-implicit Euler adds gravity to
    // velocity before position, overshooting the closed form by g*dt*t/2.
    const g = 1.625;
    const dt = 1 / 120;
    const t = 5;
    const analytic = 40 - 0.5 * g * t * t - 0.5 * g * dt * t;
    const samples = trajectory(600);
    const finalY = samples[samples.length - 2] as number;
    expect(Math.abs(finalY - analytic)).toBeLessThan(0.01);
  });

  it("rotates, which is what makes the patched call sites matter", () => {
    const samples = trajectory(600);
    const firstAngle = samples[2] as number;
    const lastAngle = samples[samples.length - 1] as number;
    expect(Math.abs(lastAngle - firstAngle)).toBeGreaterThan(1);
  });
});
