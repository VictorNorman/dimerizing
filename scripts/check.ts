/**
 * Headless sanity checks and a benchmark for the gas.
 *
 *   npm run check
 *
 * Verifies the invariants the model is supposed to hold -- conserved mass,
 * particles inside the box, disks not left interpenetrating, momentum and
 * energy conserved when no reaction happens, and Kc reaching the same
 * equilibrium from opposite starting compositions -- then reports throughput.
 */

import { Sim, SPECIES_A, type Params } from '../src/sim.ts';

const BASE: Params = {
  initialA: 200,
  initialB: 0,
  particleMass: 2,
  temperature: 8,
  dimerizationChance: 40,
  dissociationChance: 4,
  boxSize: 100,
};

let failures = 0;

function check(name: string, ok: boolean, detail = ''): void {
  if (ok) {
    console.log(`  ok    ${name}${detail ? '  ' + detail : ''}`);
  } else {
    failures++;
    console.log(`  FAIL  ${name}${detail ? '  ' + detail : ''}`);
  }
}

/** Reaches into the sim's private storage; this is a test, not a client. */
function peek(sim: Sim) {
  return sim as unknown as {
    x: Float64Array;
    y: Float64Array;
    vx: Float64Array;
    vy: Float64Array;
    mass: Float64Array;
    radius: Float64Array;
    species: Uint8Array;
    alive: Uint8Array;
    top: number;
    boxEdge: number;
  };
}

// ---------------------------------------------------------------------
console.log('\ninvariants over 20,000 steps (200 A, reactions on)');
{
  const sim = new Sim({ ...BASE });
  const expectedMass = BASE.initialA + 2 * BASE.initialB;

  let massOk = true;
  let insideOk = true;
  let sawB = false;
  let sawDissociation = false;
  let peakB = 0;

  for (let s = 0; s < 20000; s++) {
    sim.step();
    const st = sim.stats();
    if (st.totalAEquivalent !== expectedMass) {
      massOk = false;
    }
    if (st.countB > 0) {
      sawB = true;
    }
    if (st.countB > peakB) {
      peakB = st.countB;
    }
    if (sawB && st.countB < peakB) {
      sawDissociation = true;
    }

    if (s % 200 === 0) {
      const p = peek(sim);
      const inner = p.boxEdge - 0.5;
      for (let i = 0; i < p.top; i++) {
        if (!p.alive[i]) {
          continue;
        }
        const reach = inner - p.radius[i];
        // A hair of slack: positions are clamped exactly to `reach`.
        if (Math.abs(p.x[i]) > reach + 1e-9 || Math.abs(p.y[i]) > reach + 1e-9) {
          insideOk = false;
        }
      }
    }
  }

  const st = sim.stats();
  check('A + 2B conserved exactly', massOk, `= ${st.totalAEquivalent} (expected ${expectedMass})`);
  check('every particle inside the box', insideOk);
  check('dimerization happens', sawB, `peak B = ${peakB}`);
  check('dissociation happens', sawDissociation);
  check('temperature stays finite', Number.isFinite(st.measuredTemperature) && st.measuredTemperature > 0,
    `T = ${st.measuredTemperature.toFixed(3)}`);
}

// ---------------------------------------------------------------------
console.log('\ndisk overlap audit (every pair, every step, 12,000 steps)');
{
  // Sampling positions at discrete times always leaves a little residual
  // overlap: a pair can enter contact between two samples, and a reaction
  // replaces a collision rather than following it (exactly as the Python
  // model specifies), so a contact that makes a B dissociate leaves the
  // particle that hit it sitting among the fragments with nothing bounced.
  //
  // Depth is therefore bounded by how far a pair can close in one step,
  // 2 * MAX_STEP_DISPLACEMENT = 0.4 patches, comfortably under the 0.5 that
  // would mean two A's had swapped sides.
  //
  // How long an overlap lasts is not by itself interesting -- two slow
  // particles that start overlapped and drift apart at 0.1 patches/tick
  // legitimately take twenty steps to separate. What would be a real defect
  // is an overlap that *stops resolving*: a jammed pair whose disks stay
  // interpenetrated or sink deeper because a collision is being missed. So
  // that is what this asserts.
  const sim = new Sim({ ...BASE });
  const p = peek(sim);

  interface Track {
    first: number;
    last: number;
    age: number;
  }
  const live = new Map<number, Track>();
  let events = 0;
  let deepEvents = 0;
  let particleSteps = 0;
  let worst = 0;
  let longest = 0;
  let stuck = 0;

  /** A pair still no further apart than it started, long after contact. */
  const retire = (t: Track) => {
    if (t.age > longest) {
      longest = t.age;
    }
    if (t.age >= 5 && t.last >= t.first) {
      stuck++;
    }
  };

  for (let s = 0; s < 12000; s++) {
    sim.step();
    particleSteps += sim.count;
    const seen = new Set<number>();
    for (let i = 0; i < p.top; i++) {
      if (!p.alive[i]) {
        continue;
      }
      for (let j = i + 1; j < p.top; j++) {
        if (!p.alive[j]) {
          continue;
        }
        const dx = p.x[i] - p.x[j];
        const dy = p.y[i] - p.y[j];
        const contact = p.radius[i] + p.radius[j];
        const overlap = contact - Math.hypot(dx, dy);
        if (overlap <= 0.01) {
          continue;
        }

        events++;
        if (overlap > 0.15) {
          deepEvents++;
        }
        if (overlap > worst) {
          worst = overlap;
        }

        const key = i * 100000 + j;
        seen.add(key);
        const t = live.get(key);
        if (t === undefined) {
          live.set(key, {
            first: overlap,
            last: overlap,
            age: 1,
          });
        } else {
          t.last = overlap;
          t.age++;
        }
      }
    }
    for (const [key, t] of live) {
      if (seen.has(key)) {
        continue;
      }
      retire(t);
      live.delete(key);
    }
  }
  for (const t of live.values()) {
    retire(t);
  }

  const rate = events / particleSteps;
  check('particles never pass through each other', worst < 0.5,
    `worst ${worst.toFixed(3)} patches (A+A contact distance is 0.5)`);
  check('no overlap ever stops resolving', stuck === 0,
    `${stuck} jammed pairs; longest-lived overlap ${longest} steps`);
  check('overlap is rare', rate < 1e-4, `${(rate * 1e6).toFixed(1)} per million particle-steps`);
  check('deep overlap stays occasional', deepEvents / particleSteps < 3e-5,
    `${((deepEvents / particleSteps) * 1e6).toFixed(2)} per million deeper than 0.15`);
}

// ---------------------------------------------------------------------
console.log('\npure elastic gas (reactions off): momentum and energy');
{
  const sim = new Sim({
    ...BASE,
    initialA: 300,
    initialB: 50,
    dimerizationChance: 0,
    dissociationChance: 0,
  });

  const totals = () => {
    const p = peek(sim);
    let px = 0;
    let py = 0;
    let ke = 0;
    for (let i = 0; i < p.top; i++) {
      if (!p.alive[i]) {
        continue;
      }
      px += p.mass[i] * p.vx[i];
      py += p.mass[i] * p.vy[i];
      ke += 0.5 * p.mass[i] * (p.vx[i] * p.vx[i] + p.vy[i] * p.vy[i]);
    }
    return {
      px,
      py,
      ke,
    };
  };

  const before = totals();
  for (let s = 0; s < 20000; s++) {
    sim.step();
  }
  const after = totals();

  const keDrift = Math.abs(after.ke - before.ke) / before.ke;
  check('kinetic energy conserved', keDrift < 1e-9, `drift ${(keDrift * 100).toExponential(2)} %`);
  // Momentum is NOT conserved overall -- the walls are immovable and absorb
  // it -- so this only checks that nothing exploded.
  check('momentum stays bounded', Number.isFinite(after.px) && Number.isFinite(after.py));

  const st = sim.stats();
  check('no species change without reactions', st.countA === 300 && st.countB === 50,
    `A = ${st.countA}, B = ${st.countB}`);
}

// ---------------------------------------------------------------------
console.log('\nMaxwell-Boltzmann: speeds thermalize toward the set temperature');
{
  const sim = new Sim({
    ...BASE,
    initialA: 600,
    dimerizationChance: 0,
    dissociationChance: 0,
    temperature: 20,
  });
  for (let s = 0; s < 40000; s++) {
    sim.step();
  } // let the uniform initial speeds thermalize
  const st = sim.stats();
  // Every particle starts at exactly the thermal speed; collisions spread
  // them into a distribution, but the mean energy must not drift.
  const err = Math.abs(st.measuredTemperature - 20) / 20;
  check('mean KE per particle stays at the set temperature', err < 0.02,
    `measured ${st.measuredTemperature.toFixed(3)} vs 20`);

  const p = peek(sim);
  const spreadRatio = () => {
    let n = 0;
    let mean = 0;
    let m2 = 0;
    for (let i = 0; i < p.top; i++) {
      if (!p.alive[i]) {
        continue;
      }
      const v = Math.hypot(p.vx[i], p.vy[i]);
      n++;
      const d = v - mean;
      mean += d / n;
      m2 += d * (v - mean);
    }
    return Math.sqrt(m2 / n) / mean;
  };

  // 600 particles is few enough that a single snapshot's sd/mean carries
  // about 1.5 % sampling noise, so average over many well-separated ones.
  let sum = 0;
  let samples = 0;
  for (let s = 0; s < 60000; s++) {
    sim.step();
    if (s % 500 === 0) {
      sum += spreadRatio();
      samples++;
    }
  }
  const ratio = sum / samples;
  // For a 2D Maxwell-Boltzmann distribution, sd/mean = sqrt(4/pi - 1) ~= 0.5227.
  check('speed spread matches 2D Maxwell-Boltzmann', Math.abs(ratio - 0.5227) < 0.02,
    `sd/mean = ${ratio.toFixed(3)} (expected 0.523)`);
}

// ---------------------------------------------------------------------
console.log('\ndynamic equilibrium: same Kc from opposite starting compositions');
{
  const settle = (initialA: number, initialB: number) => {
    const sim = new Sim({
      ...BASE,
      initialA,
      initialB,
    });
    for (let s = 0; s < 60000; s++) {
      sim.step();
    }
    // Average Kc over a long tail rather than trusting one noisy sample.
    let sum = 0;
    let n = 0;
    for (let s = 0; s < 120000; s++) {
      sim.step();
      if (s % 200 === 0) {
        sum += sim.stats().kc;
        n++;
      }
    }
    return sum / n;
  };

  // Both runs hold A + 2B = 400, approached from all-A and from all-B.
  const fromA = settle(400, 0);
  const fromB = settle(0, 200);
  const spread = Math.abs(fromA - fromB) / ((fromA + fromB) / 2);
  check('Kc agrees from both directions', spread < 0.15,
    `from A: ${fromA.toFixed(1)}, from B: ${fromB.toFixed(1)} (${(spread * 100).toFixed(1)} % apart)`);
}

// ---------------------------------------------------------------------
console.log('\nthroughput  (the Python original manages ~274 sim-ticks/s at 200 A, ~51 at 1000 A)');
for (const n of [200, 1000, 4000]) {
  const sim = new Sim({
    ...BASE,
    initialA: n,
  });
  for (let s = 0; s < 200; s++) {
    sim.step();
  } // warm the JIT

  const steps = 4000;
  const ticksBefore = sim.ticks;
  const t0 = performance.now();
  for (let s = 0; s < steps; s++) {
    sim.step();
  }
  const seconds = (performance.now() - t0) / 1000;

  // sim-ticks is the honest cross-implementation number: this version takes
  // shorter steps than the Python original (see MAX_STEP_DISPLACEMENT), so
  // raw steps/s would flatter it.
  const perStep = (seconds * 1e6) / steps;
  console.log(
    `  ${String(n).padStart(4)} A: ${(steps / seconds).toFixed(0).padStart(7)} steps/s` +
      `  ${((sim.ticks - ticksBefore) / seconds).toFixed(0).padStart(6)} sim-ticks/s` +
      `  (${perStep.toFixed(1)} us/step, ${((perStep * 1000) / n).toFixed(1)} ns/particle/step)`,
  );
}

console.log(failures === 0 ? '\nall checks passed\n' : `\n${failures} check(s) FAILED\n`);
process.exit(failures === 0 ? 0 : 1);
