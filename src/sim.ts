/**
 * Dimerizing Gas: a 2D hard-disk gas with the reversible reaction 2A <=> B.
 *
 * This is a port of models/dimerizing_gas.py from the netlogo-py project,
 * rewritten for speed. The model -- box, species, masses, radii, thermal
 * speeds, per-collision reaction probabilities, momentum-conserving
 * combine(), thermal-kick dissociate() -- is the same. The mechanics
 * underneath are not:
 *
 *   - Velocities are stored as (vx, vy) components instead of NetLogo's
 *     (speed, heading) polar pair. The Python version paid four trig calls
 *     per elastic collision plus a sin/cos per wall test; this version pays
 *     none. `speed` is derived with a sqrt only where it is actually needed
 *     (the histogram, the measured temperature, the timestep bound).
 *
 *   - Particles live in parallel typed arrays (a struct-of-arrays layout),
 *     not objects, so a whole pass over the population is a linear walk
 *     over contiguous memory that the JIT keeps unboxed.
 *
 *   - The neighbor search is a uniform grid rebuilt each step by counting
 *     sort into flat Int32Arrays -- no hash map, no per-cell array objects,
 *     zero allocation per step. The Python version rebuilt a dict of lists
 *     every tick, which allocated ~n small objects per tick.
 *
 *   - Collisions are resolved the way falstad.com/gas does it: when two
 *     disks are found overlapping, solve the quadratic for the exact
 *     (negative) time of contact, rewind both particles to the instant they
 *     touched, apply the elastic impulse along the true line of centers,
 *     then replay the rest of the step. The Python model instead applied
 *     the impulse along a *random* axis (inherited from NetLogo's GasLab)
 *     and needed a `last_collision` bookkeeping field to stop overlapping
 *     pairs from re-colliding forever. Exact contact resolution leaves no
 *     residual overlap, so that whole mechanism is gone.
 *
 * Everything the model is actually about -- Maxwell-Boltzmann speeds,
 * dynamic equilibrium, a steady Kc = B/A^2 -- is unchanged.
 */

export const SPECIES_A = 0;
export const SPECIES_B = 1;

/** World half-width, matching `resize_world(-40, 40, -40, 40)`. */
export const MAX_PXCOR = 40;
/**
 * Particle radii in patches. These are physical, not decorative: two
 * particles collide when their disks genuinely overlap.
 */
export const A_RADIUS = 0.25;
/**
 * A B is one A's worth of matter twice over, so it takes up twice the space:
 * cbrt(2) ~ 1.26 times the radius, the ratio that makes the two *volumes*
 * 1:2. (The disks on screen are a 2D window onto 3D molecules -- scaling the
 * drawn areas 1:2 instead would need sqrt(2), and would mean a dimer whose
 * density differs from its monomers'.)
 */
export const B_RADIUS = A_RADIUS * Math.cbrt(2);

/** The widest two centers can be and still touch (B + B). */
const MAX_CONTACT_DISTANCE = 2 * B_RADIUS;

/** Upper bound on the timestep, from the Python model. */
const MAX_TICK_DELTA = 0.1073;
/**
 * Upper bound on how far a particle may travel in one step, in patches.
 * Contact is detected by testing for overlap at discrete positions, so the
 * step has to be short enough that an approaching pair cannot pass through
 * each other between two samples. The tightest case is A+A, whose disks
 * only overlap within 0.5 patches of center separation and which can close
 * at twice this speed, so 0.2 leaves better than a 2x margin. The Python
 * model used 1/ceiling(max speed) here -- up to a full patch per step,
 * which lets fast A pairs tunnel straight through one another.
 */
const MAX_STEP_DISPLACEMENT = 0.2;

/**
 * How many random axes a splitting B may try before the reaction is declined
 * for lack of room. A B wedged against one neighbor usually still has a free
 * line to come apart along, so a single draw would refuse far more splits
 * than the geometry actually forbids.
 */
const SPLIT_ORIENTATION_TRIES = 6;

// =======================================================================
// Absolute units
//
// Fixing one length turns the reduced-unit model into a dimensional one.
// Read the temperature slider as real kelvin and the mass slider as real
// amu, and every other unit follows -- including the tick, which is not a
// free choice once the other three are made.
// =======================================================================

/** The box interior, wall face to wall face, at 100% operating area. */
export const BOX_SIDE_NM = 100;

/** Boltzmann's constant, J/K. */
const K_B = 1.380649e-23;
/** One atomic mass unit, kg. */
const AMU_KG = 1.66053906660e-27;

/**
 * Patch length in metres. The wall centerline sits at MAX_PXCOR and its
 * inner face half a patch inside that, so the interior spans
 * 2 * (MAX_PXCOR - 0.5) patches when the box is at full size.
 */
export const METRES_PER_PATCH = (BOX_SIDE_NM * 1e-9) / (2 * (MAX_PXCOR - 0.5));

/**
 * One patch per tick in m/s. A reduced speed is sqrt(2T/m) on the sliders'
 * bare numbers; the same speed in SI is sqrt(2 k_B T / (m * amu)). The ratio
 * is sqrt(k_B / amu) -- which depends on neither slider, so one tick has a
 * single fixed duration no matter how the gas is set up.
 */
const METRES_PER_SECOND_PER_PATCH_TICK = Math.sqrt(K_B / AMU_KG);

/** Seconds in one tick: about 13.9 ps for a 100 nm box. */
export const SECONDS_PER_TICK = METRES_PER_PATCH / METRES_PER_SECOND_PER_PATCH_TICK;

/**
 * Reduced pressure to N/m. A 2D pressure is force per unit *length*, so in
 * reduced units it comes out as mass/time^2 -- amu per tick squared.
 */
export const PRESSURE_TO_N_PER_M = AMU_KG / (SECONDS_PER_TICK * SECONDS_PER_TICK);

/** Patches to nanometres. */
export const NM_PER_PATCH = METRES_PER_PATCH * 1e9;

export interface Params {
  /** Number of A particles created by setup. */
  initialA: number;
  /** Number of B particles created by setup. */
  initialB: number;
  /** Mass of one A particle, in amu (B is twice this). */
  particleMass: number;
  /** Temperature in reduced units: average kinetic energy per particle. */
  temperature: number;
  /** Chance (%) that an A+A contact fuses instead of bouncing. */
  dimerizationChance: number;
  /** Chance (%) that a B falls apart when it touches anything, wall included. */
  dissociationChance: number;
  /**
   * Box size as a percentage of the world's *area*, not of its half-width.
   * Read only by setup. Because the box is two-dimensional, a given
   * percentage of its area shrinks each side by sqrt(percent/100), not by
   * percent/100 -- see setup() and boxEdge's sqrt().
   */
  boxSize: number;
}

export interface Stats {
  countA: number;
  countB: number;
  /** count_a + 2 * count_b -- conserved exactly, whatever else happens. */
  totalAEquivalent: number;
  /**
   * Mean kinetic energy per particle. With k_B = 1 and two translational
   * degrees of freedom this *is* the temperature in kelvin, and since no
   * reaction creates or destroys energy any more, it is a result of the run
   * rather than a setting -- the slider only fixes where it starts.
   */
  measuredTemperature: number;
  /** [B] / [A]^2 in nm^2, from real 2D concentrations. */
  kc: number;
  /** Interior area in nm^2. */
  areaNm2: number;
  /** Concentrations, particles per nm^2. */
  concA: number;
  concB: number;
  /** Total kinetic energy in reduced units; conserved, so a leak shows here. */
  kineticEnergy: number;
  /** Net drift the gas has picked up, in patches/tick. */
  driftSpeed: number;
  /** N k_B T / A in N/m -- what an ideal 2D gas would push with. */
  idealPressure: number;
}

/**
 * The gas. One instance owns all particle storage and is reused across
 * setups: `setup()` reallocates only when the population outgrows what is
 * already there.
 */
export class Sim {
  params: Params;

  // ---- Particle storage (struct of arrays, indexed by slot) -------------
  private capacity = 0;
  private x = new Float64Array(0);
  private y = new Float64Array(0);
  private vx = new Float64Array(0);
  private vy = new Float64Array(0);
  private mass = new Float64Array(0);
  private radius = new Float64Array(0);
  private species = new Uint8Array(0);
  private alive = new Uint8Array(0);
  /**
   * The step during which each slot's particle was created. A particle born
   * mid-step is skipped by that step's collision sweep, so a fresh
   * dissociation cannot re-fuse before it has moved -- the job the Python
   * model gave to `last_collision`.
   */
  private bornAt = new Int32Array(0);

  /**
   * Slots are stable: killing a particle marks it dead and pushes its slot
   * onto this free stack rather than compacting the arrays, so the grid's
   * indices stay valid for the rest of the step.
   */
  private freeSlots = new Int32Array(0);
  private freeCount = 0;
  /** Slots [0, top) have been handed out at least once. */
  private top = 0;
  /** Living particles. */
  count = 0;

  // ---- Uniform grid ----------------------------------------------------
  private gridN = 0;
  private gridOrigin = 0;
  /**
   * Cell width. Never below the widest contact distance -- that is what makes
   * the 3x3 neighborhood scan complete -- but widened toward one particle per
   * cell when the gas is sparse. Rebuilding the grid costs one pass over the
   * cells as well as one over the particles, so a 200-particle run in a
   * fixed 1x1 grid would spend most of its time clearing 5,625 empty cells.
   */
  private cellSize = MAX_CONTACT_DISTANCE;
  private cellStart = new Int32Array(0);
  private cellCursor = new Int32Array(0);
  private cellItems = new Int32Array(0);
  private cellOf = new Int32Array(0);

  // ---- Clock -----------------------------------------------------------
  ticks = 0;
  private stepCounter = 0;
  private tickDelta = MAX_TICK_DELTA;
  /** Fastest particle seen during the last move pass; sets the next step. */
  private maxSpeedSq = 0;

  // ---- Wall-impulse pressure gauge --------------------------------------
  /**
   * Momentum delivered to the walls since the last drain, per species, in
   * reduced units (amu * patch / tick). A 2D pressure is force per unit
   * length, so dividing this by the perimeter and by the elapsed time gives
   * a partial pressure directly -- pressure as molecular impacts, with no
   * equation of state assumed anywhere.
   */
  private wallImpulse = [0, 0];
  /** Simulated time the accumulated impulse covers. */
  private wallImpulseTicks = 0;

  // ---- Geometry --------------------------------------------------------
  /** Wall centerline, in patches from the origin. */
  boxEdge = 0;
  /** Inner face of the wall: particle centers stay within boxEdge - 0.5. */
  private innerBound = 0;

  constructor(params: Params) {
    this.params = params;
    this.setup();
  }

  // =====================================================================
  // Setup
  // =====================================================================

  setup(): void {
    const p = this.params;
    // sqrt because boxSize is a fraction of the world's *area*, not of its
    // half-width -- the box is two-dimensional, so halving its area shrinks
    // each side by sqrt(2), not by 2, exactly as the Python model does.
    this.boxEdge = Math.round(MAX_PXCOR * Math.sqrt(p.boxSize / 100));
    this.innerBound = this.boxEdge - 0.5;

    // A + 2B is conserved, so the population can never exceed its
    // all-A value -- that is exactly how many slots this run can need.
    const needed = p.initialA + 2 * p.initialB + 4;
    if (needed > this.capacity) {
      this.allocate(needed);
    }

    this.alive.fill(0);
    this.freeCount = 0;
    this.top = 0;
    this.count = 0;
    this.ticks = 0;
    this.stepCounter = 0;
    this.maxSpeedSq = 0;

    // A + 2B is conserved, so the all-A population is this run's density
    // ceiling and sizing the grid to it once is enough.
    this.buildGridIndex(p.initialA + 2 * p.initialB);

    for (let i = 0; i < p.initialA; i++) {
      this.placeRandomly(this.spawn(SPECIES_A, p.particleMass, 0, 0, 0, 0));
    }
    for (let i = 0; i < p.initialB; i++) {
      this.placeRandomly(this.spawn(SPECIES_B, p.particleMass * 2, 0, 0, 0, 0));
    }
    this.tickDelta = this.computeTickDelta();
  }

  /** Reallocates every array to `n` slots, preserving what is already there. */
  private grow(n: number): void {
    const old = {
      x: this.x, y: this.y, vx: this.vx, vy: this.vy,
      mass: this.mass, radius: this.radius,
      species: this.species, alive: this.alive, bornAt: this.bornAt,
      freeSlots: this.freeSlots,
    };
    this.allocate(n);
    this.x.set(old.x);
    this.y.set(old.y);
    this.vx.set(old.vx);
    this.vy.set(old.vy);
    this.mass.set(old.mass);
    this.radius.set(old.radius);
    this.species.set(old.species);
    this.alive.set(old.alive);
    this.bornAt.set(old.bornAt);
    this.freeSlots.set(old.freeSlots);
  }

  private allocate(n: number): void {
    this.capacity = n;
    this.x = new Float64Array(n);
    this.y = new Float64Array(n);
    this.vx = new Float64Array(n);
    this.vy = new Float64Array(n);
    this.mass = new Float64Array(n);
    this.radius = new Float64Array(n);
    this.species = new Uint8Array(n);
    this.alive = new Uint8Array(n);
    this.bornAt = new Int32Array(n);
    this.freeSlots = new Int32Array(n);
    this.cellItems = new Int32Array(n);
    this.cellOf = new Int32Array(n);
  }

  /** A slot for a new particle, with a thermal speed in a random direction. */
  private spawn(species: number, mass: number, px: number, py: number, vx: number, vy: number): number {
    let i: number;
    if (this.freeCount > 0) {
      i = this.freeSlots[--this.freeCount];
    } else {
      // A + 2B is conserved, so setup() already sized this correctly and
      // this branch should be unreachable. Grow anyway rather than write
      // past the end if some future change breaks that invariant.
      if (this.top >= this.capacity) {
        this.grow(Math.max(8, this.capacity * 2));
      }
      i = this.top++;
    }
    this.alive[i] = 1;
    this.bornAt[i] = this.stepCounter;
    this.species[i] = species;
    this.mass[i] = mass;
    this.radius[i] = species === SPECIES_A ? A_RADIUS : B_RADIUS;
    this.x[i] = px;
    this.y[i] = py;
    this.vx[i] = vx;
    this.vy[i] = vy;
    this.count++;
    return i;
  }

  private kill(i: number): void {
    this.alive[i] = 0;
    this.freeSlots[this.freeCount++] = i;
    this.count--;
  }

  private placeRandomly(i: number): void {
    const reach = this.innerBound - this.radius[i];
    this.x[i] = (Math.random() * 2 - 1) * reach;
    this.y[i] = (Math.random() * 2 - 1) * reach;
    const speed = this.thermalSpeed(this.mass[i]);
    const theta = Math.random() * 2 * Math.PI;
    this.vx[i] = speed * Math.cos(theta);
    this.vy[i] = speed * Math.sin(theta);
  }

  /**
   * The speed at which a particle of this mass carries kinetic energy equal
   * to `temperature` -- equipartition with k_B = 1, which is what the
   * temperature slider means in this model.
   */
  private thermalSpeed(mass: number): number {
    return Math.sqrt((2 * this.params.temperature) / mass);
  }

  // =====================================================================
  // Stepping
  // =====================================================================

  /** Advances the gas by one timestep and reports how much simulated time passed. */
  step(): number {
    const dt = this.tickDelta;
    this.stepCounter++;
    this.moveAndBounce(dt);
    this.buildGrid();
    this.resolveCollisions();
    this.ticks += dt;
    this.tickDelta = this.computeTickDelta();
    return dt;
  }

  /**
   * The next step length: short enough that nothing outruns the contact
   * test, and never longer than the model's own ceiling.
   */
  private computeTickDelta(): number {
    if (this.maxSpeedSq <= 0) {
      // Nothing has moved yet -- fall back to the population's thermal speed
      // so the very first step is not taken at the maximum length.
      let fastest = 0;
      for (let i = 0; i < this.top; i++) {
        if (!this.alive[i]) {
          continue;
        }
        const s = this.vx[i] * this.vx[i] + this.vy[i] * this.vy[i];
        if (s > fastest) {
          fastest = s;
        }
      }
      this.maxSpeedSq = fastest;
    }
    if (this.maxSpeedSq <= 0) {
      return MAX_TICK_DELTA;
    }
    return Math.min(MAX_TICK_DELTA, MAX_STEP_DISPLACEMENT / Math.sqrt(this.maxSpeedSq));
  }

  /**
   * Advances every particle and reflects it off the box walls. A B particle
   * that reaches a wall may fall apart there instead of simply bouncing.
   *
   * Unlike the Python version, which reflected *before* moving based on a
   * predicted position, this moves first and then clamps -- so a particle
   * can never end a step buried in (or beyond) a wall, no matter how the
   * step length worked out.
   */
  private moveAndBounce(dt: number): void {
    const { x, y, vx, vy, alive, radius, species, mass, bornAt } = this;
    const dissociationChance = this.params.dissociationChance;
    const inner = this.innerBound;
    const impulse = this.wallImpulse;
    let maxSpeedSq = 0;
    const top = this.top;
    this.wallImpulseTicks += dt;

    for (let i = 0; i < top; i++) {
      if (!alive[i]) {
        continue;
      }
      // A wall dissociation earlier in this same pass can drop a fragment
      // into a recycled slot ahead of the cursor; it moves next step.
      if (bornAt[i] === this.stepCounter) {
        continue;
      }

      let px = x[i] + vx[i] * dt;
      let py = y[i] + vy[i] * dt;
      // Reach is measured from the disk's edge, so a particle never
      // visibly buries itself in the yellow wall before turning around.
      const reach = inner - radius[i];
      let hitWall = false;

      // Mirror the overshoot back into the box rather than clamping to the
      // wall. Reflection is what a hard wall actually does, and clamping is
      // a projection: two particles that reach the same corner in one step
      // both land on the identical point, which leaves them exactly
      // coincident and takes several steps of collision response to undo.
      // Each reflection reverses one velocity component, so the wall takes
      // 2 * m * |v| of momentum from that axis. Summing those is the whole
      // pressure measurement.
      const sp = species[i];
      if (px > reach) {
        px = 2 * reach - px;
        impulse[sp] += 2 * mass[i] * Math.abs(vx[i]);
        vx[i] = -Math.abs(vx[i]);
        hitWall = true;
      } else if (px < -reach) {
        px = -2 * reach - px;
        impulse[sp] += 2 * mass[i] * Math.abs(vx[i]);
        vx[i] = Math.abs(vx[i]);
        hitWall = true;
      }
      if (py > reach) {
        py = 2 * reach - py;
        impulse[sp] += 2 * mass[i] * Math.abs(vy[i]);
        vy[i] = -Math.abs(vy[i]);
        hitWall = true;
      } else if (py < -reach) {
        py = -2 * reach - py;
        impulse[sp] += 2 * mass[i] * Math.abs(vy[i]);
        vy[i] = Math.abs(vy[i]);
        hitWall = true;
      }

      // A mirror can only leave the box if the overshoot exceeded the box
      // itself, which the timestep bound rules out -- but a particle that
      // starts outside (a shrunken box, a product born against a wall) would
      // stay outside without this.
      if (px > reach) {
        px = reach;
      } else if (px < -reach) {
        px = -reach;
      }
      if (py > reach) {
        py = reach;
      } else if (py < -reach) {
        py = -reach;
      }

      x[i] = px;
      y[i] = py;

      const s = vx[i] * vx[i] + vy[i] * vy[i];
      if (s > maxSpeedSq) {
        maxSpeedSq = s;
      }

      if (hitWall && species[i] === SPECIES_B && Math.random() * 100 < dissociationChance) {
        // Bounce first, then break apart: the fragments inherit a velocity
        // that already points away from the wall. A declined split (no room
        // for the fragments) just leaves the B bounced.
        //
        // The room check here reads a grid built before this pass moved
        // anything, so it is judging positions up to MAX_STEP_DISPLACEMENT
        // stale. That is a slightly blurred veto, not a wrong one -- the
        // exact test runs in resolveCollisions() once the grid is current.
        this.dissociate(i);
      }
    }
    this.maxSpeedSq = maxSpeedSq;
  }

  // =====================================================================
  // Uniform grid
  // =====================================================================

  /**
   * Sizes the grid to the current box and to the largest population this run
   * can reach. Dimerization only ever reduces the particle count from there,
   * so a one-time sizing stays valid (and slightly conservative) throughout.
   */
  private buildGridIndex(maxPopulation: number): void {
    const span = 2 * this.innerBound;
    // Two competing costs: sweeping the cells (falls as cellSize^2 grows)
    // and scanning each cell's 3x3 candidate block (rises with it). A little
    // under one particle per cell balances them; 0.6 is where a sweep of
    // 200 to 4,000 particles measured best. Below MAX_CONTACT_DISTANCE the
    // neighborhood scan would stop being complete, so that is the floor --
    // which is where a dense gas ends up anyway.
    const ideal = maxPopulation > 0 ? 0.6 * Math.sqrt((span * span) / maxPopulation) : span;
    this.cellSize = Math.max(MAX_CONTACT_DISTANCE, ideal);
    this.gridN = Math.max(1, Math.ceil(span / this.cellSize));
    this.gridOrigin = -this.innerBound;
    const cells = this.gridN * this.gridN;
    if (this.cellStart.length < cells + 1) {
      this.cellStart = new Int32Array(cells + 1);
      this.cellCursor = new Int32Array(cells);
    }
  }

  /**
   * Buckets every living particle by cell with a counting sort. Two passes
   * over the population and one over the cells, all into preallocated flat
   * arrays -- no hashing, no per-cell arrays, nothing allocated per step.
   */
  private buildGrid(): void {
    const { cellOf, cellItems, alive, x, y } = this;
    const n = this.gridN;
    const cells = n * n;
    const cellStart = this.cellStart;
    const cursor = this.cellCursor;
    const origin = this.gridOrigin;
    const invCell = 1 / this.cellSize;
    const top = this.top;

    cellStart.fill(0, 0, cells + 1);

    // Pass 1: count per cell, remembering each particle's cell.
    for (let i = 0; i < top; i++) {
      if (!alive[i]) {
        continue;
      }
      let cx = ((x[i] - origin) * invCell) | 0;
      let cy = ((y[i] - origin) * invCell) | 0;
      if (cx < 0) {
        cx = 0;
      } else if (cx >= n) {
        cx = n - 1;
      }
      if (cy < 0) {
        cy = 0;
      } else if (cy >= n) {
        cy = n - 1;
      }
      const c = cx + cy * n;
      cellOf[i] = c;
      cellStart[c + 1]++;
    }

    // Prefix sum turns the counts into slice starts.
    for (let c = 0; c < cells; c++) {
      cellStart[c + 1] += cellStart[c];
      cursor[c] = cellStart[c];
    }

    // Pass 2: scatter particles into their slices.
    for (let i = 0; i < top; i++) {
      if (!alive[i]) {
        continue;
      }
      cellItems[cursor[cellOf[i]]++] = i;
    }
  }

  // =====================================================================
  // Collisions and reactions
  // =====================================================================

  /**
   * Tests every pair that could possibly be touching and resolves it.
   *
   * Each cell is checked against itself and four of its eight neighbors --
   * the half of the neighborhood that has not already been visited from the
   * other side -- so every pair is considered exactly once without any
   * "only the higher index tests" bookkeeping.
   */
  private resolveCollisions(): void {
    const n = this.gridN;
    const cellStart = this.cellStart;
    const items = this.cellItems;

    for (let cy = 0; cy < n; cy++) {
      for (let cx = 0; cx < n; cx++) {
        const c = cx + cy * n;
        const startA = cellStart[c];
        const endA = cellStart[c + 1];
        if (startA === endA) {
          continue;
        }

        // Within the cell: every unordered pair, once.
        for (let a = startA; a < endA; a++) {
          const i = items[a];
          if (!this.alive[i]) {
            continue;
          }
          for (let b = a + 1; b < endA; b++) {
            const j = items[b];
            if (!this.alive[j]) {
              continue;
            }
            this.interact(i, j);
            // i may have just been consumed by a reaction; its slot can even
            // have been handed straight back out to a product particle.
            if (!this.alive[i] || this.bornAt[i] === this.stepCounter) {
              break;
            }
          }
        }

        // Forward half of the neighborhood: E, NW, N, NE.
        this.crossCell(c, cx + 1, cy, n);
        this.crossCell(c, cx - 1, cy + 1, n);
        this.crossCell(c, cx, cy + 1, n);
        this.crossCell(c, cx + 1, cy + 1, n);
      }
    }
  }

  private crossCell(c: number, ox: number, oy: number, n: number): void {
    if (ox < 0 || oy < 0 || ox >= n || oy >= n) {
      return;
    }
    const other = ox + oy * n;
    const cellStart = this.cellStart;
    const items = this.cellItems;
    const startA = cellStart[c];
    const endA = cellStart[c + 1];
    const startB = cellStart[other];
    const endB = cellStart[other + 1];
    if (startB === endB) {
      return;
    }

    for (let a = startA; a < endA; a++) {
      const i = items[a];
      if (!this.alive[i]) {
        continue;
      }
      for (let b = startB; b < endB; b++) {
        const j = items[b];
        if (!this.alive[j]) {
          continue;
        }
        this.interact(i, j);
        if (!this.alive[i] || this.bornAt[i] === this.stepCounter) {
          break;
        }
      }
    }
  }

  /**
   * Resolves one candidate pair: rejects it if the disks are not really
   * overlapping, otherwise reacts or bounces.
   */
  private interact(i: number, j: number): void {
    // A particle created during this step has not moved yet. It may still
    // bounce -- otherwise a product born touching a neighbor would stay
    // interpenetrating for a step -- but it may not react again until it
    // has, so a fresh dissociation cannot turn straight back into a dimer.
    const fresh = this.bornAt[i] === this.stepCounter || this.bornAt[j] === this.stepCounter;

    const sx = this.x[i] - this.x[j];
    const sy = this.y[i] - this.y[j];
    const contact = this.radius[i] + this.radius[j];
    // Cheap rejects before the multiply: most candidate pairs fail here.
    if (sx > contact || sx < -contact || sy > contact || sy < -contact) {
      return;
    }
    const distSq = sx * sx + sy * sy;
    if (distSq >= contact * contact) {
      return;
    }

    const ux = this.vx[i] - this.vx[j];
    const uy = this.vy[i] - this.vy[j];
    // Already separating: they were dealt with on approach (or were born
    // overlapping). Touching them again would inject energy from nowhere.
    if (sx * ux + sy * uy >= 0) {
      return;
    }

    const aIsA = this.species[i] === SPECIES_A;
    const bIsA = this.species[j] === SPECIES_A;

    if (fresh) {
      this.bounce(i, j, sx, sy, ux, uy, contact);
      return;
    }

    if (aIsA && bIsA) {
      // combine() declines if there is no room for the dimer, in which case
      // the pair just bounces -- an excluded-volume veto, not a lost event.
      if (Math.random() * 100 < this.params.dimerizationChance && this.combine(i, j, sx, sy, ux, uy, contact)) {
        return;
      }
    } else {
      // At least one B. Each B rolls independently, so a B+B contact can
      // break either, both, or neither.
      // dissociate() declines when the fragments have nowhere to go, in
      // which case the contact falls through to an ordinary bounce.
      let reacted = false;
      if (!aIsA && Math.random() * 100 < this.params.dissociationChance && this.dissociate(i)) {
        reacted = true;
      }
      if (!bIsA && this.alive[j] && Math.random() * 100 < this.params.dissociationChance && this.dissociate(j)) {
        reacted = true;
      }
      if (reacted) {
        return;
      }
    }

    this.bounce(i, j, sx, sy, ux, uy, contact);
  }

  /**
   * Solves for the (negative) time at which two overlapping disks were
   * exactly touching. Returns 0 if the geometry is degenerate -- deeply
   * overlapping or barely moving -- in which case the caller falls back to
   * the current line of centers, which is still the right normal direction.
   */
  private contactTime(sx: number, sy: number, ux: number, uy: number, contact: number): number {
    const a = ux * ux + uy * uy;
    if (a <= 0) {
      return 0;
    }
    const b = 2 * (sx * ux + sy * uy);
    const c = sx * sx + sy * sy - contact * contact;
    const disc = b * b - 4 * a * c;
    if (disc < 0) {
      return 0;
    }
    const root = Math.sqrt(disc);
    const t1 = (-b - root) / (2 * a);
    const t2 = (-b + root) / (2 * a);
    const t = Math.abs(t1) < Math.abs(t2) ? t1 : t2;
    if (!Number.isFinite(t) || t > 0) {
      return 0;
    }
    // A rewind longer than a whole step means the pair was already deeply
    // tangled; replaying it would fling them apart. Take the normal as-is.
    if (-t > this.tickDelta) {
      return 0;
    }
    return t;
  }

  /**
   * An exactly elastic collision between two hard disks: rewind to contact,
   * reverse the relative velocity along the line of centers in the
   * center-of-mass frame, replay the remainder of the step. Momentum and
   * kinetic energy are both conserved exactly.
   */
  private bounce(
    i: number,
    j: number,
    sx: number,
    sy: number,
    ux: number,
    uy: number,
    contact: number,
  ): void {
    const t = this.contactTime(sx, sy, ux, uy, contact);

    // Rewind both particles to the instant of contact.
    this.x[i] += t * this.vx[i];
    this.y[i] += t * this.vy[i];
    this.x[j] += t * this.vx[j];
    this.y[j] += t * this.vy[j];

    // The true contact normal, from j to i.
    let nx = this.x[i] - this.x[j];
    let ny = this.y[i] - this.y[j];
    let len = Math.sqrt(nx * nx + ny * ny);
    if (len <= 1e-12) {
      // Exactly coincident centers: any direction is as good as another.
      const theta = Math.random() * 2 * Math.PI;
      nx = Math.cos(theta);
      ny = Math.sin(theta);
      len = 1;
    }
    nx /= len;
    ny /= len;

    const mi = this.mass[i];
    const mj = this.mass[j];
    // Relative velocity along the normal. Negative means approaching.
    const vn = ux * nx + uy * ny;
    if (vn < 0) {
      // Reversing the normal component of the relative velocity, split
      // between the two by inverse mass, is the elastic impulse.
      const jImp = (-2 * vn) / (1 / mi + 1 / mj);
      this.vx[i] += (jImp / mi) * nx;
      this.vy[i] += (jImp / mi) * ny;
      this.vx[j] -= (jImp / mj) * nx;
      this.vy[j] -= (jImp / mj) * ny;
    }

    // Replay the part of the step that happened after contact.
    this.x[i] -= t * this.vx[i];
    this.y[i] -= t * this.vy[i];
    this.x[j] -= t * this.vx[j];
    this.y[j] -= t * this.vy[j];

    this.clampIntoBox(i);
    this.clampIntoBox(j);

    // A degenerate rewind (contactTime() gave up and returned 0 on a deep
    // overlap) leaves the pair still interpenetrating. Measure the gap that
    // is actually left and push it out along the normal, split by inverse
    // mass, so the next step starts clean. After a normal rewind/replay
    // there is nothing left to do and this costs one distance test.
    const rx = this.x[i] - this.x[j];
    const ry = this.y[i] - this.y[j];
    const gap = contact - Math.sqrt(rx * rx + ry * ry);
    if (gap > 0) {
      const total = mi + mj;
      this.x[i] += nx * gap * (mj / total);
      this.y[i] += ny * gap * (mj / total);
      this.x[j] -= nx * gap * (mi / total);
      this.y[j] -= ny * gap * (mi / total);
      this.clampIntoBox(i);
      this.clampIntoBox(j);
    }
  }

  private clampIntoBox(i: number): void {
    const reach = this.innerBound - this.radius[i];
    if (this.x[i] > reach) {
      this.x[i] = reach;
    } else if (this.x[i] < -reach) {
      this.x[i] = -reach;
    }
    if (this.y[i] > reach) {
      this.y[i] = reach;
    } else if (this.y[i] < -reach) {
      this.y[i] = -reach;
    }
  }

  /**
   * 2A -> B. The dimer carries away the pair's entire kinetic energy: there
   * is no bond energy and no thermal bath, so a reaction moves energy around
   * without creating or destroying any.
   *
   * Momentum is *not* conserved here, and cannot be. Momentum conservation
   * alone fixes the dimer's velocity at the pair's center-of-mass velocity,
   * whose kinetic energy is short of the reactants' by exactly the relative
   * term 1/2 mu v_rel^2 -- strictly positive for any pair close enough to
   * react. Energy and momentum cannot both survive a two-body association,
   * which is why real recombination needs a third body. This model keeps the
   * energy and lets the momentum go: the dimer leaves along the direction
   * the pair's momentum pointed, at whatever speed carries the full energy.
   */
  private combine(
    i: number,
    j: number,
    sx: number,
    sy: number,
    ux: number,
    uy: number,
    contact: number,
  ): boolean {
    // Fuse at the point where they actually touched, not where the
    // discrete step happened to leave them overlapping.
    const t = this.contactTime(sx, sy, ux, uy, contact);
    const xi = this.x[i] + t * this.vx[i];
    const yi = this.y[i] + t * this.vy[i];
    const xj = this.x[j] + t * this.vx[j];
    const yj = this.y[j] + t * this.vy[j];

    const mi = this.mass[i];
    const mj = this.mass[j];
    const massB = mi + mj;

    // All of the reactants' kinetic energy, carried across whole.
    const energy = 0.5 * (mi * (this.vx[i] * this.vx[i] + this.vy[i] * this.vy[i])
      + mj * (this.vx[j] * this.vx[j] + this.vy[j] * this.vy[j]));
    const speed = Math.sqrt((2 * energy) / massB);

    // Heading: where the pair's momentum was already pointing, so the dimer
    // continues the motion the eye was following. A head-on pair of equals
    // has no such direction, so that case picks one at random.
    let dirX = mi * this.vx[i] + mj * this.vx[j];
    let dirY = mi * this.vy[i] + mj * this.vy[j];
    const dirLen = Math.hypot(dirX, dirY);
    if (dirLen > 1e-12) {
      dirX /= dirLen;
      dirY /= dirLen;
    } else {
      const theta = Math.random() * 2 * Math.PI;
      dirX = Math.cos(theta);
      dirY = Math.sin(theta);
    }
    const dimerVx = speed * dirX;
    const dimerVy = speed * dirY;
    // Mass-weighted midpoint: the center of mass, which is where the dimer
    // belongs and what keeps the pair's motion continuous.
    let cx = (mi * xi + mj * xj) / massB;
    let cy = (mi * yi + mj * yj) / massB;
    // Nudge the dimer off a wall before testing for room, so a reaction is
    // not vetoed by a position that clamping was going to fix anyway.
    const reach = this.innerBound - B_RADIUS;
    if (cx > reach) {
      cx = reach;
    } else if (cx < -reach) {
      cx = -reach;
    }
    if (cy > reach) {
      cy = reach;
    } else if (cy < -reach) {
      cy = -reach;
    }

    // A B is wider than the A's it replaces, so it can stick out
    // into a bystander that neither A was touching. Rather than let it be
    // born interpenetrating -- which looks wrong and can only be undone by
    // the two drifting back apart -- decline the reaction when the space
    // is not free. Physically this is excluded volume: no room, no dimer.
    if (this.wouldOverlap(cx, cy, B_RADIUS, i, j)) {
      return false;
    }

    this.kill(i);
    this.kill(j);
    this.spawn(SPECIES_B, massB, cx, cy, dimerVx, dimerVy);
    return true;
  }

  /**
   * Is a disk of this radius at this point clear of every living particle
   * except the two being consumed? Scans the 3x3 block of grid cells around
   * the point, which covers everything that could reach it.
   */
  private wouldOverlap(px: number, py: number, r: number, skipI: number, skipJ: number): boolean {
    const n = this.gridN;
    const invCell = 1 / this.cellSize;
    let cx = ((px - this.gridOrigin) * invCell) | 0;
    let cy = ((py - this.gridOrigin) * invCell) | 0;
    if (cx < 0) {
      cx = 0;
    } else if (cx >= n) {
      cx = n - 1;
    }
    if (cy < 0) {
      cy = 0;
    } else if (cy >= n) {
      cy = n - 1;
    }

    const cellStart = this.cellStart;
    const items = this.cellItems;
    for (let gy = cy - 1; gy <= cy + 1; gy++) {
      if (gy < 0 || gy >= n) {
        continue;
      }
      for (let gx = cx - 1; gx <= cx + 1; gx++) {
        if (gx < 0 || gx >= n) {
          continue;
        }
        const c = gx + gy * n;
        for (let k = cellStart[c]; k < cellStart[c + 1]; k++) {
          const o = items[k];
          if (o === skipI || o === skipJ || !this.alive[o]) {
            continue;
          }
          const dx = this.x[o] - px;
          const dy = this.y[o] - py;
          const contact = this.radius[o] + r;
          if (dx * dx + dy * dy < contact * contact) {
            return true;
          }
        }
      }
    }
    return false;
  }

  /**
   * B -> 2A. The dimer's kinetic energy is split at random between the two
   * fragments, each of which then flies off in its own random direction.
   * Nothing is drawn from a bath and nothing is spent on a bond, so the
   * gas's total kinetic energy is exactly what it was.
   *
   * As in combine(), momentum does not survive this: fixing the two
   * fragments' momenta to sum to the parent's would force their relative
   * kinetic energy to zero, leaving a pair that never separates. Energy is
   * kept and momentum is allowed to float.
   *
   * Returns false when there is nowhere to put the fragments, in which case
   * nothing changes and the caller should treat the contact as an ordinary
   * bounce -- the mirror of combine()'s excluded-volume veto.
   */
  private dissociate(b: number): boolean {
    const childMass = this.mass[b] / 2;
    const energy = 0.5 * this.mass[b] * (this.vx[b] * this.vx[b] + this.vy[b] * this.vy[b]);
    const x = this.x[b];
    const y = this.y[b];

    // Born touching rather than coincident, offset along the kick axis, so
    // the pair starts out exactly at contact and already moving apart.
    const offset = A_RADIUS + 1e-6;
    const reach = this.innerBound - A_RADIUS;

    // Two A's born touching span 2*A_RADIUS from the parent's center, but a
    // B only reaches A_RADIUS*cbrt(2) -- the radius that puts the dimer at
    // twice a monomer's volume. A split therefore claims space the dimer
    // never occupied, which a bystander is entitled to be sitting in, so the
    // room has to be checked rather than assumed. The axis is random, and a
    // B blocked along one line usually has room along another, so decline
    // only after several draws have all failed.
    let offX = 0;
    let offY = 0;
    let room = false;
    for (let tries = 0; tries < SPLIT_ORIENTATION_TRIES; tries++) {
      const theta = Math.random() * 2 * Math.PI;
      offX = Math.cos(theta) * offset;
      offY = Math.sin(theta) * offset;
      const x1 = x + offX;
      const y1 = y + offY;
      const x2 = x - offX;
      const y2 = y - offY;
      // Reject an axis that would put a fragment through a wall rather than
      // clamping it back: clamping is what would slide it into a neighbor.
      if (x1 > reach || x1 < -reach || y1 > reach || y1 < -reach
        || x2 > reach || x2 < -reach || y2 > reach || y2 < -reach) {
        continue;
      }
      if (!this.wouldOverlap(x1, y1, A_RADIUS, b, b) && !this.wouldOverlap(x2, y2, A_RADIUS, b, b)) {
        room = true;
        break;
      }
    }
    if (!room) {
      return false;
    }

    // Split the parent's energy at random between the fragments, then send
    // each off on its own heading. The placement axis above and these two
    // headings are independent, so a pair can be born heading back into one
    // another; the freshness guard makes that an ordinary elastic bounce,
    // which costs nothing in energy.
    const share = Math.random();
    const speed1 = Math.sqrt((2 * share * energy) / childMass);
    const speed2 = Math.sqrt((2 * (1 - share) * energy) / childMass);
    const theta1 = Math.random() * 2 * Math.PI;
    const theta2 = Math.random() * 2 * Math.PI;

    this.kill(b);
    const a1 = this.spawn(SPECIES_A, childMass, x + offX, y + offY,
      speed1 * Math.cos(theta1), speed1 * Math.sin(theta1));
    const a2 = this.spawn(SPECIES_A, childMass, x - offX, y - offY,
      speed2 * Math.cos(theta2), speed2 * Math.sin(theta2));
    this.clampIntoBox(a1);
    this.clampIntoBox(a2);
    return true;
  }

  // =====================================================================
  // Readouts
  // =====================================================================

  /** Interior side, wall face to wall face, in patches. */
  get interiorSide(): number {
    return 2 * this.innerBound;
  }

  /** Interior area in nm^2 -- the "volume" a 2D concentration divides by. */
  get areaNm2(): number {
    const side = this.interiorSide * NM_PER_PATCH;
    return side * side;
  }

  /**
   * Partial pressures in N/m, one per species, averaged over the simulated
   * time since the previous call -- which also resets the accumulator, so
   * successive calls tile the run rather than overlapping.
   *
   * Returns null when no time has passed, so a caller polling faster than
   * the sim steps gets "no reading yet" instead of a divide by zero.
   */
  drainWallPressure(): { a: number; b: number } | null {
    const elapsed = this.wallImpulseTicks;
    if (elapsed <= 0) {
      return null;
    }
    // Force per unit length: impulse / time / perimeter.
    const perimeter = 4 * this.interiorSide;
    const scale = PRESSURE_TO_N_PER_M / (elapsed * perimeter);
    const out = {
      a: this.wallImpulse[SPECIES_A] * scale,
      b: this.wallImpulse[SPECIES_B] * scale,
    };
    this.wallImpulse[SPECIES_A] = 0;
    this.wallImpulse[SPECIES_B] = 0;
    this.wallImpulseTicks = 0;
    return out;
  }

  stats(): Stats {
    let countA = 0;
    let countB = 0;
    let energy = 0;
    let px = 0;
    let py = 0;
    for (let i = 0; i < this.top; i++) {
      if (!this.alive[i]) {
        continue;
      }
      if (this.species[i] === SPECIES_A) {
        countA++;
      } else {
        countB++;
      }
      const m = this.mass[i];
      const vx = this.vx[i];
      const vy = this.vy[i];
      energy += 0.5 * m * (vx * vx + vy * vy);
      px += m * vx;
      py += m * vy;
    }
    const n = countA + countB;
    const area = this.areaNm2;
    // Concentrations are per nm^2, so Kc = [B]/[A]^2 comes out in nm^2 --
    // a real, dimensional constant that no longer has the box size baked
    // into it the way a bare count ratio did.
    const concA = countA / area;
    const concB = countB / area;
    const temperature = n > 0 ? energy / n : 0;
    return {
      countA,
      countB,
      totalAEquivalent: countA + 2 * countB,
      measuredTemperature: temperature,
      kc: concA > 0 ? concB / (concA * concA) : 0,
      areaNm2: area,
      concA,
      concB,
      kineticEnergy: energy,
      // Reactions do not conserve momentum in this model, so this wanders
      // instead of sitting at zero. Reported as a speed the whole gas would
      // have if its momentum were shared out, which is directly comparable
      // to a thermal speed and so says how much of the pressure is drift.
      driftSpeed: n > 0 ? Math.hypot(px, py) / this.totalMass() : 0,
      // What an ideal 2D gas would push with at this count, temperature and
      // area: Pi = N k_B T / A. The wall gauge is measured independently,
      // so the two agreeing is a result rather than a definition.
      idealPressure: area > 0
        ? ((n * temperature) / (this.interiorSide * this.interiorSide)) * PRESSURE_TO_N_PER_M
        : 0,
    };
  }

  /** Total mass of the living gas, in amu. */
  private totalMass(): number {
    let m = 0;
    for (let i = 0; i < this.top; i++) {
      if (this.alive[i]) {
        m += this.mass[i];
      }
    }
    return m;
  }

  /**
   * Bins the two species' speeds into `bins` buckets spanning [0, hiSpeed],
   * adding into the caller's arrays. Returns the fastest particle so the
   * caller can keep the axis sized sensibly.
   */
  accumulateSpeeds(outA: Float64Array, outB: Float64Array, hiSpeed: number): number {
    const bins = outA.length;
    const scale = bins / hiSpeed;
    let fastest = 0;
    for (let i = 0; i < this.top; i++) {
      if (!this.alive[i]) {
        continue;
      }
      const speed = Math.sqrt(this.vx[i] * this.vx[i] + this.vy[i] * this.vy[i]);
      if (speed > fastest) {
        fastest = speed;
      }
      let bin = (speed * scale) | 0;
      if (bin >= bins) {
        bin = bins - 1;
      }
      if (this.species[i] === SPECIES_A) {
        outA[bin]++;
      } else {
        outB[bin]++;
      }
    }
    return fastest;
  }

  /** Read-only views for the renderer. Valid until the next setup(). */
  view() {
    return {
      x: this.x,
      y: this.y,
      radius: this.radius,
      species: this.species,
      alive: this.alive,
      top: this.top,
    };
  }
}
