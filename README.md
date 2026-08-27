# Dimerizing Gas — 2A ⇌ B

A 2D hard-disk gas in a box with a reversible dimerization reaction, running on
a canvas at interactive speed. This is a rewrite of `models/dimerizing_gas.py`
from [netlogo-py](../netlogo-py), which is faithful to the model but roughly
**15× faster per unit of simulated time** (and ~43× per step) while resolving
collisions *more* accurately than the original.

Inspired by [Paul Falstad's gas applet](https://www.falstad.com/gas/gas.html),
whose collision handling this borrows.

## Running it

```sh
npm install
npm run dev      # http://localhost:5173
```

```sh
npm run build    # typecheck + production bundle into dist/
npm test         # physics checks, then build, then a headless UI smoke test
```

`dist/` is a plain static bundle — no server needed, open `index.html` from
anywhere.

## The model

Two species share a box: **A** (blue, mass `particle-mass`, radius 0.25) and
**B** (red, twice the mass and twice the volume — so ∛2 ≈ 1.26 times the
radius, and it gets hit more often).

- Two A particles whose disks actually touch fuse into one B with probability
  `dimerization probability`. The dimer carries away the pair's **entire**
  kinetic energy, leaving along the direction their momentum pointed.
- A B that touches anything — a particle or a wall — falls apart into two A's
  with probability `dissociation probability`. The parent's kinetic energy is
  split at random between the fragments, each leaving on its own random
  heading.
- Everything else is an ordinary elastic collision, exactly conserving both
  momentum and kinetic energy.

There is no bond energy and no thermal bath: **a reaction moves energy around
without creating or destroying any**. Kinetic energy holds to 1 part in 10¹⁵
over 60,000 steps of reacting, which `npm run check` asserts. Mass is
conserved exactly and always: `A + 2B` never changes.

Momentum is *not* conserved at a reaction, and cannot be. Momentum
conservation alone pins the dimer's velocity to the pair's centre-of-mass
velocity, whose kinetic energy falls short of the reactants' by exactly
½μv_rel² — strictly positive for any pair close enough to react. Energy and
momentum cannot both survive a two-body association, which is why real
recombination needs a third body or a photon. This model keeps the energy and
lets the momentum go. The resulting drift stays small, a few percent of a
thermal speed, because the walls bound it.

Temperature is therefore a *result*, not a setting — the slider only says
where it starts. With energy fixed and T = E/N, pairing A's up **raises** the
temperature and splitting B's **lowers** it: from all-A at 298 K a run settles
near 500 K, from all-B near 257 K.

The point of the model is **dynamic equilibrium**. Individual molecules keep
reacting in both directions forever, but the populations settle into a stable
balance where the forward and reverse rates match, and `K_c = [B]/[A]²`
holds steady however you got there. (The vessel's size is part of the
constant because a mass-action K_c is built from concentrations, and
concentration is what changes when you resize the box at a fixed particle
count — that's why the slider is a percentage of the maximum *area*, the
volume term a K_c actually contains.) Try `initial number of B = 0` against a
run that starts with B already present — they land on the same K_c.

### Units

Fixing one length makes the whole model dimensional. The box interior at
100% `operating area` is defined as **100 nm** on a side, the temperature
slider is read as real kelvin and the mass slider as real amu, and every
other unit then follows — including the tick, which is not a free choice
once the other three are made:

| quantity | value |
|---|---|
| 1 patch | 1.266 nm |
| A diameter | 0.633 nm |
| B diameter | 0.797 nm |
| 1 tick | 13.9 ps |
| system area at 100% | 10,000 nm² |

The tick works out as √(amu / k_B) × 1 patch and depends on neither slider,
so it has one fixed duration however the gas is set up. A 400-tick run is
about 5.5 ns — real molecular-dynamics timescales, and A's diameter lands on
a real molecular size.

That makes the derived quantities real too. Concentrations are particles per
nm², so `K_c = [B]/[A]²` comes out in nm² — a dimensional constant with the
box size no longer baked into it, unlike the bare count ratio it replaces.

**Pressure** is measured mechanically. In 2D it is force per unit *length*
(N/m, the units of surface tension), and the gauge sums the momentum each
wall absorbs — 2m|v⊥| per reflection — then divides by perimeter and elapsed
time. Nothing about an equation of state is assumed, so comparing it against
N·k_B·T/A is a genuine test. Each species is accumulated separately, giving
partial pressures that add to the total by Dalton's law.

The ideal gas law does hold in two dimensions, to within the excluded-volume
correction: at 1000 particles the disks cover ~3% of the box and the measured
pressure runs about **7% above** N·k_B·T/A, close to the 2D hard-disk second
virial term 1 + 2φ. `npm run check` asserts that.

## What makes it fast

The Python original manages ~274 simulated ticks/s at 200 particles and ~51 at
1000. This version, single-threaded:

| particles | steps/s | sim-ticks/s | per particle per step |
|----------:|--------:|------------:|----------------------:|
|       200 | 109,000 |       4,400 |                 46 ns |
|     1,000 |  21,000 |         700 |                 47 ns |
|     4,000 |   3,900 |         120 |                 64 ns |

Five things account for it:

1. **Velocity components, not polar coordinates.** NetLogo stores motion as
   `(speed, heading)`, so every elastic collision costs four trig calls and
   every wall test a `sin`/`cos`. Storing `(vx, vy)` costs none; `speed` is
   derived with a `sqrt` only where it is genuinely needed — the histogram, the
   measured temperature, the timestep bound.

2. **Struct of arrays.** Particles live in parallel `Float64Array`s rather than
   objects, so a pass over the population is a linear walk over contiguous
   memory that the JIT keeps unboxed.

3. **A counting-sort uniform grid.** Neighbors are found through flat
   `Int32Array`s rebuilt each step in two linear passes — no hash map, no
   per-cell arrays, zero allocation per step. The Python version rebuilt a dict
   of lists every tick. Cells are sized to the gas's density (a little under
   one particle per cell, floored at the widest contact distance) instead of a
   fixed 1×1, which alone is worth 3× at 200 particles: a fixed grid spends
   most of a sparse run clearing 5,625 empty cells.

4. **Exact contact resolution.** When two disks are found overlapping, solve
   the quadratic for the exact time they touched, rewind both to that instant,
   apply the impulse along the true line of centers, then replay the rest of
   the step. This is Falstad's approach. The Python model instead applied the
   impulse along a *random* axis (inherited from NetLogo's GasLab) and needed a
   `last_collision` field on every particle to stop overlapping pairs from
   re-colliding forever — bookkeeping that exact resolution makes unnecessary.

5. **Batched rendering.** One canvas path per species, filled once, instead of
   a fill call per particle; sub-pixel particles are drawn as rects, which
   rasterize several times faster than arcs at that size.

## Where it deliberately differs from the Python model

- **Shorter timesteps.** The original stepped up to a full patch per tick,
  which lets two fast A particles (0.5 patches across) pass straight through
  one another undetected. This caps displacement at 0.2 patches per step, a 2×
  margin on the tightest case. That is ~5× more steps per unit simulated time —
  already priced into the table above.
- **Line-of-centers impulses**, as described above, rather than a random
  collision axis.
- **Excluded volume on both reactions.** A product can be born sticking into a
  bystander that no reactant was touching, so each direction checks for room
  first and declines the reaction when there is none — the contact just
  bounces. Forming a dimer tests the one point the B would occupy. Splitting
  one tests both fragment positions, and tries several random axes before
  giving up, since a B blocked along one line usually has room along another.
  The split needs the check more than the fusion does: two A's born touching
  span `2·A_RADIUS` from the parent's center while the B only reached
  `A_RADIUS·∛2`, so a split always claims ground the dimer never held.
  Without this, products get born interpenetrating and can only untangle by
  drifting apart.
- **The speed histogram** pools the newest 10 ticks' speeds into a fixed bin
  grid and shows each bar as the *mean* count per bin across that window — the
  Python model's `speed_average_ticks` time average, over twice its window of
  5 (rather than the original's pooled window whose bin edges moved every
  tick).
- **Larger population limits** (up to 4,000 A), since the speed allows it.
- **A simulation-speed control**, in simulated ticks per frame. Each frame
  spends at most 11 ms stepping, so a population big enough to blow that budget
  slows down smoothly instead of locking up the page.

`operating area`, `initial number of A`, `initial number of B` and
`particle-mass` are read only by Initiate, exactly as in the original.
Everything else, `initial temperature` included, takes effect live.

## What the tests check

`npm run check` runs the model headlessly and verifies:

- `A + 2B` is conserved exactly over 20,000 steps
- no particle ever leaves the box
- particles never pass through one another (disk overlap stays below a full
  contact distance), overlap is rare (~10 per million particle-steps), and no
  overlapping pair ever stops resolving — nothing jams
- kinetic energy is conserved to 1 part in 10¹⁵ *through reactions*, from
  both starting compositions, and temperature rises from all-A while falling
  from all-B
- wall-impulse pressure sits within the excluded-volume correction of
  N·k_B·T/A — the 2D ideal gas law, measured rather than assumed
- with reactions off, kinetic energy is conserved to 1 part in 10¹⁴
- speeds relax to a 2D Maxwell-Boltzmann distribution — mean kinetic energy
  holds at the set temperature, and sd/mean converges on √(4/π − 1) ≈ 0.523
- K_c reaches the same equilibrium starting from all-A and from all-B

`npm run smoke` loads the built bundle against a stub DOM and drives the real
animation loop, checking that every element resolves, the buttons work, and the
monitors and plots update.

## Layout

```
index.html        page shell and controls
src/sim.ts        the gas: storage, grid, collisions, reactions
src/render.ts     particle view
src/plots.ts      scrolling line plot + histogram, on canvas
src/main.ts       control wiring and the animation loop
src/styles.css
scripts/check.ts  headless physics checks and benchmark
scripts/smoke.mjs headless UI smoke test of the built bundle
```

No runtime dependencies.
