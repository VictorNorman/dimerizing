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
**B** (red, twice the mass and twice the radius, so it gets hit more often).

- Two A particles whose disks actually touch fuse into one B with probability
  `dimerization-chance`. Momentum is conserved exactly; the pair's kinetic
  energy *relative* to their shared center of mass is not carried forward —
  forming a bond is exothermic, and that energy went into the bond.
- A B that touches anything — a particle or a wall — falls apart into two A's
  with probability `dissociation-chance`. Momentum is conserved, and the
  fragments fly apart with a random-direction relative velocity sized off the
  current `temperature`: bond-breaking draws its energy from an ambient thermal
  bath.
- Everything else is an ordinary elastic collision, exactly conserving both
  momentum and kinetic energy.

So the gas is quasi-isothermal (coupled to a bath at `temperature`) rather than
energy-closed — which is the right physical picture for what "B/A² should be
constant at equilibrium" describes: a constant-temperature reaction vessel, not
an insulated one. Mass is conserved exactly and always: `A + 2B` never changes.

The point of the model is **dynamic equilibrium**. Individual molecules keep
reacting in both directions forever, but the populations settle into a stable
balance where the forward and reverse rates match, and `Kc = B/A²`
holds steady however you got there. (The vessel's size is part of the
constant because a mass-action Kc is built from concentrations, and
concentration is what changes when you resize the box at a fixed particle
count — that's why the slider is a percentage of the world's *area*, the
volume term a Kc actually contains.) Try `initial-B = 0` against a run that
starts with B already present — they land on the same Kc.

### Units

`temperature` and `particle-mass` are labeled K and amu because that is the
intuitive way to think about them, but this is not a dimensionally real physics
model. There is no Boltzmann constant: temperature *is* average kinetic energy
per particle (k_B = 1, the usual reduced-units convention in molecular
simulation). Box size is a genuine dimensionless percentage of the world's *area*
-- the two-dimensional box shrinks each side by sqrt(percent/100), so the
vessel's area lands on the requested percentage.

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
- **Excluded volume on dimerization.** A B is twice the radius of the two A's
  it replaces, so it can be born sticking into a bystander neither A was
  touching. When there is no room, the reaction is declined and the pair simply
  bounces. Without this, dimers get born interpenetrating and can only
  untangle by drifting apart.
- **The speed histogram** pools the newest 5 ticks' speeds into a fixed bin
  grid and shows each bar as the *mean* count per bin across that window,
  matching the Python model's `speed_average_ticks` time average (rather than
  the original's pooled window whose bin edges moved every tick).
- **Larger population limits** (up to 4,000 A), since the speed allows it.
- **A simulation-speed control**, in simulated ticks per frame. Each frame
  spends at most 11 ms stepping, so a population big enough to blow that budget
  slows down smoothly instead of locking up the page.

`box size`, `initial A`, `initial B` and `particle-mass` are read only by
Setup, exactly as in the original; the panel marks them "on setup". Everything
else takes effect live.

## What the tests check

`npm run check` runs the model headlessly and verifies:

- `A + 2B` is conserved exactly over 20,000 steps
- no particle ever leaves the box
- particles never pass through one another (disk overlap stays below a full
  contact distance), overlap is rare (~30 per million particle-steps), and no
  overlapping pair ever stops resolving — nothing jams
- with reactions off, kinetic energy is conserved to 1 part in 10¹⁴
- speeds relax to a 2D Maxwell-Boltzmann distribution — mean kinetic energy
  holds at the set temperature, and sd/mean converges on √(4/π − 1) ≈ 0.523
- Kc reaches the same equilibrium starting from all-A and from all-B

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
