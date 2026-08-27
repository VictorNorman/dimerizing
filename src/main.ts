import './styles.css';
import { Sim, A_RADIUS, B_RADIUS, type Params, type Stats } from './sim';
import { Renderer, COLOR_A, COLOR_B } from './render';
import { HistogramPlot, TimeSeriesPlot } from './plots';

// =======================================================================
// Controls
// =======================================================================

interface SliderDef {
  key: keyof Params | 'ticksPerSecond';
  label: string;
  min: number;
  max: number;
  step: number;
  value: number;
  units?: string;
  /** Read only by setup(); drags are held in `pending` until the next Initiate. */
  deferred?: boolean;
  format?: (v: number) => string;
}

const SLIDERS: SliderDef[] = [
  {
    key: 'initialA',
    label: 'initial number of A',
    min: 0,
    max: 4000,
    step: 10,
    value: 1000,
    deferred: true,
  },
  {
    key: 'initialB',
    label: 'initial number of B',
    min: 0,
    max: 2000,
    step: 5,
    value: 0,
    deferred: true,
  },
  {
    key: 'particleMass',
    label: 'particle-mass (A)',
    min: 2,
    max: 100,
    step: 1,
    value: 44,
    units: 'amu',
    deferred: true,
  },
  {
    key: 'temperature',
    label: 'initial temperature',
    min: 1,
    max: 1000,
    step: 1,
    value: 298,
    units: 'K',
  },
  {
    key: 'dimerizationChance',
    label: 'dimerization probability',
    min: 0,
    max: 100,
    step: 1,
    value: 40,
    units: '%',
  },
  {
    key: 'dissociationChance',
    label: 'dissociation probability',
    min: 0,
    max: 100,
    step: 1,
    value: 4,
    units: '%',
  },
  {
    key: 'boxSize',
    label: 'operating area',
    min: 5,
    max: 100,
    step: 1,
    value: 100,
    units: '% of maximum area',
    deferred: true,
  },
  {
    key: 'ticksPerSecond',
    label: 'simulation speed',
    min: 5,
    max: 1200,
    step: 5,
    value: 25,
    units: 'ticks/s',
  },
];

const params: Params = {
  initialA: 1000,
  initialB: 0,
  particleMass: 44,
  temperature: 298,
  dimerizationChance: 40,
  dissociationChance: 4,
  boxSize: 100,
};

/** Values of the deferred sliders as dragged, applied at the next Initiate. */
const pending: Partial<Params> = {};
let ticksPerSecond = 25;

/**
 * Every slider that describes the experiment rather than the viewing of it.
 * Changing one of these mid-run either does nothing until the next Initiate
 * or silently redefines the conditions a run is already partway through, so
 * they are locked while the sim is going. `simulation speed` stays live: it
 * only sets how fast you watch.
 */
const lockedRows: { row: HTMLElement; input: HTMLInputElement }[] = [];

function buildSliders(host: HTMLElement): void {
  for (const def of SLIDERS) {
    const row = document.createElement('div');
    row.className = 'slider-row';

    const head = document.createElement('div');
    head.className = 'slider-head';
    const name = document.createElement('span');
    name.className = 'name';
    name.textContent = def.label;
    const val = document.createElement('span');
    val.className = 'val';
    head.append(name, val);

    const input = document.createElement('input');
    input.type = 'range';
    input.min = String(def.min);
    input.max = String(def.max);
    input.step = String(def.step);
    input.value = String(def.value);

    const show = (v: number) => {
      const text = def.format ? def.format(v) : String(v);
      if (!def.units) {
        val.textContent = text;
        return;
      }
      // A percent sign binds to its number -- "100% of maximum area", not
      // "100 % of ...". Every other unit reads as a separate word.
      const gap = def.units.startsWith('%') ? '' : ' ';
      val.textContent = `${text}${gap}${def.units}`;
    };
    show(def.value);

    input.addEventListener('input', () => {
      const v = Number(input.value);
      show(v);
      if (def.key === 'ticksPerSecond') {
        ticksPerSecond = v;
      } else if (def.deferred) {
        pending[def.key] = v;
      } else {
        (params[def.key] as number) = v;
      }
    });

    row.append(head, input);
    host.append(row);
    if (def.key !== 'ticksPerSecond') {
      lockedRows.push({
        row,
        input,
      });
    }
  }
}

/** Locks the experiment sliders while the sim is running; see `lockedRows`. */
function setSlidersLocked(locked: boolean): void {
  for (const { row, input } of lockedRows) {
    input.disabled = locked;
    row.classList.toggle('locked', locked);
  }
}

// =======================================================================
// Wiring
// =======================================================================

const el = <T extends HTMLElement>(id: string): T => {
  const node = document.getElementById(id);
  if (!node) {
    throw new Error(`missing #${id}`);
  }
  return node as T;
};

buildSliders(el('sliders'));

/**
 * Diameter of the smaller species' swatch. The other is scaled off the real
 * radii, so the two dots beside the counts stand in the same 1 : cbrt(2)
 * proportion as the particles on the canvas.
 */
const DOT_A_PX = 10;

function paintDot(id: string, color: string, radius: number): void {
  const dot = el(id);
  const px = DOT_A_PX * (radius / A_RADIUS);
  dot.style.width = `${px}px`;
  dot.style.height = `${px}px`;
  dot.style.background = color;
}

paintDot('dot-a', COLOR_A, A_RADIUS);
paintDot('dot-b', COLOR_B, B_RADIUS);

const sim = new Sim(params);
const renderer = new Renderer(el<HTMLCanvasElement>('world'));

const SPEED_BINS = 48;
/**
 * How many ticks the Speed Distribution histogram averages over. Each sample
 * bins one tick's speeds into the shared 48-bin grid and the displayed bars
 * are the *mean* count per bin over the newest SPEED_AVERAGE_TICKS snapshots
 * -- a true pooled time average, as the Python update_variables() does,
 * rather than a per-tick count that would jump around tick to tick. The
 * Python model's speed_average_ticks was 5; 10 buys a visibly steadier
 * histogram at the cost of lagging a change in conditions twice as long.
 */
const SPEED_AVERAGE_TICKS = 10;

const popPlot = new TimeSeriesPlot(el<HTMLCanvasElement>('plot-pop'), {
  xLabel: 'time (ticks)',
  yLabel: 'count',
  pens: [
    {
      label: 'A',
      color: COLOR_A,
    },
    {
      label: 'B',
      color: COLOR_B,
    },
  ],
});

const kcPlot = new TimeSeriesPlot(el<HTMLCanvasElement>('plot-kc'), {
  xLabel: 'time (ticks)',
  yLabel: '[B] / [A]² (nm²)',
  pens: [{
    label: 'K',
    sub: 'c',
    color: '#9b82d8',
  }],
  includeZero: false,
});

const speedPlot = new HistogramPlot(el<HTMLCanvasElement>('plot-speed'), {
  xLabel: 'speed',
  yLabel: '# molecules',
  pens: [
    {
      label: 'A speed',
      color: 'rgba(87, 182, 232, 0.72)',
    },
    {
      label: 'B speed',
      color: 'rgba(224, 85, 97, 0.72)',
    },
  ],
  bins: SPEED_BINS,
});

const tmpA = new Float64Array(SPEED_BINS);
const tmpB = new Float64Array(SPEED_BINS);
/**
 * The last SPEED_AVERAGE_TICKS per-tick speed histograms, oldest first --
 * what the Speed Distribution plot averages over. Each entry is a freshly
 * allocated pair of bin counts for one tick.
 */
const speedHistory: { a: Float64Array; b: Float64Array }[] = [];
let hiSpeed = 1;

let running = false;
let lastStats: Stats = sim.stats();
/** The next integer tick at which the plots and monitors take a sample. */
let nextSampleTick = 1;

// ---- monitors ----------------------------------------------------------

const monA = el('mon-a');
const monB = el('mon-b');
const monTotal = el('mon-total');
const monTemp = el('mon-temp');
const monArea = el('mon-area');
const monKc = el('mon-kc');
const monPressureA = el('mon-pressure-a');
const monPressureB = el('mon-pressure-b');
const monPressure = el('mon-pressure');
const monPressureIdeal = el('mon-pressure-ideal');
const monPressureRatio = el('mon-pressure-ratio');
const monTicks = el('mon-ticks');
const perfParticles = el('perf-particles');

/**
 * How many ticks of wall impacts the pressure readout averages over. A
 * single tick's impulses are far too sporadic to read: pressure is a mean
 * over many collisions, so the gauge needs a window the way the speed
 * histogram does.
 */
const PRESSURE_AVERAGE_TICKS = 20;

const pressureHistory: { a: number; b: number }[] = [];

/** Mean partial pressures over the window, in N/m. */
function meanPressure(): { a: number; b: number } {
  if (pressureHistory.length === 0) {
    return {
      a: 0,
      b: 0,
    };
  }
  let a = 0;
  let b = 0;
  for (const p of pressureHistory) {
    a += p.a;
    b += p.b;
  }
  return {
    a: a / pressureHistory.length,
    b: b / pressureHistory.length,
  };
}

/** Three significant figures, matching the Python monitor's sig_figs=3. */
function sig3(v: number): string {
  if (!Number.isFinite(v)) {
    return '—';
  }
  if (v === 0) {
    return '0';
  }
  return Number(v.toPrecision(3)).toString();
}

function refreshMonitors(): void {
  monA.textContent = String(lastStats.countA);
  monB.textContent = String(lastStats.countB);
  monTotal.textContent = String(lastStats.totalAEquivalent);
  monTemp.textContent = sig3(lastStats.measuredTemperature);
  monArea.textContent = sig3(lastStats.areaNm2);
  monKc.textContent = sig3(lastStats.kc);

  // Pressures are shown in mN/m: a dilute 2D gas at these densities sits
  // around a few tenths, where a monolayer runs to tens.
  const p = meanPressure();
  const total = p.a + p.b;
  const ideal = lastStats.idealPressure;
  monPressureA.textContent = sig3(p.a * 1e3);
  monPressureB.textContent = sig3(p.b * 1e3);
  monPressure.textContent = sig3(total * 1e3);
  monPressureIdeal.textContent = sig3(ideal * 1e3);
  monPressureRatio.textContent = ideal > 0 && total > 0 ? (total / ideal).toFixed(3) : '—';

  monTicks.textContent = sim.ticks.toFixed(1);
}

// ---- sampling ----------------------------------------------------------

/** One plot sample: the same once-per-whole-tick cadence the Python used. */
function sample(): void {
  lastStats = sim.stats();
  popPlot.push(sim.ticks, [lastStats.countA, lastStats.countB]);
  kcPlot.push(sim.ticks, [lastStats.kc]);

  // Drain the wall gauge on the same once-per-tick cadence, so the window
  // below covers PRESSURE_AVERAGE_TICKS of simulated time.
  const p = sim.drainWallPressure();
  if (p) {
    pressureHistory.push(p);
    while (pressureHistory.length > PRESSURE_AVERAGE_TICKS) {
      pressureHistory.shift();
    }
  }

  tmpA.fill(0);
  tmpB.fill(0);
  const fastest = sim.accumulateSpeeds(tmpA, tmpB, hiSpeed);

  // Pool this tick into the rolling window, keeping only the newest
  // SPEED_AVERAGE_TICKS snapshots -- the Python model's del
  // speed_history[:-speed_average_ticks].
  speedHistory.push({ a: tmpA.slice(), b: tmpB.slice() });
  while (speedHistory.length > SPEED_AVERAGE_TICKS) {
    speedHistory.shift();
  }

  const bars = speedPlot.binArrays;
  // Each bar is the *mean* count per bin across the window -- scale=1/len
  // -- so early in a run, when the window is short, the bars still read as
  // real molecule counts rather than starting out a windowful too short.
  const windowScale = 1 / speedHistory.length;
  for (let i = 0; i < SPEED_BINS; i++) {
    let sumA = 0;
    let sumB = 0;
    for (const snap of speedHistory) {
      sumA += snap.a[i];
      sumB += snap.b[i];
    }
    bars[0][i] = sumA * windowScale;
    bars[1][i] = sumB * windowScale;
  }

  // Keep the axis just above the fastest particle. It snaps outward at once
  // so nothing piles into the last bin, and creeps back in slowly so the
  // bars are not sliding across a moving grid every tick.
  const target = Math.max(fastest * 1.05, 1e-6);
  hiSpeed = target > hiSpeed ? target : hiSpeed + (target - hiSpeed) * 0.02;
  speedPlot.setRange(hiSpeed);
}

function resetDerived(): void {
  popPlot.clear();
  kcPlot.clear();
  speedPlot.clear();
  speedHistory.length = 0;
  pressureHistory.length = 0;
  sim.drainWallPressure();
  nextSampleTick = 1;
  hiSpeed = 3 * Math.sqrt((2 * params.temperature) / params.particleMass);
  sample();
  refreshMonitors();
  drawPlots();
}

function drawPlots(): void {
  popPlot.draw();
  kcPlot.draw();
  speedPlot.draw();
}

// ---- buttons -----------------------------------------------------------

const goButton = el<HTMLButtonElement>('btn-go');

function setRunning(on: boolean): void {
  running = on;
  goButton.textContent = on ? 'Stop' : 'Go';
  goButton.classList.toggle('running', on);
  setSlidersLocked(on);
}

el('btn-setup').addEventListener('click', () => {
  // The sim holds this same object, which is how the live sliders take
  // effect mid-run; the deferred ones only land here.
  Object.assign(params, pending);
  sim.setup();
  resetDerived();
});

goButton.addEventListener('click', () => setRunning(!running));

el('btn-step').addEventListener('click', () => {
  setRunning(false);
  sim.step();
  if (sim.ticks >= nextSampleTick) {
    sample();
    nextSampleTick = Math.floor(sim.ticks) + 1;
  }
  refreshMonitors();
  drawPlots();
});

// ---- animation loop ----------------------------------------------------

/**
 * Share of one display frame that stepping may spend. The rest belongs to
 * rendering and to the browser, so a population big enough to blow the
 * budget slows down smoothly instead of locking the page up.
 */
const STEP_BUDGET_SHARE = 0.7;

/**
 * Floor and ceiling on that budget. The ceiling is what it has always been,
 * sized for a 60 Hz frame; the floor keeps a very fast panel from starving
 * the sim of steps entirely.
 */
const MIN_STEP_BUDGET_MS = 3;
const MAX_STEP_BUDGET_MS = 11;

/**
 * A frame gap longer than this is a stall -- a backgrounded tab, a blocked
 * main thread -- not elapsed simulation time. Advancing by the true gap
 * would make the sim lurch to catch up, so cap what one frame may claim.
 */
const MAX_FRAME_MS = 100;

/**
 * Estimated display refresh interval. The fastest recent frame is the best
 * evidence: a loaded frame can only run longer than the panel allows, never
 * shorter. It drifts back up slowly so moving the window to a 50 Hz monitor
 * is picked up within a second or two.
 */
let refreshMs = 1000 / 60;

let lastFrameMs = 0;
let measureStart = performance.now();

// The timestamp is always supplied by requestAnimationFrame; the default
// keeps the loop honest under a harness that calls it bare.
function frame(now: number = performance.now()): void {
  requestAnimationFrame(frame);

  const gap = lastFrameMs === 0 ? 0 : now - lastFrameMs;
  lastFrameMs = now;
  if (gap > 0 && gap < MAX_FRAME_MS) {
    refreshMs = gap < refreshMs ? gap : refreshMs + (gap - refreshMs) * 0.01;
  }

  if (running) {
    const budget = Math.min(
      MAX_STEP_BUDGET_MS,
      Math.max(MIN_STEP_BUDGET_MS, refreshMs * STEP_BUDGET_SHARE),
    );
    const deadline = performance.now() + budget;
    // Ticks are paced against the clock, not the frame, so the slider means
    // the same speed on a 120 Hz laptop panel and a 50 Hz external monitor.
    const target = sim.ticks + ticksPerSecond * (Math.min(gap, MAX_FRAME_MS) / 1000);
    let sinceCheck = 0;
    while (sim.ticks < target) {
      sim.step();
      if (sim.ticks >= nextSampleTick) {
        sample();
        nextSampleTick = Math.floor(sim.ticks) + 1;
      }
      // performance.now() is not free, and with a small population a step
      // can be cheaper than the clock read.
      if (++sinceCheck >= 8) {
        sinceCheck = 0;
        if (performance.now() >= deadline) {
          break;
        }
      }
    }
    refreshMonitors();
    drawPlots();
  }

  renderer.draw(sim);

  // The population changes only through reactions, so refreshing it on a
  // 500 ms window rather than every frame keeps the text from flickering.
  if (now - measureStart >= 500) {
    perfParticles.textContent = `${sim.count.toLocaleString()} particles`;
    measureStart = now;
  }
}

// ---- resize ------------------------------------------------------------

let resizeTimer = 0;
window.addEventListener('resize', () => {
  window.clearTimeout(resizeTimer);
  resizeTimer = window.setTimeout(() => {
    renderer.resize();
    popPlot.resize();
    kcPlot.resize();
    speedPlot.resize();
    drawPlots();
    renderer.draw(sim);
  }, 100);
});

resetDerived();
requestAnimationFrame(frame);
