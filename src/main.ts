import './styles.css';
import { Sim, type Params, type Stats } from './sim';
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
  /** Read only by setup(), so the UI says so rather than silently ignoring drags. */
  deferred?: boolean;
  format?: (v: number) => string;
}

const SLIDERS: SliderDef[] = [
  {
    key: 'initialA',
    label: 'initial number of A',
    min: 0,
    max: 500,
    step: 2,
    value: 200,
    deferred: true,
  },
  {
    key: 'initialB',
    label: 'initial number of B',
    min: 0,
    max: 250,
    step: 1,
    value: 0,
    deferred: true,
  },
  {
    key: 'particleMass',
    label: 'particle-mass (A)',
    min: 1,
    max: 10,
    step: 1,
    value: 2,
    units: 'amu',
    deferred: true,
  },
  {
    key: 'temperature',
    label: 'temperature',
    min: 1,
    max: 500,
    step: 1,
    value: 8,
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
    label: 'box size',
    min: 5,
    max: 100,
    step: 1,
    value: 100,
    units: '% of world area',
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
  initialA: 200,
  initialB: 0,
  particleMass: 2,
  temperature: 8,
  dimerizationChance: 40,
  dissociationChance: 4,
  boxSize: 100,
};

/** Values of the deferred sliders as dragged, applied at the next Setup. */
const pending: Partial<Params> = {};
let ticksPerSecond = 25;

function buildSliders(host: HTMLElement): void {
  for (const def of SLIDERS) {
    const row = document.createElement('div');
    row.className = 'slider-row' + (def.deferred ? ' deferred' : '');

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
      val.textContent = def.units ? `${text} ${def.units}` : text;
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

const sim = new Sim(params);
const renderer = new Renderer(el<HTMLCanvasElement>('world'));

const SPEED_BINS = 48;
/**
 * How many ticks the Speed Distribution histogram averages over, matching
 * the Python model's speed_average_ticks = 5. Each sample bins one tick's
 * speeds into the shared 48-bin grid and the displayed bars are the *mean*
 * count per bin over the newest SPEED_AVERAGE_TICKS snapshots -- a true
 * pooled time average, exactly what the Python update_variables() does,
 * rather than a per-tick count that would jump around tick to tick.
 */
const SPEED_AVERAGE_TICKS = 5;

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
  yLabel: 'B / A² (×10⁴)',
  pens: [{
    label: 'Kc',
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
const monKc = el('mon-kc');
const monTicks = el('mon-ticks');
const perfSteps = el('perf-steps');
const perfFps = el('perf-fps');
const perfParticles = el('perf-particles');

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
  monKc.textContent = sig3(lastStats.kc);
  monTicks.textContent = sim.ticks.toFixed(1);
}

// ---- sampling ----------------------------------------------------------

/** One plot sample: the same once-per-whole-tick cadence the Python used. */
function sample(): void {
  lastStats = sim.stats();
  popPlot.push(sim.ticks, [lastStats.countA, lastStats.countB]);
  kcPlot.push(sim.ticks, [lastStats.kc]);

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
  // real molecule counts rather than starting out five times too short.
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
let stepsSinceMeasure = 0;
let framesSinceMeasure = 0;
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
      stepsSinceMeasure++;
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

  framesSinceMeasure++;
  const elapsed = now - measureStart;
  if (elapsed >= 500) {
    const sps = (stepsSinceMeasure * 1000) / elapsed;
    perfSteps.textContent = `${Math.round(sps).toLocaleString()} steps/s`;
    perfFps.textContent = `${Math.round((framesSinceMeasure * 1000) / elapsed)} fps`;
    perfParticles.textContent = `${sim.count.toLocaleString()} particles`;
    stepsSinceMeasure = 0;
    framesSinceMeasure = 0;
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
