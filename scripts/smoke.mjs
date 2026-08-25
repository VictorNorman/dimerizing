/**
 * Runs the *built* bundle (dist/) against a stub DOM and canvas, to prove the
 * UI wiring works end to end: every element id resolves, the sliders get
 * built, the animation loop steps the gas, the monitors and plots update, the
 * buttons do what they say, and nothing throws.
 *
 *   npm run build && npm run smoke
 *
 * This is deliberately not a browser test -- it cannot tell you the page
 * *looks* right -- but it catches every way the wiring can break without one.
 */
import { readdirSync } from 'node:fs';
import { pathToFileURL } from 'node:url';

const calls = new Map();
const ctxStub = new Proxy(
  {
    canvas: null,
    measureText: () => ({ width: 10 }),
  },
  {
    get(t, k) {
      if (k in t) {
        return t[k];
      }
      if (typeof k === 'string') {
        calls.set(k, (calls.get(k) ?? 0) + 1);
        return () => {};
      }
      return undefined;
    },
    set(t, k, v) { t[k] = v; return true; },
  },
);

const madeNodes = [];
function makeNode(tag = 'div') {
  const node = {
    tagName: tag.toUpperCase(),
    style: {},
    _text: '',
    children: [],
    listeners: {},
    classList: {
      add() {},
      remove() {},
      toggle() {},
    },
    className: '',
    get textContent() { return this._text; },
    set textContent(v) { this._text = String(v); },
    append(...kids) { this.children.push(...kids); },
    addEventListener(ev, fn) { (this.listeners[ev] ||= []).push(fn); },
    getBoundingClientRect: () => ({
      width: 640,
      height: 200,
      top: 0,
      left: 0,
    }),
    getContext: () => ctxStub,
    width: 0,
    height: 0,
  };
  madeNodes.push(node);
  return node;
}

const ids = ['sliders', 'world', 'plot-pop', 'plot-kc', 'plot-speed', 'mon-a', 'mon-b',
  'mon-total', 'mon-temp', 'mon-kc', 'mon-ticks', 'perf-particles',
  'btn-setup', 'btn-go', 'btn-step', 'dot-a', 'dot-b'];
const byId = Object.fromEntries(ids.map((id) => [id, makeNode(id.startsWith('plot') || id === 'world' ? 'canvas' : 'div')]));

let rafCallback = null;
globalThis.MutationObserver = class { observe() {} disconnect() {} };
globalThis.fetch = async () => ({ ok: true });
globalThis.document = {
  getElementById: (id) => byId[id] ?? null,
  createElement: (tag) => {
    const n = makeNode(tag);
    if (tag === 'link') {
      n.relList = { supports: () => true };
    }
    return n;
  },
  querySelectorAll: () => [],
};
globalThis.window = {
  devicePixelRatio: 2,
  addEventListener: () => {},
  setTimeout: () => 0,
  clearTimeout: () => {},
};
globalThis.requestAnimationFrame = (fn) => { rafCallback = fn; return 1; };
globalThis.performance = performance;

const assetDir = new URL('../dist/assets/', import.meta.url);
const bundle = readdirSync(assetDir).find((f) => f.endsWith('.js'));
if (!bundle) {
  throw new Error('no built bundle found -- run `npm run build` first');
}
await import(pathToFileURL(new URL(bundle, assetDir).pathname).href);

let failures = 0;
const check = (name, ok, detail = '') => {
  console.log(`  ${ok ? 'ok  ' : 'FAIL'}  ${name}${detail ? '  ' + detail : ''}`);
  if (!ok) {
    failures++;
  }
};

check('module loaded and requested a frame', rafCallback !== null);
check('slider rows built', byId.sliders.children.length === 8, `${byId.sliders.children.length} rows`);
check('monitors initialised', byId['mon-a'].textContent === '200', `# of A = ${byId['mon-a'].textContent}`);
check('conserved total initialised', byId['mon-total'].textContent === '200');
check('species swatches sized from the real radii',
  byId['dot-a'].style.background === '#57b6e8' && byId['dot-b'].style.background === '#e05561'
  && Math.abs(parseFloat(byId['dot-b'].style.width) / parseFloat(byId['dot-a'].style.width) - Math.cbrt(2)) < 1e-9,
  `A ${byId['dot-a'].style.width}, B ${byId['dot-b'].style.width}`);

// Press Go, then run frames.
byId['btn-go'].listeners.click[0]();
// The perf readout only refreshes on a 500 ms wall-clock window, so drive
// frames until real time has actually passed rather than a fixed count.
const t0 = performance.now();
let frames = 0;
while (performance.now() - t0 < 700 && frames < 100000) {
  const fn = rafCallback;
  rafCallback = null;
  fn();
  frames++;
}

const a = Number(byId['mon-a'].textContent);
const b = Number(byId['mon-b'].textContent);
check('simulation advanced', Number(byId['mon-ticks'].textContent) > 0, `ticks = ${byId['mon-ticks'].textContent}`);
check('dimerization visible in the monitors', b > 0, `A = ${a}, B = ${b}`);
check('A + 2B still 200', a + 2 * b === 200, `= ${a + 2 * b}`);
check('measured temperature is a number', Number.isFinite(Number(byId['mon-temp'].textContent)),
  `T = ${byId['mon-temp'].textContent}`);
check('Kc monitor populated', byId['mon-kc'].textContent !== '0', `Kc = ${byId['mon-kc'].textContent}`);
check('particle count readout populated', byId['perf-particles'].textContent.includes('particles'),
  byId['perf-particles'].textContent);
check('particles were drawn', (calls.get('fill') ?? 0) + (calls.get('fillRect') ?? 0) > 0);
check('plots were drawn', (calls.get('stroke') ?? 0) > 0 && (calls.get('fillText') ?? 0) > 0);

// Setup button must reset cleanly.
byId['btn-setup'].listeners.click[0]();
check('Setup resets the run', byId['mon-ticks'].textContent === '0.0' && byId['mon-total'].textContent === '200',
  `ticks = ${byId['mon-ticks'].textContent}, total = ${byId['mon-total'].textContent}`);

// Step button while paused.
byId['btn-step'].listeners.click[0]();
check('Step advances one step', Number(byId['mon-ticks'].textContent) > 0);

console.log(failures === 0 ? '\nUI smoke test passed\n' : `\n${failures} failure(s)\n`);
process.exit(failures ? 1 : 0);
