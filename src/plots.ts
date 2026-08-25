/**
 * Small canvas plotting helpers: a scrolling multi-pen line plot and a
 * two-species histogram. Both are deliberately plain -- no library, no DOM
 * per data point -- because they redraw next to a particle simulation that
 * wants the rest of the frame budget.
 */

const AXIS = '#3a4254';
const GRID = '#232936';
const TEXT = '#8b95a9';
const LABEL = '#aab3c5';
const FONT = '11px ui-sans-serif, system-ui, -apple-system, sans-serif';
const SUB_FONT = '8px ui-sans-serif, system-ui, -apple-system, sans-serif';

/** Left edge of the rotated y-axis caption; see drawFrame(). */
const Y_CAPTION_X = 1;

/** A "nice" tick step (1, 2 or 5 times a power of ten) near range/targetTicks. */
function niceStep(range: number, targetTicks: number): number {
  if (!(range > 0)) {
    return 1;
  }
  const raw = range / targetTicks;
  const mag = Math.pow(10, Math.floor(Math.log10(raw)));
  const norm = raw / mag;
  const step = norm <= 1 ? 1 : norm <= 2 ? 2 : norm <= 5 ? 5 : 10;
  return step * mag;
}

function formatTick(v: number, step: number): string {
  if (v === 0) {
    return '0';
  }
  const abs = Math.abs(v);
  // Only genuinely tiny values get exponent form: "2e-4" beats "0.0002" at
  // axis size, but "1e+4" is worse than "10,000" for a tick count or a
  // molecule count, which is all a large value on these axes ever is.
  if (abs < 0.001) {
    return v.toExponential(0);
  }
  const decimals = Math.max(0, -Math.floor(Math.log10(step)));
  return v.toLocaleString(undefined, {
    minimumFractionDigits: 0,
    maximumFractionDigits: Math.min(decimals, 4),
  });
}

interface PenSpec {
  label: string;
  /** Drawn small and low against the label, for names like K_c. */
  sub?: string;
  color: string;
}

/** Shared canvas plumbing: device-pixel sizing, plot frame, axis labels. */
abstract class PlotBase {
  protected canvas: HTMLCanvasElement;
  protected ctx: CanvasRenderingContext2D;
  protected width = 0;
  protected height = 0;
  /**
   * Left gutter: the rotated caption's ~13px of thickness at the very edge,
   * then room for the widest y tick label, which ends 6px short of the plot.
   * Sized for five characters ("4,000", "1,000") now that large values are
   * spelled out rather than given as exponents.
   */
  protected padL = 54;
  protected padR = 10;
  protected padT = 20;
  protected padB = 30;
  protected xLabel: string;
  protected yLabel: string;
  protected pens: PenSpec[];

  constructor(canvas: HTMLCanvasElement, xLabel: string, yLabel: string, pens: PenSpec[]) {
    this.canvas = canvas;
    const ctx = canvas.getContext('2d');
    if (!ctx) {
      throw new Error('2d context unavailable');
    }
    this.ctx = ctx;
    this.xLabel = xLabel;
    this.yLabel = yLabel;
    this.pens = pens;
    this.resize();
  }

  /** Matches the backing store to the element's CSS size and the DPR. */
  resize(): void {
    const dpr = window.devicePixelRatio || 1;
    const rect = this.canvas.getBoundingClientRect();
    const w = Math.max(1, Math.round(rect.width));
    const h = Math.max(1, Math.round(rect.height));
    this.canvas.width = Math.round(w * dpr);
    this.canvas.height = Math.round(h * dpr);
    this.ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    this.width = w;
    this.height = h;
  }

  protected get plotW(): number {
    return this.width - this.padL - this.padR;
  }

  protected get plotH(): number {
    return this.height - this.padT - this.padB;
  }

  protected clearCanvas(): void {
    this.ctx.clearRect(0, 0, this.width, this.height);
  }

  /**
   * Draws grid lines, tick labels and the two axis captions, and returns the
   * value-to-pixel mapping the caller should use for its data.
   */
  protected drawFrame(x0: number, x1: number, y0: number, y1: number): void {
    const { ctx } = this;
    const L = this.padL;
    const T = this.padT;
    const W = this.plotW;
    const H = this.plotH;

    ctx.font = FONT;
    ctx.lineWidth = 1;

    // Horizontal grid + y tick labels.
    const yStep = niceStep(y1 - y0, 4);
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    for (let v = Math.ceil(y0 / yStep) * yStep; v <= y1 + yStep * 1e-9; v += yStep) {
      const py = Math.round(T + H - ((v - y0) / (y1 - y0)) * H) + 0.5;
      ctx.strokeStyle = GRID;
      ctx.beginPath();
      ctx.moveTo(L, py);
      ctx.lineTo(L + W, py);
      ctx.stroke();
      ctx.fillStyle = TEXT;
      ctx.fillText(formatTick(v, yStep), L - 6, py);
    }

    // Vertical grid + x tick labels.
    const xStep = niceStep(x1 - x0, 5);
    ctx.textAlign = 'center';
    ctx.textBaseline = 'top';
    for (let v = Math.ceil(x0 / xStep) * xStep; v <= x1 + xStep * 1e-9; v += xStep) {
      const px = Math.round(L + ((v - x0) / (x1 - x0)) * W) + 0.5;
      ctx.strokeStyle = GRID;
      ctx.beginPath();
      ctx.moveTo(px, T);
      ctx.lineTo(px, T + H);
      ctx.stroke();
      ctx.fillStyle = TEXT;
      ctx.fillText(formatTick(v, xStep), px, T + H + 6);
    }

    // Axis lines.
    ctx.strokeStyle = AXIS;
    ctx.beginPath();
    ctx.moveTo(L + 0.5, T);
    ctx.lineTo(L + 0.5, T + H + 0.5);
    ctx.lineTo(L + W, T + H + 0.5);
    ctx.stroke();

    // Captions.
    ctx.fillStyle = LABEL;
    ctx.textAlign = 'right';
    ctx.textBaseline = 'bottom';
    ctx.fillText(this.xLabel, L + W, this.height - 2);
    ctx.save();
    // The rotated caption's *thickness* runs rightward from this x, so it
    // occupies roughly [Y_CAPTION_X, Y_CAPTION_X + line height] and the y
    // tick labels end at L - 6. Kept hard against the left edge to leave
    // that gutter for the numbers, which are the wider of the two.
    ctx.translate(Y_CAPTION_X, T);
    ctx.rotate(-Math.PI / 2);
    ctx.textAlign = 'right';
    ctx.textBaseline = 'top';
    ctx.fillText(this.yLabel, 0, 0);
    ctx.restore();
  }

  /** A compact colored-square legend along the top of the plot area. */
  protected drawLegend(): void {
    const { ctx } = this;
    ctx.font = FONT;
    ctx.textBaseline = 'middle';
    ctx.textAlign = 'left';
    let x = this.padL + 4;
    const y = 10;
    for (const pen of this.pens) {
      ctx.fillStyle = pen.color;
      ctx.fillRect(x, y - 4, 8, 8);
      ctx.fillStyle = LABEL;
      ctx.fillText(pen.label, x + 12, y);
      let width = ctx.measureText(pen.label).width;
      if (pen.sub) {
        ctx.font = SUB_FONT;
        ctx.fillText(pen.sub, x + 12 + width, y + 3);
        width += ctx.measureText(pen.sub).width;
        ctx.font = FONT;
      }
      x += 12 + width + 14;
    }
  }
}

/**
 * A scrolling line plot over a fixed-capacity ring buffer. Old samples fall
 * off the left rather than accumulating forever, which also keeps the y axis
 * autoscaled to what is currently on screen -- the reason a big opening
 * transient in Kc stops squashing the plateau once it scrolls away.
 */
export class TimeSeriesPlot extends PlotBase {
  private capacity: number;
  private ts: Float64Array;
  private values: Float64Array[];
  private head = 0;
  private length = 0;
  /** Forces y = 0 into range, so counts read against a real baseline. */
  private includeZero: boolean;

  constructor(
    canvas: HTMLCanvasElement,
    opts: { xLabel: string; yLabel: string; pens: PenSpec[]; capacity?: number; includeZero?: boolean },
  ) {
    super(canvas, opts.xLabel, opts.yLabel, opts.pens);
    this.capacity = opts.capacity ?? 1200;
    this.includeZero = opts.includeZero ?? true;
    this.ts = new Float64Array(this.capacity);
    this.values = opts.pens.map(() => new Float64Array(this.capacity));
  }

  clear(): void {
    this.head = 0;
    this.length = 0;
  }

  push(t: number, values: number[]): void {
    const i = this.head;
    this.ts[i] = t;
    for (let p = 0; p < this.values.length; p++) {
      this.values[p][i] = values[p];
    }
    this.head = (i + 1) % this.capacity;
    if (this.length < this.capacity) {
      this.length++;
    }
  }

  /** Ring index of the k-th oldest visible sample. */
  private at(k: number): number {
    return (this.head - this.length + k + this.capacity * 2) % this.capacity;
  }

  draw(): void {
    this.clearCanvas();
    if (this.length < 1) {
      this.drawFrame(0, 1, 0, 1);
      this.drawLegend();
      return;
    }

    const x0 = this.ts[this.at(0)];
    const x1 = Math.max(this.ts[this.at(this.length - 1)], x0 + 1e-9);

    let lo = this.includeZero ? 0 : Infinity;
    let hi = -Infinity;
    for (const series of this.values) {
      for (let k = 0; k < this.length; k++) {
        const v = series[this.at(k)];
        if (v < lo) {
          lo = v;
        }
        if (v > hi) {
          hi = v;
        }
      }
    }
    if (!Number.isFinite(lo)) {
      lo = 0;
    }
    if (!Number.isFinite(hi)) {
      hi = 1;
    }
    if (hi <= lo) {
      hi = lo + 1;
    }
    // A little headroom so a pen riding the maximum is not clipped to the frame.
    hi += (hi - lo) * 0.08;

    this.drawFrame(x0, x1, lo, hi);

    const { ctx } = this;
    const L = this.padL;
    const T = this.padT;
    const W = this.plotW;
    const H = this.plotH;
    const sx = W / (x1 - x0);
    const sy = H / (hi - lo);

    ctx.save();
    ctx.beginPath();
    ctx.rect(L, T, W, H);
    ctx.clip();
    ctx.lineWidth = 1.5;
    ctx.lineJoin = 'round';
    for (let p = 0; p < this.values.length; p++) {
      const series = this.values[p];
      ctx.strokeStyle = this.pens[p].color;
      ctx.beginPath();
      for (let k = 0; k < this.length; k++) {
        const idx = this.at(k);
        const px = L + (this.ts[idx] - x0) * sx;
        const py = T + H - (series[idx] - lo) * sy;
        if (k === 0) {
          ctx.moveTo(px, py);
        } else {
          ctx.lineTo(px, py);
        }
      }
      ctx.stroke();
    }
    ctx.restore();

    this.drawLegend();
  }
}

/**
 * Two overlaid translucent bar series on a shared bin grid -- the speed
 * distribution, one pen per species.
 */
export class HistogramPlot extends PlotBase {
  private bins: number;
  private data: Float64Array[];
  private hiSpeed = 1;

  constructor(canvas: HTMLCanvasElement, opts: { xLabel: string; yLabel: string; pens: PenSpec[]; bins: number }) {
    super(canvas, opts.xLabel, opts.yLabel, opts.pens);
    this.bins = opts.bins;
    this.data = opts.pens.map(() => new Float64Array(opts.bins));
  }

  get binArrays(): Float64Array[] {
    return this.data;
  }

  setRange(hiSpeed: number): void {
    this.hiSpeed = Math.max(hiSpeed, 1e-6);
  }

  clear(): void {
    for (const d of this.data) {
      d.fill(0);
    }
    this.hiSpeed = 1;
  }

  draw(): void {
    this.clearCanvas();
    let hi = 0;
    for (const d of this.data) {
      for (let i = 0; i < d.length; i++) {
        if (d[i] > hi) {
          hi = d[i];
        }
      }
    }
    if (hi <= 0) {
      hi = 1;
    }
    hi *= 1.1;

    this.drawFrame(0, this.hiSpeed, 0, hi);

    const { ctx } = this;
    const L = this.padL;
    const T = this.padT;
    const W = this.plotW;
    const H = this.plotH;
    const barW = W / this.bins;

    ctx.save();
    ctx.beginPath();
    ctx.rect(L, T, W, H);
    ctx.clip();
    for (let p = 0; p < this.data.length; p++) {
      const d = this.data[p];
      ctx.fillStyle = this.pens[p].color;
      for (let i = 0; i < this.bins; i++) {
        const v = d[i];
        if (v <= 0) {
          continue;
        }
        const h = (v / hi) * H;
        ctx.fillRect(L + i * barW, T + H - h, Math.max(1, barW - 0.5), h);
      }
    }
    ctx.restore();

    this.drawLegend();
  }
}
