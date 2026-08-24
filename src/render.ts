/**
 * Draws the gas. Everything is batched: one path per species, filled once,
 * rather than a fill call per particle -- the difference between a few
 * hundred and a few thousand particles staying at 60 fps.
 */

import { MAX_PXCOR, SPECIES_A, SPECIES_B, type Sim } from './sim';

export const COLOR_A = '#57b6e8';
export const COLOR_B = '#e05561';
const COLOR_WALL = '#e2c94a';
const COLOR_FLOOR = '#0b0d12';

/** World span in patches, matching resize_world(-40, 40, -40, 40). */
const WORLD_SPAN = 2 * MAX_PXCOR + 1;

export class Renderer {
  private canvas: HTMLCanvasElement;
  private ctx: CanvasRenderingContext2D;
  private size = 0;

  constructor(canvas: HTMLCanvasElement) {
    this.canvas = canvas;
    const ctx = canvas.getContext('2d', { alpha: false });
    if (!ctx) {
      throw new Error('2d context unavailable');
    }
    this.ctx = ctx;
    this.resize();
  }

  resize(): void {
    const dpr = window.devicePixelRatio || 1;
    const rect = this.canvas.getBoundingClientRect();
    const size = Math.max(1, Math.round(Math.min(rect.width, rect.height)));
    this.canvas.width = Math.round(size * dpr);
    this.canvas.height = Math.round(size * dpr);
    this.ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    this.size = size;
  }

  draw(sim: Sim): void {
    const { ctx } = this;
    const size = this.size;
    // Patches per pixel. y is flipped: NetLogo's +y is up, the canvas's is down.
    const scale = size / WORLD_SPAN;
    const half = size / 2;

    ctx.fillStyle = COLOR_FLOOR;
    ctx.fillRect(0, 0, size, size);

    // The box wall: one patch wide, centered on |coord| == boxEdge, which is
    // exactly what stroking a rect of that half-width with lineWidth = one
    // patch produces.
    ctx.strokeStyle = COLOR_WALL;
    ctx.lineWidth = scale;
    const edge = sim.boxEdge * scale;
    ctx.strokeRect(half - edge, half - edge, 2 * edge, 2 * edge);

    const { x, y, radius, species, alive, top } = sim.view();

    for (const target of [SPECIES_A, SPECIES_B]) {
      ctx.fillStyle = target === SPECIES_A ? COLOR_A : COLOR_B;
      // Below about a pixel and a half across an arc is indistinguishable
      // from a rect and several times more expensive to rasterize. Radius is
      // a per-species constant, so the first particle settles the choice.
      let useRects = false;
      for (let i = 0; i < top; i++) {
        if (alive[i] && species[i] === target) {
          useRects = radius[i] * scale < 1.5;
          break;
        }
      }

      if (useRects) {
        for (let i = 0; i < top; i++) {
          if (!alive[i] || species[i] !== target) {
            continue;
          }
          const r = radius[i] * scale;
          ctx.fillRect(half + x[i] * scale - r, half - y[i] * scale - r, 2 * r, 2 * r);
        }
      } else {
        ctx.beginPath();
        for (let i = 0; i < top; i++) {
          if (!alive[i] || species[i] !== target) {
            continue;
          }
          const px = half + x[i] * scale;
          const py = half - y[i] * scale;
          const r = radius[i] * scale;
          ctx.moveTo(px + r, py);
          ctx.arc(px, py, r, 0, Math.PI * 2);
        }
        ctx.fill();
      }
    }
  }
}
