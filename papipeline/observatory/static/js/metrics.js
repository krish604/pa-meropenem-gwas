// The bottom metrics zone.
//
// Everything here is measured, not assumed. The sparklines plot a real
// rolling window of host samples; a metric that has not been measurable
// twice yet shows a dash and an empty plot rather than a flat line at
// zero, because a flat line at zero is a claim and this interface does not
// make claims it cannot support.

import {
  bytesPerSecond, duration, gb, minutes, na, percent, rate, EM_DASH,
} from './format.js';
import { byId, Series } from './util.js';

const WINDOW = 60;

const series = {
  cpu: new Series(WINDOW),
  ram: new Series(WINDOW),
  disk: new Series(WINDOW),
  uptime: new Series(WINDOW),
};

let lastSnapshot = null;

export function update(snapshot) {
  lastSnapshot = snapshot;
  const host = snapshot.host || {};
  const tput = snapshot.throughput || {};
  const counts = snapshot.counts || {};

  // ── host series ────────────────────────────────────────────────
  series.cpu.push(host.cpu_percent ?? null);
  series.ram.push(host.memory_percent ?? null);
  series.disk.push((host.disk_read_bps ?? 0) + (host.disk_write_bps ?? 0) || null);
  series.uptime.push(host.uptime_seconds ?? null);

  drawSpark('cpu', series.cpu, '#22c55e');
  drawSpark('ram', series.ram, '#3b82f6');
  drawSpark('disk', series.disk, '#f59e0b');
  drawSpark('uptime', series.uptime, '#a855f7');

  setMetric('s-cpu', percent(host.cpu_percent), host.cpu_percent);
  setMetric('s-ram', host.memory_used_mb !== null && host.memory_used_mb !== undefined
    ? `${gb(host.memory_used_mb)} / ${gb(host.memory_total_mb)}` : EM_DASH,
  host.memory_percent);

  const io = (host.disk_read_bps ?? null) || (host.disk_write_bps ?? null);
  setMetric('s-disk', io === null ? EM_DASH : bytesPerSecond(io), io);
  setMetric('s-uptime', duration(host.uptime_seconds), host.uptime_seconds);

  // ── counters ───────────────────────────────────────────────────
  setText('c-completed', counts.SUCCEEDED || 0);
  setText('c-running', counts.RUNNING || 0);
  setText('c-queued', queued(snapshot));
  setText('c-failed', counts.FAILED || 0);
  setText('c-invalid', counts.INVALID || 0);
  setText('c-incomplete', counts.INCOMPLETE || 0);

  // ETA and throughput come with the basis they were derived from; if there
  // is no basis, there is no number.
  setText('c-eta', tput.eta_minutes === null || tput.eta_minutes === undefined
    ? EM_DASH : minutes(tput.eta_minutes));
  setText('c-throughput', tput.tasks_per_minute === null || tput.tasks_per_minute === undefined
    ? EM_DASH : rate(tput.tasks_per_minute, 'tasks/min'));

  // ── cohort progress ────────────────────────────────────────────
  const total = snapshot.total_expected || 0;
  const done = (counts.SUCCEEDED || 0);
  const failedish = (counts.FAILED || 0) + (counts.INVALID || 0) + (counts.INCOMPLETE || 0);
  const finished = done + failedish;
  const pct = total > 0 ? Math.round((finished / total) * 100) : 0;
  const fill = byId('cohort-fill');
  if (fill) fill.style.width = `${total > 0 ? Math.min(100, pct) : 0}%`;
  setText('cohort-fraction', total > 0 ? `${finished} / ${total} genomes (${pct}%)` : '—');

  const note = byId('cohort-note');
  if (note) {
    if (total === 0) {
      note.textContent = 'no genomes recorded in the execution store';
    } else {
      const basis = tput.basis ? ` · ${tput.basis}` : '';
      note.textContent = `${done} validated success${done === 1 ? '' : 'es'}, `
        + `${failedish} needing attention${basis}`;
    }
  }
}

function queued(snapshot) {
  // A task is queued when the store has a row for it that has not started.
  // Absent rows are not counted: a stage that has never been planned is
  // not a queued task.
  return (snapshot.counts || {}).PENDING || 0;
}

function setText(id, value) {
  const node = byId(id);
  if (node) node.textContent = value === null || value === undefined ? EM_DASH : String(value);
}

function setMetric(id, text, sample) {
  const node = byId(id);
  if (!node) return;
  if (na(text) || text === EM_DASH) {
    node.textContent = EM_DASH;
    node.classList.add('na');
    node.title = 'not measurable on this host right now';
  } else {
    node.textContent = text;
    node.classList.remove('na');
    node.title = '';
  }
}

/** A sparkline over the rolling window. No data means no line. */
function drawSpark(name, data, color) {
  const canvas = document.querySelector(`.spark-canvas[data-series="${name}"]`);
  if (!canvas) return;
  const ratio = window.devicePixelRatio || 1;
  const rect = canvas.getBoundingClientRect();
  const w = Math.max(1, rect.width);
  const h = Math.max(1, rect.height);
  if (canvas.width !== Math.round(w * ratio) || canvas.height !== Math.round(h * ratio)) {
    canvas.width = Math.round(w * ratio);
    canvas.height = Math.round(h * ratio);
  }
  const ctx = canvas.getContext('2d');
  ctx.setTransform(ratio, 0, 0, ratio, 0, 0);
  ctx.clearRect(0, 0, w, h);

  const values = data.present;
  if (values.length < 2) return;   // one sample is not a trend

  const min = Math.min(...values);
  const max = Math.max(...values);
  const span = (max - min) || 1;
  const stepX = w / Math.max(1, data.values.length - 1);
  const yFor = (v) => h - ((v - min) / span) * (h - 3) - 1.5;

  // Area under the line, additive so it reads as light rather than paint.
  ctx.save();
  ctx.globalCompositeOperation = 'lighter';
  ctx.beginPath();
  let started = false;
  data.values.forEach((v, i) => {
    if (v === null || v === undefined || !isFinite(v)) return;
    const x = i * stepX;
    const y = yFor(v);
    if (!started) { ctx.moveTo(x, h); ctx.lineTo(x, y); started = true; }
    else ctx.lineTo(x, y);
  });
  if (started) {
    ctx.lineTo(w, h);
    ctx.closePath();
    const grad = ctx.createLinearGradient(0, 0, 0, h);
    grad.addColorStop(0, `${color}44`);
    grad.addColorStop(1, `${color}00`);
    ctx.fillStyle = grad;
    ctx.fill();
  }
  ctx.restore();

  ctx.beginPath();
  started = false;
  data.values.forEach((v, i) => {
    if (v === null || v === undefined || !isFinite(v)) return;
    const x = i * stepX;
    const y = yFor(v);
    if (!started) { ctx.moveTo(x, y); started = true; } else ctx.lineTo(x, y);
  });
  ctx.strokeStyle = color;
  ctx.lineWidth = 1.4;
  ctx.stroke();

  // Mark the most recent real sample.
  const lastX = (data.values.length - 1) * stepX;
  const lastY = yFor(data.last);
  ctx.fillStyle = color;
  ctx.beginPath();
  ctx.arc(lastX, lastY, 1.7, 0, Math.PI * 2);
  ctx.fill();
}

export function metricsSnapshot() { return lastSnapshot; }
