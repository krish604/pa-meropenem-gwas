// Formatting. One rule everywhere: an absent measurement renders as an
// explicit dash or the words "not reported", never as 0 and never as a
// plausible-looking substitute. A metric that reads 0% but means "unknown"
// is worse than no metric, because it is believed.

export const EM_DASH = '—';

/** Seconds as H:MM:SS, or the dash when there is nothing to show. */
export function duration(seconds) {
  if (seconds === null || seconds === undefined || !isFinite(seconds)) return EM_DASH;
  const s = Math.max(0, Math.floor(seconds));
  const h = Math.floor(s / 3600);
  const m = Math.floor((s % 3600) / 60);
  const sec = s % 60;
  const pad = (n) => String(n).padStart(2, '0');
  return h > 0 ? `${h}:${pad(m)}:${pad(sec)}` : `${pad(m)}:${pad(sec)}`;
}

/** Seconds as "MM:SS (Ns)" — the detail panel's elapsed style. */
export function durationLong(seconds) {
  if (seconds === null || seconds === undefined || !isFinite(seconds)) return EM_DASH;
  return `${duration(seconds)} (${Math.round(seconds)} s)`;
}

export function minutes(value) {
  if (value === null || value === undefined || !isFinite(value)) return EM_DASH;
  if (value === 0) return '0 min';
  if (value < 1) return '<1 min';
  if (value < 60) return `${Math.round(value)} min`;
  const h = Math.floor(value / 60);
  const m = Math.round(value % 60);
  return m ? `${h} h ${m} min` : `${h} h`;
}

export function rate(value, unit) {
  if (value === null || value === undefined || !isFinite(value)) return EM_DASH;
  if (value === 0) return `0 ${unit}`;
  return `${value < 0.01 ? value.toFixed(4) : value.toFixed(2)} ${unit}`;
}

export function percent(value, digits = 0) {
  if (value === null || value === undefined || !isFinite(value)) return EM_DASH;
  return `${value.toFixed(digits)}%`;
}

export function gb(mb) {
  if (mb === null || mb === undefined || !isFinite(mb)) return EM_DASH;
  return `${(mb / 1024).toFixed(1)} GB`;
}

export function mb(mbValue) {
  if (mbValue === null || mbValue === undefined || !isFinite(mbValue)) return EM_DASH;
  return `${mbValue.toFixed(1)} MB`;
}

export function bytesPerSecond(value) {
  if (value === null || value === undefined || !isFinite(value)) return EM_DASH;
  const units = ['B/s', 'KB/s', 'MB/s', 'GB/s'];
  let v = value;
  let i = 0;
  while (v >= 1024 && i < units.length - 1) { v /= 1024; i += 1; }
  return `${v < 10 && i > 0 ? v.toFixed(1) : Math.round(v)} ${units[i]}`;
}

export function count(value) {
  return value === null || value === undefined ? EM_DASH : String(value);
}

/** An ISO stamp from the store as local HH:MM:SS. */
export function clockTime(stamp) {
  if (!stamp) return EM_DASH;
  const d = new Date(stamp);
  if (isNaN(d.getTime())) return String(stamp);
  return d.toLocaleTimeString('en-GB', { hour12: false });
}

export function localStamp(epochSeconds) {
  const d = new Date((epochSeconds || 0) * 1000);
  return d.toLocaleTimeString('en-GB', { hour12: false });
}

/** Text nodes only: never interpolate untrusted strings as HTML. */
export function el(tag, className, text) {
  const node = document.createElement(tag);
  if (className) node.className = className;
  if (text !== undefined && text !== null) node.textContent = String(text);
  return node;
}

export function clear(node) {
  while (node.firstChild) node.removeChild(node.firstChild);
  return node;
}

/** Mark a cell as genuinely unavailable, so styling can grey it. */
export function na(value) {
  return value === null || value === undefined || value === '' || value === EM_DASH;
}
