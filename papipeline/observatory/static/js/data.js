// Talking to the server, and holding the live connection.
//
// The store is authoritative, so a snapshot always wins over an event: an
// event says something changed, and the snapshot says what it is now. On
// reconnect the client re-reads rather than assuming its old picture is
// still true, because a pipeline that ran for an hour while the tab was
// closed has moved on.

import { byId } from './util.js';

const RECONNECT_BASE_MS = 700;
const RECONNECT_MAX_MS = 15000;

export async function getJSON(path) {
  const res = await fetch(path, { headers: { Accept: 'application/json' } });
  if (!res.ok) throw new Error(`${res.status} ${res.statusText} for ${path}`);
  return res.json();
}

export class Stream extends EventTarget {
  constructor({ runKey, onSnapshot, onEvent, onStatus, onHello }) {
    super();
    this.runKey = runKey;
    this.onSnapshot = onSnapshot || (() => {});
    this.onEvent = onEvent || (() => {});
    this.onStatus = onStatus || (() => {});
    this.onHello = onHello || (() => {});
    this.source = null;
    this.attempt = 0;
    this.sequence = 0;
    this.closed = false;
    this.lastMessageAt = 0;
  }

  connect() {
    this.closed = false;
    this._open();
  }

  _open() {
    if (this.closed) return;
    const url = `/api/stream?since=${this.sequence}`;
    // EventSource reconnects on its own, but it will not tell us why, and
    // we want visible status plus a backoff we control.
    this.source = new EventSource(url);

    this.source.addEventListener('open', () => {
      this.attempt = 0;
      this.onStatus('live');
    });

    this.source.addEventListener('hello', (e) => {
      const data = JSON.parse(e.data);
      this.lastMessageAt = Date.now();
      if (typeof data.sequence === 'number') this.sequence = data.sequence;
      this.onHello(data);
      this.onStatus('live');
    });

    this.source.addEventListener('snapshot', (e) => {
      this.lastMessageAt = Date.now();
      this.onSnapshot(JSON.parse(e.data));
    });

    this.source.addEventListener('event', (e) => {
      this.lastMessageAt = Date.now();
      const data = JSON.parse(e.data);
      this.onEvent(data);
    });

    this.source.addEventListener('error', () => {
      this.onStatus('down');
      // EventSource retries on its own; close it so our backoff governs and
      // we can also refresh state rather than only resuming the stream.
      try { this.source.close(); } catch { /* already closed */ }
      if (this.closed) return;
      const delay = Math.min(RECONNECT_BASE_MS * 2 ** this.attempt, RECONNECT_MAX_MS);
      this.attempt += 1;
      setTimeout(() => {
        // State may have moved on while we were disconnected.
        this.refresh().finally(() => this._open());
      }, delay);
    });
  }

  /** Re-read authoritative state. Used on connect and on reconnect. */
  async refresh() {
    try {
      const snapshot = await getJSON(`/api/snapshot?run=${encodeURIComponent(this.runKey)}`);
      this.onSnapshot(snapshot);
      return snapshot;
    } catch (err) {
      this.onStatus('down', String(err));
      return null;
    }
  }

  close() {
    this.closed = true;
    if (this.source) {
      try { this.source.close(); } catch { /* already closed */ }
      this.source = null;
    }
  }
}

/** The clock in the top bar. Local time, labelled as such. */
export function startClock() {
  const node = byId('clock');
  if (!node) return () => {};
  const tick = () => {
    const now = new Date();
    const stamp = now.toLocaleString('en-GB', { hour12: false }).replace(',', '');
    node.textContent = stamp;
  };
  tick();
  const handle = setInterval(tick, 1000);
  return () => clearInterval(handle);
}

/** Fallback stage list. The server's /api/structure is authoritative; this
 *  only keeps the canvas renderable before the first response lands. */
export const STAGE_LIST = [
  'validation', 'annotation', 'mlst', 'amr', 'regulators',
  'structural_variants', 'mechanisms', 'virulence', 'pangenome',
  'phylogeny', 'phenotype', 'gwas', 'convergence', 'cooccurrence',
  'integration', 'reporting',
];

export async function getTask(stage, subject) {
  const path = `/api/tasks/${encodeURIComponent(stage)}/${encodeURIComponent(subject || '')}`;
  return getJSON(path);
}

export async function getLog(stage, subject) {
  const path = `/api/logs/${encodeURIComponent(stage)}/${encodeURIComponent(subject || '')}`;
  return getJSON(path);
}
