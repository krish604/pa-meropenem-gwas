// Small helpers with no dependencies.

export function clamp(v, lo, hi) {
  return Math.min(hi, Math.max(lo, v));
}

export function byId(id) {
  return document.getElementById(id);
}

export function on(el, event, handler) {
  if (el) el.addEventListener(event, handler);
}

/** Replace children with a fragment, in one reflow. */
export function replace(node, children) {
  if (!node) return node;
  const frag = document.createDocumentFragment();
  for (const child of children) {
    frag.appendChild(typeof child === 'string' ? document.createTextNode(child) : child);
  }
  node.textContent = '';
  node.appendChild(frag);
  return node;
}

/**
 * A fixed-length ring of samples for the sparklines.
 *
 * Fixed length matters: a growing series would rescale its own axis, so a
 * flat line and a spiking one would look identical, and an empty series
 * would be indistinguishable from a series of zeroes.
 */
export class Series {
  constructor(length = 60) {
    this.length = length;
    this.values = [];
  }

  push(value) {
    this.values.push(value);
    if (this.values.length > this.length) this.values.shift();
  }

  /** Values with nulls removed; nulls mean "not measurable", not zero. */
  get present() {
    return this.values.filter((v) => v !== null && v !== undefined && isFinite(v));
  }

  get last() {
    for (let i = this.values.length - 1; i >= 0; i -= 1) {
      const v = this.values[i];
      if (v !== null && v !== undefined && isFinite(v)) return v;
    }
    return null;
  }
}
