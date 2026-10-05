/* A minimal DOM shim for the pure-JS unit tests.
 *
 * No jsdom, no downloads: just enough of the DOM for `h()` / `createRouter()`
 * / the badge and banner builders. Every node records its children and text so
 * a test can assert what would have been rendered.
 */

class ShimNode {}

class ShimText extends ShimNode {
  constructor(text) {
    super();
    this.nodeType = 3;
    this._text = String(text);
  }
  get textContent() {
    return this._text;
  }
  set textContent(value) {
    this._text = String(value);
  }
}

class ShimElement extends ShimNode {
  constructor(tag) {
    super();
    this.nodeType = 1;
    this.tagName = String(tag).toUpperCase();
    this.children = [];
    this.attributes = {};
    this.dataset = {};
    this.style = { setProperty() {} };
    this.className = '';
    this.id = '';
    this._text = '';
    this._listeners = new Map();
  }

  get childNodes() {
    return this.children;
  }

  get firstChild() {
    return this.children.length ? this.children[0] : null;
  }

  get textContent() {
    if (this.children.length === 0) return this._text;
    return this.children.map((child) => child.textContent).join('');
  }

  set textContent(value) {
    this.children = [];
    this._text = String(value);
  }

  appendChild(child) {
    this.children.push(child);
    child.parentNode = this;
    return child;
  }

  removeChild(child) {
    const index = this.children.indexOf(child);
    if (index >= 0) this.children.splice(index, 1);
    return child;
  }

  append(...nodes) {
    for (const node of nodes) this.appendChild(node);
  }

  setAttribute(name, value) {
    this.attributes[name] = String(value);
  }

  getAttribute(name) {
    return this.attributes[name];
  }

  addEventListener(type, fn) {
    if (!this._listeners.has(type)) this._listeners.set(type, new Set());
    this._listeners.get(type).add(fn);
  }

  removeEventListener(type, fn) {
    const set = this._listeners.get(type);
    if (set) set.delete(fn);
  }

  dispatch(type, event) {
    const set = this._listeners.get(type);
    if (set) for (const fn of set) fn(event);
  }

  /* -- APIs the canvas page modules need, added additively ------------- */

  get clientWidth() {
    return this._clientWidth ?? 800;
  }

  get clientHeight() {
    return this._clientHeight ?? 600;
  }

  getBoundingClientRect() {
    return {
      width: this.clientWidth,
      height: this.clientHeight,
      top: 0,
      left: 0,
      right: this.clientWidth,
      bottom: this.clientHeight,
      x: 0,
      y: 0,
    };
  }

  getContext() {
    return canvasContext();
  }

  querySelector() {
    return null;
  }

  querySelectorAll() {
    return [];
  }

  closest() {
    return null;
  }

  remove() {
    if (this.parentNode && this.parentNode.removeChild) {
      this.parentNode.removeChild(this);
    }
  }

  click() {
    this.dispatch('click', { type: 'click' });
  }

  focus() {}

  get classList() {
    const self = this;
    const classes = () => String(self.className || '').split(/\s+/).filter(Boolean);
    return {
      add: (c) => { self.className = [...new Set([...classes(), c])].join(' '); },
      remove: (c) => { self.className = classes().filter((x) => x !== c).join(' '); },
      contains: (c) => classes().includes(c),
      toggle: (c) => {
        if (classes().includes(c)) self.className = classes().filter((x) => x !== c).join(' ');
        else self.className = [...classes(), c].join(' ');
      },
    };
  }
}

/** A no-op 2D context: every method is a function, every property is settable. */
function canvasContext() {
  const store = {};
  return new Proxy(store, {
    get(target, prop) {
      if (prop in target) return target[prop];
      return () => {};
    },
    set(target, prop, value) {
      target[prop] = value;
      return true;
    },
  });
}

const documentElement = new ShimElement('html');

const documentShim = {
  documentElement,
  body: new ShimElement('body'),
  createElement: (tag) => new ShimElement(tag),
  createElementNS: (_ns, tag) => new ShimElement(tag),
  createTextNode: (text) => new ShimText(text),
  getElementById: () => null,
  querySelector: () => null,
  querySelectorAll: () => [],
  addEventListener() {},
  removeEventListener() {},
};

const locationShim = { hash: '', pathname: '/', search: '' };
const historyShim = { replaceState() {} };

globalThis.Node = ShimNode;
globalThis.document = documentShim;
globalThis.location = locationShim;
globalThis.history = historyShim;
globalThis.window = globalThis;
globalThis.devicePixelRatio = 1;
globalThis.matchMedia = () => ({
  matches: false,
  media: '',
  addEventListener() {},
  removeEventListener() {},
  addListener() {},
  removeListener() {},
});
globalThis.getComputedStyle = () => ({
  getPropertyValue: () => '',
});
// A no-op rAF: the mount-contract test asserts the *synchronous* setup and
// teardown, and a page whose data never arrives must not schedule a draw that
// would run against a null tree after the test has finished.
let _rafId = 0;
globalThis.requestAnimationFrame = () => {
  _rafId += 1;
  return _rafId;
};
globalThis.cancelAnimationFrame = () => {};

class ShimObserver {
  observe() {}
  unobserve() {}
  disconnect() {}
  takeRecords() {
    return [];
  }
}
globalThis.ResizeObserver = globalThis.ResizeObserver || ShimObserver;
globalThis.MutationObserver = globalThis.MutationObserver || ShimObserver;
globalThis.IntersectionObserver = globalThis.IntersectionObserver || ShimObserver;

class ShimEventSource {
  constructor() {
    this.readyState = 0;
  }
  addEventListener() {}
  removeEventListener() {}
  close() {}
}
globalThis.EventSource = globalThis.EventSource || ShimEventSource;

if (typeof globalThis.addEventListener !== 'function') {
  const listeners = new Map();
  globalThis.addEventListener = (type, fn) => {
    if (!listeners.has(type)) listeners.set(type, new Set());
    listeners.get(type).add(fn);
  };
  globalThis.removeEventListener = (type, fn) => {
    const set = listeners.get(type);
    if (set) set.delete(fn);
  };
}

export { ShimElement, ShimText, documentShim, locationShim, canvasContext };
