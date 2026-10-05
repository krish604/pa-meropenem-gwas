/* router.js — the hash router, and the only thing that mounts a page.
 *
 * A page module is an ES module with exactly one export the router calls:
 *
 *     export function mount(container, ctx) { ...; return teardown }
 *
 * `teardown` is `() => void` and must remove every listener, abort every
 * in-flight AbortController and close every EventSource the module opened
 * (DESIGN §7). **A module that returns nothing is a leak, and this router
 * refuses it**: the teardown call is skipped, the violation is logged loudly,
 * the page is replaced by a panel that says the contract was broken, and the
 * violation is recorded on `window.__paDashboard.contractViolations` so a test
 * can assert it without scraping the console.
 *
 * Page modules are loaded with a dynamic `import()` rather than a static one,
 * and the specifier is built from `import.meta.url`. Two consequences, both
 * wanted: a module another owner has not written yet degrades to a named
 * "module not yet present" panel instead of a blank page, and nothing is ever
 * fetched from an origin other than this one (UI-D4).
 */

/** Thrown when a page module is not on disk, or fails to parse. */
export class ModuleMissingError extends Error {
  constructor(specifier, cause) {
    super(`page module not present: ${specifier}`);
    this.name = 'ModuleMissingError';
    this.specifier = specifier;
    this.cause = cause;
  }
}

/** Thrown when a page module mounts and returns no teardown. */
export class TeardownContractError extends Error {
  constructor(id, returned) {
    super(
      `the page module for "${id}" returned ${describe(returned)} instead of a ` +
        `teardown function. Its listeners, in-flight requests and EventSource ` +
        `connections cannot be released, so the router refuses the mount. ` +
        `DESIGN §7: every page module returns () => void.`
    );
    this.name = 'TeardownContractError';
    this.routeId = id;
    this.returned = returned;
  }
}

function describe(value) {
  if (value === undefined) return 'undefined';
  if (value === null) return 'null';
  if (typeof value === 'object') return `an object (${Object.prototype.toString.call(value)})`;
  return `a ${typeof value}`;
}

/** Where the router records contract violations. Read by QA; never read by a page. */
const VIOLATIONS_KEY = '__paDashboard';

function recordViolation(entry) {
  const root = (globalThis[VIOLATIONS_KEY] = globalThis[VIOLATIONS_KEY] || {});
  root.contractViolations = root.contractViolations || [];
  root.contractViolations.push(entry);
}

function matchPath(pattern, path) {
  const patternParts = pattern.split('/').filter((s) => s.length > 0);
  const pathParts = path.split('/').filter((s) => s.length > 0);
  if (patternParts.length !== pathParts.length) return null;
  const params = {};
  for (let i = 0; i < patternParts.length; i += 1) {
    const p = patternParts[i];
    if (p.startsWith(':')) {
      params[p.slice(1)] = decodeURIComponent(pathParts[i]);
    } else if (p !== pathParts[i]) {
      return null;
    }
  }
  return params;
}

/** `#/stages/mlst` -> `/stages/mlst`; no hash -> `/`. */
export function parseHash(hash) {
  const raw = String(hash || '').replace(/^#/, '');
  const cleaned = raw.replace(/^\/+/, '');
  return `/${cleaned}`.replace(/\/+$/, '') || '/';
}

/**
 * @param {object} options
 * @param {Array}  options.routes  `[{id, path, label, group, load, exportName}]`
 * @param {HTMLElement} options.outlet  the element pages are mounted into
 * @param {object} options.ctx      the DESIGN §7 context, frozen per mount
 */
export function createRouter({ routes, outlet, ctx }) {
  const listeners = new Map(); // event -> Set<fn>
  const moduleCache = new Map(); // id -> module
  const failed = new Map(); // id -> Error, so a known-missing module is not retried silently
  let currentTeardown = null;
  let currentRoute = null;
  let started = false;
  let mountToken = 0;
  let lastResolvedPath = null;

  function emit(event, payload) {
    const set = listeners.get(event);
    if (!set) return;
    for (const fn of Array.from(set)) {
      try {
        fn(payload);
      } catch (error) {
        // A subscriber that throws must not stop the others, and must not
        // leave the router half-mounted.
        console.error(`[router] a "${event}" subscriber threw`, error);
      }
    }
  }

  function current() {
    return currentRoute;
  }

  function on(event, fn) {
    if (!listeners.has(event)) listeners.set(event, new Set());
    listeners.get(event).add(fn);
    return () => off(event, fn);
  }

  function off(event, fn) {
    const set = listeners.get(event);
    if (set) set.delete(fn);
  }

  /**
   * `go(path)` sets the hash and routes, and returns the promise for the
   * route. The hash change fires `hashchange` in every browser, but relying on
   * that alone would make `go()` fire-and-forget, so `resolve()` runs here too
   * and the duplicate `hashchange` is recognised and ignored.
   */
  function go(path, { replace = false } = {}) {
    const target = `#${String(path).startsWith('/') ? path : `/${path}`}`;
    if (location.hash === target) return resolve();
    if (replace) {
      const url = `${location.pathname}${location.search}${target}`;
      history.replaceState(null, '', url);
      return resolve();
    }
    location.hash = target;
    return resolve();
  }

  function findRoute(path) {
    for (const route of routes) {
      const params = matchPath(route.path, path);
      if (params) return { route, params };
    }
    return null;
  }

  async function loadModule(route) {
    if (moduleCache.has(route.id)) return moduleCache.get(route.id);
    if (failed.has(route.id)) throw failed.get(route.id);
    try {
      const mod = await route.load();
      if (!mod || typeof mod[route.exportName] !== 'function') {
        throw new Error(
          `${route.specifier} does not export a function named ` +
            `"${route.exportName}". The §7 contract is ` +
            `export function mount(container, ctx) { ...; return teardown }.`
        );
      }
      moduleCache.set(route.id, mod);
      emit('modulePresent', { route: route.id, specifier: route.specifier });
      return mod;
    } catch (error) {
      const wrapped =
        error instanceof ModuleMissingError
          ? error
          : new ModuleMissingError(route.specifier, error);
      failed.set(route.id, wrapped);
      emit('moduleMissing', { route: route.id, specifier: route.specifier, error: wrapped });
      throw wrapped;
    }
  }

  /** Drop the module cache so a page another owner just wrote can be picked up. */
  function invalidate(id) {
    if (id === undefined) {
      moduleCache.clear();
      failed.clear();
    } else {
      moduleCache.delete(id);
      failed.delete(id);
    }
    // A retry has to be allowed to route the hash it is already on.
    lastResolvedPath = null;
  }

  function releaseCurrent() {
    if (typeof currentTeardown !== 'function') {
      currentTeardown = null;
      return;
    }
    const fn = currentTeardown;
    currentTeardown = null;
    try {
      fn();
    } catch (error) {
      // A teardown that throws has still released what it managed to; say so
      // rather than pretending the page was clean.
      console.error(`[router] teardown of "${currentRoute?.id}" threw`, error);
      recordViolation({ kind: 'teardown-threw', route: currentRoute?.id ?? null });
    }
  }

  function clearOutlet() {
    while (outlet.firstChild) outlet.removeChild(outlet.firstChild);
  }

  function panel(node) {
    clearOutlet();
    outlet.appendChild(node);
  }

  function text(tag, className, content) {
    const el = document.createElement(tag);
    if (className) el.className = className;
    if (content !== undefined && content !== null) el.textContent = String(content);
    return el;
  }

  function renderMissing(route, error) {
    const box = text('div', 'pending-module');
    box.appendChild(text('h3', null, 'Module not yet present'));
    box.appendChild(
      text(
        'p',
        null,
        `This page is declared by the router and its endpoint exists, but the ` +
          `script that draws it has not been written yet. That is an honest ` +
          `absence, not an error and not an empty result.`
      )
    );
    const dl = document.createElement('dl');
    const add = (k, v) => {
      const dt = text('dt', null, k);
      const dd = text('dd', 'mono', v);
      dl.append(dt, dd);
    };
    add('page', route.path);
    add('module', route.specifier);
    add('expected export', `${route.exportName}(container, ctx)`);
    add(
      'reason',
      error && error.cause && error.cause.message
        ? `${error.cause.name || 'Error'}: ${error.cause.message}`
        : error
        ? error.message
        : 'unknown'
    );
    box.appendChild(dl);
    const actions = text('p', null, '');
    const retry = text('button', null, 'Retry loading this module');
    retry.addEventListener('click', () => {
      invalidate(route.id);
      void resolve();
    });
    actions.appendChild(retry);
    box.appendChild(actions);
    panel(box);
  }

  function renderContractViolation(route, error) {
    console.error(`[router] TEARDOWN CONTRACT VIOLATION — ${error.message}`);
    recordViolation({
      kind: 'no-teardown',
      route: route.id,
      specifier: route.specifier,
      returned: typeof error.returned,
      at: new Date().toISOString(),
    });
    const box = text('div', 'message');
    box.dataset.severity = 'error';
    box.appendChild(text('h3', null, 'Page module broke the teardown contract'));
    box.appendChild(text('p', null, error.message));
    box.appendChild(
      text(
        'p',
        'panel-note',
        'The mount was refused rather than trusted: whatever the module ' +
          'registered before returning cannot be released from here. The ' +
          'violation is on window.__paDashboard.contractViolations and in the ' +
          'browser console.'
      )
    );
    panel(box);
  }

  function renderNoRoute(path) {
    const box = text('div', 'message');
    box.dataset.severity = 'warn';
    box.appendChild(text('h3', null, 'No such page'));
    box.appendChild(text('p', null, `The address "${path}" matches no declared page.`));
    const list = text('ul', null, '');
    for (const route of routes) {
      const item = text('li');
      const link = text('a', null, route.path);
      link.href = `#${route.path}`;
      item.appendChild(link);
      list.appendChild(item);
    }
    box.appendChild(list);
    panel(box);
  }

  async function resolve() {
    const token = (mountToken += 1);
    const path = parseHash(location.hash);
    // `go()` already routed this hash; the `hashchange` it provokes is the
    // same navigation, not a second one.
    if (path === lastResolvedPath) return;
    lastResolvedPath = path;
    const found = findRoute(path);

    emit('beforeNavigate', { path });
    releaseCurrent();

    if (!found) {
      currentRoute = null;
      renderNoRoute(path);
      emit('navigated', { path, route: null, params: null });
      return;
    }

    const { route, params } = found;
    currentRoute = { id: route.id, path, params: params || {} };
    emit('navigated', { path, route: route.id, params: params || {} });

    let mod;
    try {
      mod = await loadModule(route);
    } catch (error) {
      if (token !== mountToken) return; // a newer navigation won
      renderMissing(route, error);
      return;
    }
    if (token !== mountToken) return;

    const mount = mod[route.exportName];
    const container = document.createElement('div');
    container.className = 'page-body';
    container.dataset.page = route.id;
    clearOutlet();
    outlet.appendChild(container);

    let returned;
    try {
      returned = mount(container, ctx);
    } catch (error) {
      console.error(`[router] "${route.id}" threw while mounting`, error);
      const box = text('div', 'message');
      box.dataset.severity = 'error';
      box.appendChild(text('h3', null, 'This page failed to render'));
      box.appendChild(text('p', null, String((error && error.message) || error)));
      panel(box);
      return;
    }

    if (typeof returned !== 'function') {
      renderContractViolation(route, new TeardownContractError(route.id, returned));
      return;
    }

    currentTeardown = () => {
      clearOutlet();
      returned();
    };
    emit('mounted', { path, route: route.id, params: params || {} });
  }

  function onHashChange() {
    void resolve();
  }

  function start() {
    if (started) return;
    started = true;
    window.addEventListener('hashchange', onHashChange);
    void resolve();
  }

  function destroy() {
    if (!started) return;
    started = false;
    window.removeEventListener('hashchange', onHashChange);
    releaseCurrent();
    listeners.clear();
    moduleCache.clear();
    failed.clear();
    currentRoute = null;
  }

  return { start, destroy, go, current, on, off, invalidate, parseHash, routes };
}
