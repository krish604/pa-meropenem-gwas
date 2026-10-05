import './dom-shim.mjs';

import test from 'node:test';
import assert from 'node:assert/strict';

import {
  createRouter,
  parseHash,
  ModuleMissingError,
  TeardownContractError,
} from '../../static/router.js';

test('parseHash normalises a hash to a leading-slash path', () => {
  assert.equal(parseHash('#/stages/mlst'), '/stages/mlst');
  assert.equal(parseHash('/stages/mlst'), '/stages/mlst');
  assert.equal(parseHash('#'), '/');
  assert.equal(parseHash(''), '/');
  assert.equal(parseHash('#/'), '/');
});

function makeOutlet() {
  return document.createElement('div');
}

test('the router calls the teardown of the page it leaves', async () => {
  let torn = 0;
  const routes = [
    { id: 'a', path: '/a', label: 'A', group: 'g', exportName: 'mount', load: async () => ({ mount: () => () => { torn += 1; } }) },
    { id: 'b', path: '/b', label: 'B', group: 'g', exportName: 'mount', load: async () => ({ mount: () => () => { torn += 1; } }) },
  ];
  const router = createRouter({ routes, outlet: makeOutlet(), ctx: {} });
  router.start();
  await router.go('/a');
  assert.equal(torn, 0);
  await router.go('/b');
  assert.equal(torn, 1, 'navigating away must release the previous page');
  router.destroy();
  assert.equal(torn, 2, 'destroy must release the current page');
});

test('a page that returns no teardown is refused and recorded', async () => {
  delete globalThis.__paDashboard;
  const routes = [
    { id: 'leaky', path: '/leaky', label: 'L', group: 'g', exportName: 'mount', load: async () => ({ mount: () => undefined }) },
  ];
  const outlet = makeOutlet();
  const router = createRouter({ routes, outlet, ctx: {} });
  await router.go('/leaky');
  const violations = globalThis.__paDashboard?.contractViolations || [];
  assert.equal(violations.length, 1);
  assert.equal(violations[0].kind, 'no-teardown');
  assert.equal(violations[0].route, 'leaky');
  // The page is replaced by a panel that says the contract was broken, never
  // left as a blank page.
  const text = outlet.textContent;
  assert.match(text, /teardown contract/);
  assert.notEqual(text.trim(), '');
  router.destroy();
});

test('a module that is not on disk degrades to a named panel', async () => {
  const routes = [
    { id: 'missing', path: '/missing', label: 'M', group: 'g', exportName: 'mount', specifier: './pages/missing.js', load: async () => { throw new Error('not found'); } },
  ];
  const outlet = makeOutlet();
  const router = createRouter({ routes, outlet, ctx: {} });
  await router.go('/missing');
  const text = outlet.textContent;
  assert.match(text, /Module not yet present/);
  assert.match(text, /\.\/pages\/missing\.js/);
  router.destroy();
});

test('an unknown path renders "No such page", never a blank', async () => {
  const routes = [{ id: 'a', path: '/a', label: 'A', group: 'g', exportName: 'mount', load: async () => ({ mount: () => () => {} }) }];
  const outlet = makeOutlet();
  const router = createRouter({ routes, outlet, ctx: {} });
  await router.go('/nope');
  assert.match(outlet.textContent, /No such page/);
  router.destroy();
});

test('a teardown that throws is recorded, not swallowed', async () => {
  delete globalThis.__paDashboard;
  const routes = [
    { id: 'thrower', path: '/thrower', label: 'T', group: 'g', exportName: 'mount', load: async () => ({ mount: () => () => { throw new Error('boom'); } }) },
    { id: 'next', path: '/next', label: 'N', group: 'g', exportName: 'mount', load: async () => ({ mount: () => () => {} }) },
  ];
  const router = createRouter({ routes, outlet: makeOutlet(), ctx: {} });
  await router.go('/thrower');
  await router.go('/next');
  const violations = globalThis.__paDashboard?.contractViolations || [];
  assert.equal(violations.length, 1);
  assert.equal(violations[0].kind, 'teardown-threw');
  router.destroy();
});

test('the error classes name what went wrong', () => {
  assert.equal(new ModuleMissingError('./x.js', new Error('e')).name, 'ModuleMissingError');
  const error = new TeardownContractError('a', undefined);
  assert.equal(error.name, 'TeardownContractError');
  assert.match(error.message, /undefined/);
  assert.match(error.message, /returned/);
});
