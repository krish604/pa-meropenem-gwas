/* The router/page integration contract (DESIGN §7).
 *
 * Every route the shell declares in `app.js#ROUTES` is loaded from disk and
 * mounted against a mock context. The DESIGN §7 contract is:
 *
 *     export function mount(container, ctx) { ...; return teardown }
 *
 * so each mount must return a function, and calling it must not throw. This is
 * the test that would have caught a route whose module landed but whose export
 * name, path parameter or teardown was wrong — a class of defect the pure
 * router unit tests cannot see because they use their own fake modules.
 *
 * The mock `api.get` never resolves: this asserts the *synchronous* mount and
 * teardown contract, and keeps the async data-loading path (which needs a real
 * server) out of a unit test.
 */

import './dom-shim.mjs';

import test from 'node:test';
import assert from 'node:assert/strict';

import { ROUTES } from '../../static/app.js';

/* --------------------------------------------------------------- the ctx */

function mockBadge() {
  return { key: 'completed', label: 'completed', glyph: '✓', tone: 'ok', reason: '' };
}

function makeCtx() {
  const never = () => new Promise(() => {});
  return {
    api: { get: never, post: never },
    router: {
      go() {},
      current: () => ({ path: '/', params: {} }),
      on: () => () => {},
      off() {},
    },
    theme: { tokens: {} },
    banner: {
      mode: null,
      power: () => ({ n: null, flag: null, underpowered: null }),
      underpowered: () => null,
      text: () => 'not reported',
      mount() {},
    },
    badges: {
      state: () => mockBadge(),
      tone: () => 'ok',
      render: () => document.createElement('span'),
      adopt() {},
      notProduced: 'not produced',
      notAssessed: 'not assessed',
      notReported: 'not reported',
      timings: 'timings not recorded',
    },
    on: { runState: () => () => {}, delta: () => () => {} },
  };
}

/* ------------------------------------------------------------- the test */

for (const route of ROUTES) {
  test(`route ${route.id} (${route.path}) mounts and returns a teardown`, async () => {
    const mod = await route.load();
    assert.equal(
      typeof mod[route.exportName],
      'function',
      `${route.specifier} must export a function named ${route.exportName}`
    );

    const container = document.createElement('div');
    const ctx = makeCtx();
    const teardown = mod[route.exportName](container, ctx);

    assert.equal(
      typeof teardown,
      'function',
      `${route.specifier}#${route.exportName} returned ${typeof teardown} instead ` +
        `of a teardown function (DESIGN §7)`
    );

    assert.doesNotThrow(() => teardown(), `teardown for ${route.id} threw`);
  });
}

/* ------------------------------------- the isolate-detail path parameter */

test('the isolate-detail route resolves a sample_id parameter', async () => {
  const route = ROUTES.find((r) => r.id === 'isolate');
  assert.ok(route, 'the isolate detail route must be declared');
  assert.equal(route.path, '/isolates/:sample_id');

  // `parseHash` + the route matcher must carry the sample id through. This is
  // the wiring DATA reported as uncertain: a param route is only useful if the
  // hash path actually matches it.
  const { createRouter } = await import('../../static/router.js');
  const outlet = document.createElement('div');
  const seen = [];
  const router = createRouter({
    routes: [route],
    outlet,
    ctx: makeCtx(),
  });
  router.on('mounted', ({ params }) => seen.push(params));
  await router.go('/isolates/TEST_PA_001');
  assert.deepEqual(seen, [{ sample_id: 'TEST_PA_001' }]);
  router.destroy();
});
