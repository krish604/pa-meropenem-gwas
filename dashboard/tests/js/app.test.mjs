import './dom-shim.mjs';

import test from 'node:test';
import assert from 'node:assert/strict';

import { ABSENT, basisOf, createBadges, createBanner, value } from '../../static/app.js';
import { cell, isAbsenceWord } from '../../static/pages/tables.js';

/* ----------------------------------------------------------------- badges */

test('a missing badge state is never blank, 0 or a default badge', () => {
  const badges = createBadges();
  const missing = badges.state(null);
  assert.notEqual(missing.label.trim(), '');
  assert.equal(missing.label, ABSENT.NOT_REPORTED);
  assert.notEqual(missing.label, '0');
  assert.ok(missing.reason.trim().length > 0, 'a missing state must carry a reason');
});

test('an unrecognised state is reported, not coerced onto one of the six', () => {
  const badges = createBadges();
  const unknown = badges.state('definitely_not_a_badge');
  assert.match(unknown.label, /not recognised/);
  assert.ok(unknown.reason.includes('definitely_not_a_badge'));
  assert.notEqual(unknown.label, 'completed');
});

test('the six badges carry label, glyph and tone together', () => {
  const badges = createBadges();
  for (const key of ['completed', 'running', 'failed', 'refused', 'not_assessed', 'not_run']) {
    const badge = badges.state(key);
    assert.equal(badge.key, key);
    assert.ok(badge.label.trim());
    assert.ok(badge.glyph.trim());
    assert.ok(badge.tone.trim());
  }
  assert.equal(badges.state('not_assessed').label, 'not assessed');
  assert.equal(badges.state('not_run').label, 'not run');
});

test('a server badge payload wins over the fallback', () => {
  const badges = createBadges();
  badges.adopt({ completed: { label: 'ran', glyph: '✓', tone: 'ok' } });
  assert.equal(badges.state('completed').label, 'ran');
  assert.equal(badges.source(), 'server');
});

/* ----------------------------------------------------------------- banner */

test('the D3 banner takes N from basis.n, never a page count', () => {
  const banner = createBanner();
  const fromBasis = banner.power(null, { n: 7, min_samples: 20 });
  assert.equal(fromBasis.n, 7);
  assert.equal(fromBasis.underpowered, true);
  assert.match(fromBasis.flag, /n=7, underpowered, not a finding/);
});

test('the flag drops "underpowered" when the cohort is at or above the minimum', () => {
  const banner = createBanner();
  const adequate = banner.power(null, { n: 900, min_samples: 3 });
  assert.equal(adequate.underpowered, false);
  assert.equal(adequate.flag, 'n=900, not a finding');
  assert.doesNotMatch(adequate.flag, /underpowered/);
  // Below the minimum the pipeline's verbatim wording is kept.
  const low = banner.power(null, { n: 2, min_samples: 3 });
  assert.equal(low.underpowered, true);
  assert.match(low.flag, /n=2, underpowered, not a finding/);
  // With no minimum recorded, nothing is claimed either way.
  const unknown = banner.power(null, { n: 5 });
  assert.equal(unknown.underpowered, null);
  assert.match(unknown.flag, /underpowered/);
});

test('an explicit n wins over basis.n, and a missing n is not reported', () => {
  const banner = createBanner();
  assert.equal(banner.power(999, { n: 5 }).n, 999);
  const none = banner.power(null, null);
  assert.equal(none.n, null);
  assert.equal(none.flag, null);
  assert.equal(none.underpowered, null);
  assert.match(banner.text(null, null), /not reported/);
});

test('basisOf reads meta.basis, bases.<key> and power.n, and returns null otherwise', () => {
  assert.deepEqual(basisOf({ meta: { basis: { n: 3 } } }), { n: 3 });
  assert.deepEqual(basisOf({ bases: { stages: { n: 2 } } }, 'stages'), { n: 2 });
  assert.equal(basisOf({ power: { n: 4 } }).n, 4);
  assert.equal(basisOf({}), null);
  assert.equal(basisOf(null), null);
});

/* ---------------------------------------------------------------- absence */

test('value() never renders null/undefined/[] as a blank', () => {
  assert.match(value(null, 'no source').textContent, /no source/);
  assert.match(value(undefined).textContent, new RegExp(ABSENT.NOT_REPORTED));
  assert.match(value([], 'not produced', { emptyIsAbsent: true }).textContent, /not produced/);
  assert.equal(value(0).textContent, '0');
  assert.equal(value('R').textContent, 'R');
});

test('cell() distinguishes a measured zero from an absent value', () => {
  assert.equal(cell(0).textContent, '0');
  assert.match(cell([], 'present table, no row').textContent, /none recorded/);
  assert.match(cell(null, 'no file').textContent, /no file/);
  assert.match(cell('').textContent, /empty string/);
  assert.equal(cell('oprD').textContent, 'oprD');
});

test('isAbsenceWord recognises the server absence words only', () => {
  for (const word of ['not assessed', 'not produced', 'not reported', 'not run']) {
    assert.equal(isAbsenceWord(word), true);
  }
  assert.equal(isAbsenceWord('R'), false);
  assert.equal(isAbsenceWord('0'), false);
  assert.equal(isAbsenceWord(''), false);
  assert.equal(isAbsenceWord(0), false);
});
