// State vocabulary and colour, shared by the canvas, the counters and the
// detail panel so that one state is never two different colours.
//
// These seven states are exactly StageState in papipeline/execution/state.py.
// Nothing else is a state: the UI has no states of its own, and in
// particular has no "busy" or "almost done" that the engine cannot produce.

export const STATES = {
  PENDING:    { label: 'PENDING',    color: '#64748b', glow: 0.0,  order: 0 },
  RUNNING:    { label: 'RUNNING',    color: '#3b82f6', glow: 1.0,  order: 1 },
  SUCCEEDED:  { label: 'SUCCEEDED',  color: '#22c55e', glow: 0.45, order: 2 },
  RETRYING:   { label: 'RETRYING',   color: '#f59e0b', glow: 0.8,  order: 3 },
  FAILED:     { label: 'FAILED',     color: '#ef4444', glow: 0.7,  order: 4 },
  INVALID:    { label: 'INVALID',    color: '#f97316', glow: 0.7,  order: 5 },
  INCOMPLETE: { label: 'INCOMPLETE', color: '#a855f7', glow: 0.7,  order: 6 },
};

// Per-stage accents, so a stage keeps its identity across state changes —
// the colour says *which stage*, the ring and pulse say *what state*. The
// stage palette is a fixed hue per stage; only SUCCEEDED/RUNNING/FAILED
// override it, because those are the three an operator must spot instantly.
export const STAGE_HUES = {
  validation:           '#22c55e',
  annotation:           '#22d3ee',
  mlst:                 '#eab308',
  amr:                  '#f97316',
  regulators:           '#d946ef',
  structural_variants:  '#2dd4bf',
  mechanisms:           '#6366f1',
  virulence:            '#ec4899',
  pangenome:            '#a855f7',
  phylogeny:            '#14b8a6',
  phenotype:            '#38bdf8',
  gwas:                 '#8b5cf6',
  convergence:          '#10b981',
  cooccurrence:         '#f59e0b',
  integration:          '#e2e8f0',
  reporting:            '#94a3b8',
};

export function stageHue(stage) {
  return STAGE_HUES[stage] || '#64748b';
}

export function stateStyle(state) {
  return STATES[state] || STATES.PENDING;
}

export function stateColor(state) {
  return stateStyle(state).color;
}

/** True when the stage is doing work right now. Drives the idle overlay. */
export function isLive(state) {
  return state === 'RUNNING' || state === 'RETRYING';
}

/**
 * The colour a node should actually be drawn in.
 *
 * Identity colours by stage; alarm colours win, because a failed stage
 * that still looks like its own hue is easy to miss in a dense network.
 */
export function nodeColor(stage, state) {
  if (state === 'FAILED' || state === 'INVALID' || state === 'INCOMPLETE') {
    return stateColor(state);
  }
  return stageHue(stage);
}

export const EVENT_STYLE = {
  TASK_CREATED:    { color: '#94a3b8' },
  TASK_STARTED:    { color: '#22d3ee' },
  TASK_COMPLETED:  { color: '#22c55e' },
  TASK_VALIDATED:  { color: '#4ade80' },
  TASK_FAILED:     { color: '#ef4444' },
  TASK_INVALID:    { color: '#f97316' },
  TASK_INCOMPLETE: { color: '#a855f7' },
  TASK_RETRY:      { color: '#f59e0b' },
  TASK_RESUMED:    { color: '#38bdf8' },
};

export function eventColor(name) {
  return (EVENT_STYLE[name] || { color: '#94a3b8' }).color;
}

/** The seven states, for the legend. */
export const STATE_LIST = Object.keys(STATES);
