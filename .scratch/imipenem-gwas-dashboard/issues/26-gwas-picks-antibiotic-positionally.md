# 26: The GWAS picks its antibiotic positionally, with no silent default allowed

**What to build:** An explicit, validated antibiotic selection for stage 12, so
the analysed drug is a stated fact rather than a list ordering.

**Blocked by:** nothing. Stage 11 already writes one phenotype file per
configured antibiotic.

**Status:** ready-for-agent

**Do not fix opportunistically.** This is recorded so it is not lost, and so it
is fixed deliberately. It is not a blocker for ticket 14.

## The defect

`papipeline/stages/gwas.py:134`:

```python
antibiotic = config.antibiotics[0]
```

The GWAS analyses the **first** antibiotic in `config/science.yaml` →
`antibiotics:`. Nothing else in the repo behaves this way:

| component | how it picks the drug | on an unknown drug |
|---|---|---|
| `stages/phenotype.py:97-102` | reads `antibiotic` from each row | **hard error**, naming the configured set |
| `loader.py:451-460` | cross-validates `science.yaml` against `config/antibiotics.tsv` | **hard error** on a mismatch |
| `stages/gwas.py:134` | `config.antibiotics[0]` | **silently runs on a different drug** |

Reorder `antibiotics:` in `science.yaml` — a natural edit, and one that changes
nothing else in the pipeline — and the GWAS output becomes meropenem results
under an imipenem heading. No warning, no error, no test failure. The same
happens if a second drug is ever added "just to record its prevalence": the
stages that are meant to report both will silently report one.

## Why it matters: the no-silent-drug-default requirement

The pipeline has a standing rule that the analysed drug is configuration, not
assumption, and that an unresolvable drug is a **refusal** rather than a
default. `phenotype.py` and the loader both implement that. This one line is
the exception, and it is the exception that matters most, because it is the stage
that produces the headline result.

A positional pick is also unreviewable: nothing in the code, the docstring, or
the emitted output records *which* drug was analysed. The `Raises` section of
`gwas.run()` does not mention it, because there is nothing to raise.

## What to build

- An explicit target, from config rather than from list order — e.g.
  `gwas.target_antibiotic` in `science.yaml`, validated to be a member of
  `antibiotics:` and of `config/antibiotics.tsv`.
- A hard error when the target is unset, rather than falling back to `[0]`. The
  absence of a fallback is the point; a default here reintroduces exactly the
  silent behaviour being removed.
- Record the target on the run result and in the report, so a GWAS output states
  its drug in the artefact and not only in the config.
- A test that the GWAS runs on the configured target when several antibiotics are
  configured — the case that currently passes by accident.

## Red-first, and the trap

Write the test with **two or more** antibiotics configured and the target set to
the *second*. A test with one antibiotic cannot distinguish an explicit
selection from a positional default, because they agree — which is precisely why
this defect has survived. This is the same small-N lesson as the cohort and cap
tests: a fixture too small to separate the correct implementation from the
plausible wrong one proves nothing.

## Relationship to the audit that found this

From the antibiotic-selection audit, alongside finding that AST parsing is
genuinely generic and that drug selection is otherwise config-driven and
strictly validated. This was the single soft spot in that audit.
