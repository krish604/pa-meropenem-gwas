# `convergence_real/` - a synthetic cohort with a KNOWN convergence answer

Nine samples in three lineages of three, built so that every category
`stages.convergence.classify` can return is reachable and each determinant has
exactly one correct answer. Nothing here is derived from `PDC_essential.tsv`,
which is real clinical data and is never committed.

`papipeline/testing/synthetic.py` is frozen, so these are new files rather than
a change to it. They are committed, byte-stable, and **copied into `tmp_path`
before use** - a test must never write into this directory or into `results/`.

## The cohort

| lineage | samples |
| --- | --- |
| `L1` | `A1`, `A2`, `A3` |
| `L2` | `B1`, `B2`, `B3` |
| `L3` | `C1`, `C2`, `C3` |

Nine analysed samples. `convergence.widespread_fraction` is 0.75, so
`widespread_background` needs **7 or more** carriers out of these nine - a
number that is unreachable if the denominator is the 967-isolate manifest and
is reached comfortably if the denominator is the nine actually analysed. That is
the whole point of this fixture: it is sensitive to the denominator, which is
why `blaOXA-1` is carried by eight.

## The planted answers

| determinant | source | carriers | n | lineages | expected category |
| --- | --- | --- | --- | --- | --- |
| `blaOXA-1` | stage 4, `04_amr.tsv` | all but `C3` | 8 | 3 | `widespread_background` (8/9 = 0.889 > 0.75) |
| `fosA` | stage 4, `04_amr.tsv` | `B1 B2 B3` | 3 | 1 | `lineage_associated` (3 carriers clears `rare_max_samples` = 2, then fails on one lineage) |
| `oprD:A100T` | stage 6, `regulator_variants.tsv` | `A1 B1 C1` | 3 | 3 | `recurrent_convergent` (3 >= `min_independent_lineages` = 2) |
| `mexR:I24A` | stage 6, `regulator_variants.tsv` | `A1 A2 A3` | 3 | 1 | `lineage_associated` (same three carriers as `oprD:A100T`, opposite distribution) |

The last two are the pair the convergence question actually turns on: identical
carrier counts, identical rarity, opposite lineage distributions, opposite
answers. A stage that got the lineage logic wrong would collapse them into one
category, so these two rows cannot both be right by accident.

`fosA` and `mexR:I24A` carry three carriers rather than one on purpose. At
`rare_max_samples = 2`, a single-carrier determinant stops at `rare_isolated`
and never reaches the lineage test at all - so a one-carrier fixture would pass
whether or not the lineage branch worked.

## Schemas

The headers are copied verbatim from
`tests/unit/test_convergence_real_inputs.py`, which is the authoritative
description of a real-shaped intermediate. A hand-written header one column
short fails in `read_tsv` with a `DataContractError` about the *fixture*, which
says nothing about the behaviour under test.

* `04_amr.tsv` - stage 4's `stages/04_amr.tsv`. Only `determinant_type == AMR`
  rows are read; `load_amr_names` drops everything else.
* `regulator_variants.tsv` - stage 6's `regulators/regulator_variants.tsv`.
  `variant` must not be `.`: `read_tsv` maps `.` to `None` and
  `load_variant_names` skips such a row, which yields a table that looks
  populated and produces an empty mutation namespace.
* `tree_metadata.tsv` - stage 9's crosswalk, `phylogeny_dir/tree_metadata.tsv`.
  `lineage_label` is the MLST sequence type under `lineage.method: st`.