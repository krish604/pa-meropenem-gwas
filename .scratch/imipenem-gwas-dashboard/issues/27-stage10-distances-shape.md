# 27: Stage 10's `--distances` file: decide packed vs square (four sources disagree)

**What to build:** A decision on which on-disk shape `similarity.tsv` is
contractually, plus the code and doc changes that decision implies. Today four
authoritative-looking sources disagree, and the disagreement is latent only
because the one consumer that needs the square form is unwired.

**Blocked by:** None — this is a design decision the build deliberately
surfaced rather than resolved.

**Blocks:** Correct wiring of the pyseer adapter (`papipeline/gwas_real/
adapter.py`, built, NOT wired), which takes its `--distances` path from
`STAGE_TABLES["similarity"][0]`.

**Status:** ready-for-human

## The four voices

1. **Packed (what actually lands on disk).** `STAGE_TABLES["similarity"]` is
   `("similarity.tsv", ("sample_id", "distances"))`
   (`papipeline/execution/contracts.py:118-124`, with a rationale comment: a
   square header would change with the cohort, so the vector travels keyed on
   the cohort instead). `tests/unit/test_similarity_contract.py:92` pins
   `required_columns(STAGE) == (ID_COLUMN, "distances")`, and contract
   validation (`has_columns`) requires a `distances` column — a square header
   would FAIL validation.
   **Mechanism:** `stages/similarity.py::run(out_path=...)` first calls
   `write_matrix` (square), then `run.py:1630-1636` calls
   `write_tsv(table_path(stage_dir, "similarity"), similarity_rows, ...)` on
   **the same path** — the packed two-column form overwrites the square
   matrix. Last writer wins. Verified on the TEST run:
   `results/test/intermediate/stages/similarity.tsv` is 2 columns,
   21 lines (`distances` holds `A=0.0;B=0.3;...`).

2. **Square (the stage's own writer + its docstring).**
   `stages/similarity.py`'s module docstring says "Output is a square matrix,
   `sample_id` leading"; `write_matrix` writes exactly that; and
   `tests/unit/test_similarity_stage.py:270-295` pins `write_matrix`'s
   square output (calling it directly, not through dispatch).

3. **`docs/data_contract.md:462` says square — with the wrong filename.** The
   row reads `` `10_similarity.tsv` ``; the real file is `similarity.tsv`
   (there is no numbered-name convention on disk).

4. **The consumer needs square.** `pyseer/input.py::load_structure` does
   `pd.read_table(infile, index_col=0)` then
   `m.loc[intersecting_samples, intersecting_samples]`: sample IDs as BOTH
   index and columns. Empirically, running that against the dispatched packed
   file raises **KeyError** (`m.loc[rows, rows]` finds no such columns).
   `pyseer/__main__.py:220` additionally makes `--lmm --lineage` exit without
   a `--distances` matrix, so lineage-as-a-model-term needs this file too.
   The adapter's `_read_matrix_ids` (`papipeline/gwas_real/adapter.py:341-361`)
   parses a square header and, on the packed file, returns `["distances"]` —
   treating the column name as a sample ID.

## Why nothing breaks today

- The adapter is **unwired**: nothing executed reads `similarity.tsv` as a
  pyseer input. The wired path (`_build_gwas_engine` → `PyseerEngine`) passes
  `--similarity` from its own kinship matrix and never `--lineage`/`--distances`.
- The adapter's own tests pass because their helper `build_distances`
  (`tests/unit/test_gwas_real_adapter.py:166-180`) calls
  `similarity.run(out_path=path)` **directly** — receiving the square form —
  and therefore never exercises the dispatch overwrite. The tests validate a
  shape the pipeline does not produce.

## The decision

- **Keep packed:** `docs/data_contract.md:462` gets filename + shape corrected;
  the adapter must convert packed → square (pyseer cannot read packed, full
  stop). `test_similarity_stage`'s square claim about `write_matrix` stays true
  of the function, not the artifact.
- **Switch to square:** `STAGE_TABLES`, `test_similarity_contract.py:92` (and
  possibly `test_similarity_stage`) change — a test edit, which the current
  build forbids (AGENTS.md rule 3, "no test is edited on this build");
  `run.py`'s `write_tsv` overwrite goes away; the units sidecar
  (`similarity.tsv.units.json`) already matches either reading.

Related latent note: `dashboard/server/badges.py` checks `refused`
(`stage in REAL_REFUSING_STAGES and run_mode == "REAL"`) BEFORE `completed`,
so a hypothetical pyseer-resolved REAL run in which gwas completed would badge
it `refused`. REAL never runs here (`allow_real_mode` false); record only.

**Found:** 2026-10-08 during the meropenem build's verification pass, while
checking whether `stages/gwas.py`'s refusal message could drop its stage-10
clause. The clause turned out to be true of the dispatched artifact; the
contradiction behind it is this ticket.
