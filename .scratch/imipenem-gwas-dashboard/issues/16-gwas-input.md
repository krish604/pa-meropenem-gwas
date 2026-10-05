# 16: GWAS input: three feature families and unique-pattern reduction

**What to build:** The association model's actual input: all three required feature families together, with the multiple-testing threshold derived from unique variant patterns in a separate, inspectable step. This is the ticket the reviewer scrutinises hardest, because the whole analysis rests on it.

**Blocked by:** 14, 15

**Status:** ready-for-agent

The eight boxes below are ticked per-item, and the wording of a tick is part of
what it claims. Two of them are ticked as **code path present, not exercised**,
because the ticket's "verified individually" is not yet true of any run. See
*Still open* for what would change each.

- [x] The input contains SNP/indel variants from the PAO1 alignment, gene presence/absence from the pangenome, and the OprD features. All three, verified individually. **Partly: the code path handles all three and the feature-matrix types are checked per column, but nothing populates two of them yet** — see *Still open* for the OprD producer (ticket 15) and for the SNP/indel family reaching pyseer through a gene-P/A `.Rtab` rather than `--vcf`.
- [x] Gene presence/absence alone is not treated as a sufficient input.
- [x] Unique variant patterns are reduced in a separate step whose output is written to disk and can be read, and the number of survivors is reported.
- [x] pyseer's helper script name and flags are confirmed against the installed package, not from memory. **pyseer's own** flags (`--output-patterns`, `--pres`, `--lrt-pvalue`, `--filter-pvalue`, `--min-af`, `--max-af`, `--cpu`) were read from `pyseer --help` *and* from `pyseer/__main__.py`'s argparse block; the vendored helper's flags were read from its own live `--help`, since it is not installed at all (see below).
- [x] A mixed linear model is used, with the phylogeny from stage 9. **Code path present: `--lmm` is passed and kinship is built from stage 9's core-SNP alignment when it exists. Not exercised, because no pipeline code constructs `PyseerEngine` yet** — see *Still open*.
- [x] Minimum samples per group is enforced so a hit supported by two isolates cannot be reported.
- [x] A feature confounded with a single lineage is reported as lineage-linked, not as an association.
- [x] Unitigs are not claimed among the supported families.

## Progress

### Step 12a, unique-pattern reduction — DONE

`papipeline/stages/gwas.py`: `reduce_unique_patterns()`,
`count_unique_patterns()`, `tally_patterns()`,
`parse_count_patterns_output()`, `UniquePatternReduction`, and
`PyseerEngine.run()` as two passes. 75 tests in
`tests/unit/test_gwas_unique_patterns.py`, against synthetic pyseer output in
pyseer 1.1.2's own formats.

**The decision the Phase 1 finding required has been taken: the helper is
vendored, not substituted.** `scripts/gwas/count_patterns.py`, body
byte-identical to upstream, header recording source URL, retrieval date, the
SHA-256 of the upstream file and the SHA-256 of the vendored logic. The test
suite recomputes the logic digest and asserts it equals the one recorded in the
header, so an edit to the correction fails CI and an edit to the recorded digest
without the correction fails CI too. `config.gwas.unique_patterns.method`
accepts only `bonferroni`, so a second name cannot appear in configuration
without the pipeline refusing it.

Five things that were not in the ticket and had to be settled:

1. **The count is cross-checked three ways.** Against an independent tally of
   the same file, against pyseer's own printed threshold, and against the
   `%d tested variants` count pyseer writes to stderr. Any disagreement stops
   the run rather than picking one. Note the three-way check validates the
   *input* population and the helper's arithmetic — not the divisor, which is
   the distinct-pattern count from item 2. See item 5 for what this check
   cannot do.
2. **The divisor is the *unique-pattern* count; pyseer's *tested* count is the
   reduction's input.** Two numbers, easy to confuse, and confusing them
   changes the threshold — so both appear as separate rows in the audit record
   and neither is described as the other. The tested count is pyseer's own, not
   the input matrix's variant count: variants outside the frequency bounds and
   variants stopped by the prefilter gate never reach the patterns file. pyseer
   writes a pattern for every non-prefilter variant *before* evaluating its own
   p-value gate (`pyseer/__main__.py`), so `--filter-pvalue` and
   `--lrt-pvalue` do not narrow the file; they are still stated explicitly so a
   future change to pyseer's defaults cannot. What bounds the pattern set is
   the frequency bounds in the shared command; what establishes the input
   population is the stderr cross-check.

   **This section has been wrong twice, in opposite directions, and both
   corrections are recorded here rather than quietly overwritten.** It first
   claimed those flags protected the denominator — they do not, and they did
   not. It then claimed the tested count *was* the denominator — it is not; it
   is the reduction's input, and the divisor is the distinct-pattern count.
   Review history lives here, on the ticket, rather than in
   `docs/scientific_rules.md`, which states only what is true now.

   It also costs one extra full model fit, because pyseer 1.1.2 cannot report
   the number of distinct patterns without fitting.

3. **A header-only association table is now an empty result, not an error —
   but a single non-header line is neither.** Before 12a the header-only case
   was unreachable, since the stage ran at pyseer's default `--lrt-pvalue 1`.
   Now that pass 12b runs at a derived cutoff, "nothing passes" is the ordinary
   null outcome. A truncated write or a stray diagnostic presents the same way,
   though, and reading it as "no associations" would convert a broken run into a
   reported negative result — so that case is refused with the columns found.
   `read_tsv` is unchanged and still refuses header-only files, which is right
   for a table that must hold records.
4. **The helper's environmental arguments are validated, not inherited.** Its
   `--memory`/`--cores`/`--temp` defaults (1024, 1, `/tmp`) are machine facts
   and come from the overlay, with accessors that raise rather than fall back.
   Both the patterns path *and* the scratch directory are validated, because
   the helper interpolates both into one unquoted `shell=True` string; the
   memory floor is checked because the helper subtracts 10 before handing the
   value to `sort`; and `--cores > 1` is probed against `sort --help` because
   `--parallel` is not in POSIX.
5. **The stderr cross-check can only refute, never establish.** If pyseer is
   not writing that line — because a caller captured its stderr elsewhere, or a
   future pyseer stops writing it — the check logs and stands down, and the
   reduction proceeds on the other two. `PYSEER_TESTED_COUNT_PATTERN` is pinned
   against the installed package precisely because a check that fails open is
   only honest if you know it is one.

### Still open

- [ ] **No end-to-end run against the real pyseer binary exists, and adding one
      inside the test suite would cross a mode boundary.** Every
      `PyseerEngine.run` test replaces `_invoke`, so nothing has yet *observed*
      pyseer writing its patterns file, printing its header, or emitting
      `%d tested variants`. The format assumptions are instead checked against
      the installed package by **importing** it:
      `test_the_pattern_encoding_is_pyseers_not_a_packed_byte_vector` calls
      `pyseer.input.hash_pattern` and asserts the naive encoding differs (it
      does — pyseer hashes the `int64` view, 8 bytes per sample, and a first
      draft of the fixtures hashed one byte per value and was wrong);
      `test_the_synthetic_association_table_has_pyseers_own_column_count` and
      `test_the_synthetic_table_rows_come_from_pyseers_formatter` build rows
      through `pyseer.utils.format_output`. That closes the gap where the test
      file agreed only with itself, but it does not observe the tool.
      spec.md D8 defines TEST as running no real tools, so an end-to-end run
      belongs in a mode that permits tools, or behind an explicit opt-in rather
      than in `pytest`. **This is the largest remaining risk in the 12a work and
      it deserves its own ticket.**
- [ ] **Blocked on 15.** No `oprD_absent` / `oprD_LoF` producer exists, so the
      third family is not yet populated from anything real. `build_input()`
      will carry such a column the moment one appears; the ticket's
      "all three, verified individually" is only true of the code path, not of a
      run.
- [ ] The `--pres` Rtab question below: `PyseerEngine._write_rtab` already
      synthesises an `.Rtab` from the feature matrix, so pyseer is fed a
      conformant file — but that matrix is gene presence/absence only. The
      SNP/indel family reaches pyseer through that same Rtab, which is a
      narrowing of what pyseer can take (`--vcf` is the route for cohort
      variants) and has not been resolved. It is upstream of the reduction and
      out of this ticket's scope.
- [ ] **Nothing constructs `PyseerEngine` from pipeline code.**
      `gwas.py:run()` refuses any mode but `TEST` *before* it looks at an
      injected engine, so the engine-injection path is unreachable — which
      includes `scripts/run_downstream_stages.py`, whose call has always
      raised. The guard is pre-existing and correct as a refusal; the dead call
      site is not, and should be made explicit or removed.
- [ ] `gwas.py:run()` still refuses REAL. Nothing here changed that, and it
      should not: the reduction is TEST-only until a real engine is wired.

## Finding from Phase 1: the unique-pattern helper must be vendored

pyseer 1.1.2 is the highest conda-installable version (1.2.0+ require `glmnet_py`,
which is PyPI-only). It ships **no** `pyseer_unique_patterns`,
`pyseer_gwas` or `pyseer_combine_features` script, and has no `unique_patterns`
module. Its wheel's `entry_points.txt` declares exactly five console scripts —
`pyseer`, `annotate_hits_pyseer`, `scree_plot_pyseer`, `square_mash` and
`phandango_mapper` — read off the installed `dist-info`, not from memory.

Per pyseer's own documentation the unique-pattern threshold comes from
`scripts/count_patterns.py` **in the GitHub source tree**, which is never
packaged. Its entire logic is `LC_ALL=C sort -u <patterns> | wc -l` then
`alpha / count`, fed by pyseer's own `--output-patterns` flag (which 1.1.2 does
have, confirmed from `pyseer --help`).

**Resolved:** vendored verbatim into `scripts/gwas/count_patterns.py` with
attribution, and tested against a known pattern file — both against
`sort -u | wc -l` computed independently, and by running the script and reading
its own `--help` for its flag names. The correction method was not substituted.
See `docs/environment-arm64.md` §4 and `docs/scientific_rules.md` §3a.

Also: `--pres` expects a roary/piggy `.Rtab` matrix, not a TSV, so the pangenome
stage must emit an Rtab. Still open — see above.
