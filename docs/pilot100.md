# Pilot-100 analysis

A pilot analysis of **100** *Pseudomonas aeruginosa* assemblies, run against
the existing tested pipeline in `papipeline/`, `config/`, `workflow/`,
`scripts/` and `tests/`.

> **STATUS: cohort prepared, AMR detection complete, analysis run on the
> analysable subset.** See `results/pilot100/PILOT100_QC_REPORT.md` for the
> final counts.

## The selection rule

This is the single most important thing about this analysis, so it is stated
first and enforced in code.

The cohort is:

> **The first 100 genome assemblies physically present in `data/`, sorted by
> GCA accession.**

It is **not** the first 100 rows of `PDC_essential.tsv`.

`PDC_essential.tsv` is opened only *after* the 100 accessions have been
selected, and is used solely to attach metadata to them. The ordering of that
file cannot influence which genomes are analysed.

### The two selections are materially different

| Selection | Count | Overlap with this cohort |
|---|---|---|
| First 100 assemblies in `data/`, GCA order | 100 | — |
| First 100 rows of `PDC_essential.tsv` | 100 | **86 of 100** |

14 of the pilot cohort's assemblies would **not** appear in the PDC
row-order selection. The two are not interchangeable, and
`cohort.prove_selection_is_filesystem_based()` recomputes this comparison on
every run so the claim is auditable rather than asserted.

## Assembly integrity: 45 of the 100 are unusable

The first 100 were selected as specified. When the pipeline tried to read
them, **45 of the 100 turned out to have no usable sequence file**:

| Exclusion reason | Count |
|---|---|
| `corrupt_binary_data` (valid FASTA, then binary garbage) | 28 |
| `undersized` (partial stub, 322 kb - 3.9 Mb) | 16 |
| `zero_byte_file` | 1 |
| **Excluded total** | **45** |
| **Analysable** | **55** |

The 28 corrupt files all break at byte offset 4 194 304, exactly 4 MiB - a
power-of-two block boundary, the signature of a truncated transfer rather
than a FASTA quirk.

This is a property of `data/`, not of the pipeline. Across all 835
assemblies present:

| State | Count |
|---|---|
| Corrupt (binary garbage) | 426 |
| Zero bytes | 233 |
| Intact | 93 |
| Small / partial | 51 |
| Truncated at 4 MiB | 32 |

Only **88 of 835** are intact and of plausible *P. aeruginosa* size.

**Root cause.** The pre-existing `run_pdc.sh` downloaded a `--dehydrated`
NCBI package and then ran `rehydrate`. The rehydrate logged `exit: 0` and
`fna count: 835 / 835`, but that report is wrong: dehydrated packages strip
sequence and rehydration restores it, and for most accessions it evidently
did not complete.

### How the shortfall is handled

The **selection is not changed**. It remains the first 100 assemblies in
`data/` in GCA order. What changes is which of those 100 can be *analysed*:

- each of the 100 is assessed by `cohort.assess_assembly()` and gets
  `analysis_included` TRUE/FALSE plus an `exclusion_reason` in
  `results/pilot100/assembly_integrity.tsv`;
- excluded assemblies are **dropped, never substituted** - the cohort is not
  back-filled to 100 from elsewhere in the accession list;
- every output table, figure and slide reports **n = 55 analysable of 100
  selected**;
- any AMRFinderPlus report belonging to an excluded assembly is deleted
  before the sweep, so a result from a corrupt input can never be reused.

`data/` is still not modified.

## Usage

```bash
# 1. build the cohort and join metadata. Runs NO AMR analysis.
python3 scripts/run_pilot100.py --prepare-only

# 2. run the analysis on the prepared cohort
micromamba run -n pilot100 python scripts/run_pilot100.py --run
```

`--prepare-only` is safe to re-run: it is idempotent and never touches
`data/`.

## Environment

```bash
micromamba create -y -n pilot100 -c conda-forge -c bioconda \
    python=3.11 pyyaml numpy pandas scipy matplotlib-base networkx jinja2 \
    python-pptx pillow pytest snakemake-minimal mlst iqtree snp-sites \
    seqkit blast bakta pyseer
micromamba run -n pilot100 amrfinder_update --database \
    "$(micromamba run -n pilot100 python -c 'import sys,os;print(os.path.dirname(sys.executable))')/share/amrfinderplus"
```

| Tool | Version | Stage |
|---|---|---|
| AMRFinderPlus | 4.2.7 | AMR gene detection |
| AMRFinderPlus database | 2026-08-07.1 | AMR gene detection |
| mlst | 2.33.1 | available |
| IQ-TREE | 3.1.3 | available |
| snp-sites, seqkit, BLAST+ | installed | available |
| Bakta, pyseer, Snakemake | installed | available |

**AMRFinderPlus is not on conda.** The package is `ncbi-amrfinderplus`; it
ships without a database, so `amrfinder_update` must provision one out of
band. Note the v4 CLI rename: `--input` is now `--nucleotide`.

**Panaroo could not be installed on macOS ARM** (its dependency chain
requires a Perl build that does not exist for `osx-arm64`). Stages 9 and 10
are therefore not part of this pilot. The pilot's deliverables do not depend
on them.

## What is produced

### Cohort and metadata

| File | Contents |
|---|---|
| `results/pilot100/first100_assembly_manifest.tsv` | Pilot_ID, Assembly, Genome_path, Genome_file, Genome_found |
| `results/pilot100/first100_assembly_metadata.tsv` | full PDC metadata join, with `PDC_metadata_match` |
| `results/pilot100/assembly_metadata_matching_QC.tsv` | per-assembly PDC match status |
| `results/pilot100/pilot100_manifest.tsv` | final cohort manifest |
| `results/pilot100/genomes/` | 100 symlinks; `data/` is never modified |

### Analysis, per antibiotic (imipenem, meropenem)

| File | Contents |
|---|---|
| `{Antibiotic}_phenotype.tsv` | categorical phenotype per isolate |
| `{Antibiotic}_gene_summary.tsv` | Assembly, Isolate, AST_phenotype, Gene, Mechanism, Detected, Evidence |
| `{Antibiotic}_gene_frequency.tsv` | per-gene counts and frequencies, split by category |
| `{Antibiotic}_mechanism_summary.tsv` | per-mechanism counts, split by category |
| `{Antibiotic}_gene_pairs.tsv` | co-occurring gene pairs, counted once per isolate |
| `{Antibiotic}_mechanism_combinations.tsv` | observed mechanism combinations |

### Figures (`figures/`, 300 DPI)

`pilot100_{Imipenem,Meropenem}_genes.png`,
`pilot100_{Imipenem,Meropenem}_mechanisms.png`,
`pilot100_gene_heatmap.png`, `pilot100_mechanism_heatmap.png`,
`pilot100_{Imipenem,Meropenem}_gene_pairs.png`

A figure is **not** created when the underlying table is too sparse, and the
reason is recorded in `results/pilot100/figure_manifest.tsv`.

### Reports

- `results/pilot100/Pilot100_Pseudomonas_AMR_Analysis.pptx` — 16 slides
- `results/pilot100/PILOT100_QC_REPORT.md`

## Where the biology comes from

No gene-to-mechanism relationship is invented. Three separate, differently
sourced layers are kept distinct and never mixed:

| Layer | Source | Used for |
|---|---|---|
| Gene detection | AMRFinderPlus 4.2.7, database 2026-08-07.1 | which genes are present |
| Functional category | AMRFinderPlus' own `Class`/`Subclass` fields | Efflux, PDC, other beta-lactamase, other AMR |
| ESBL / MBL split | `config/gene_families.tsv` | ESBL, MBL |
| *P. aeruginosa* mechanism | `config/mechanisms.tsv` only | reduced permeability, efflux regulation, AmpC regulation |

### Why `config/gene_families.tsv` exists

AMRFinderPlus reports a single `Class` of `BETA-LACTAM` for every
beta-lactamase, with `Subclass` values of `CEPHALOSPORIN`, `BETA-LACTAM` or
`CARBAPENEM`. **It does not distinguish ESBL from MBL.** Verified on database
2026-08-07.1.

So that split has to come from somewhere. It comes from a reviewable,
editable table using standard beta-lactamase family nomenclature, and the
table says so in its own header. Consequences:

- ESBL and MBL are **nomenclature labels, not enzyme activity measurements**;
- `blaKPC` is deliberately **not** listed as an MBL. KPC is a class A serine
  carbapenemase; calling it a metallo-beta-lactamase would be a factual
  error, so it falls through to "Other beta-lactamase";
- an unlisted beta-lactamase is reported as "Other beta-lactamase", never
  guessed into ESBL or MBL.

### OprD

AMRFinderPlus does not report `oprD`, because a porin is not an AMR gene.
OprD status therefore comes from the PDC genotype field, and is labelled as
such:

| PDC genotype call | Reported as |
|---|---|
| any `oprD*` call | OprD locus reported |
| `oprD_*Ter`, `oprD_*fs*`, indel | OprD disruption → reduced permeability |
| `oprD_V359L` (missense) | variant reported; **not** treated as disruption |
| no `oprD` call | not reported; absence of a call is not evidence of absence |

An `oprD` **gene** hit, were one ever produced, would mean an **intact
locus** — never reduced permeability. This is the rule from
`docs/scientific_rules.md` #1, applied here.

## Limitations

These are substantive and are repeated in the QC report and the deck:

0. **n = 55, not 100.** 45 of the 100 selected assemblies have no usable
   sequence file and are excluded. Every count in this pilot is out of 55.
   This is the single largest limitation and it is a data problem, not an
   analytical one.
1. **Gene detection is not phenotypic resistance.** No isolate is reported as
   resistant because of a gene.
2. **No MIC or zone diameter exists or was derived.** The source provides
   categorical R/I/S only.
3. **No population-structure correction.** No GWAS, phylogenomics or
   convergence analysis was run, so nothing here is corrected for lineage or
   clonal structure. All frequencies are raw counts within this cohort.
4. **Co-occurrence is association only**, unadjusted for clonal structure. A
   pair may co-occur merely because both belong to the same clone.
5. **The cohort is not a random sample.** It is the first 100 of 835
   available assemblies in accession order. Nothing here estimates prevalence
   in the full collection.
6. **ESBL/MBL are nomenclature labels**, not activity measurements.
7. **PDC metadata is used as given**, without independent verification.
8. **Stages 9 and 10 (pan-genome, phylogenomics) were not run** — Panaroo
   is uninstallable on this platform. Lineage structure is therefore
   unknown, so lineage-linked carriage cannot be distinguished from
   resistance-associated carriage.

## Code layout

| Module | Role |
|---|---|
| `papipeline/pilot/cohort.py` | discovery, GCA sorting, selection, PDC join |
| `papipeline/pilot/pdc_fields.py` | AST phenotype and AMR genotype parsers |
| `papipeline/pilot/amr_detect.py` | AMRFinderPlus execution and v4 parsing |
| `papipeline/pilot/mechanisms.py` | category assignment, OprD rules |
| `papipeline/pilot/analysis.py` | gene/mechanism/pair/combination tables |
| `papipeline/pilot/figures.py` | 300-DPI figures with sufficiency checks |
| `papipeline/pilot/pptx.py` | PowerPoint deck |
| `papipeline/pilot/report.py` | QC report |
| `papipeline/pilot/orchestrator.py` | `run_analysis()` |
| `scripts/run_pilot100.py` | CLI: `--prepare-only`, `--run` |

The main pipeline's stage 1 (genome validation) is reused unchanged, as are
its typed record models and its `ClaimStatus` vocabulary.

## Tests

`tests/unit/test_pilot.py` — 81 tests covering selection ordering, the
refusal to substitute when assemblies are short, PDC field parsing
(including ambiguous duplicate keys and invalid categories), mechanism
categorisation, the OprD rules, and the claim ceiling.

One real bug these caught: `oprd_disrupting_genes()` compared
`name.lower().startswith("oprD")`, which never matched, so OprD disruption
was silently reported as absent for every isolate. There is now a
regression test for it.
