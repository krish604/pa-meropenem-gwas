"""The 20-sample cap must stop a real run, not merely be a function.

Ticket 04, wired. `enforce_sample_cap` had thirteen tests and zero callers: the
cap was correct as a unit and absent as a behaviour, which is how hard rule 7
stayed unmet while the ticket read "done".

This test is the evidence the unit tests could not give. Every committed
fixture cohort is exactly 20 samples, so with the 20-sample cap nothing in the
suite could tell "the cap works" apart from "the cohort was never big enough".
The cohort here is 21.

Both directions are asserted, because a cap that fires on everything is not a
cap either: 20 passes, 21 stops. And the 21-sample case must stop *before* any
stage does work.

Written before the wiring.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from papipeline.config.loader import load_config
from papipeline.errors import SampleCapExceeded
from papipeline.run import run_pipeline

PIPELINE_ROOT = Path(__file__).resolve().parents[2]
FIXTURES = PIPELINE_ROOT / "test_data"
OVER_CAP_ID = "TEST_PA_021"


def _clone_fixtures(destination: Path) -> Path:
    """A self-contained project: its own config and its own fixture cohort.

    Both are copied so a run cannot reach back into the repository - a cap that
    were only enforced against the real test_data would not be evidence of
    anything.
    """
    shutil.copytree(PIPELINE_ROOT / "config", destination / "config")
    shutil.copytree(FIXTURES, destination / "test_data")
    return destination / "test_data"


def _add_one_sample(test_data: Path) -> None:
    """Clone sample 001 into sample 021 so the cohort is internally consistent.

    Every per-sample artefact has to exist, because the stages read them: a
    manifest entry with a genome but no annotation table would fail later for an
    unrelated reason and prove nothing about the cap.
    """
    source_id = "TEST_PA_001"
    metadata = test_data / "metadata" / "sample_metadata.tsv"
    lines = metadata.read_text(encoding="utf-8").splitlines()
    header = lines[0]
    template = next(line for line in lines[1:] if line.startswith(source_id + "\t"))
    extra = template.replace(source_id, OVER_CAP_ID)
    metadata.write_text("\n".join(lines + [extra]) + "\n", encoding="utf-8")

    source_genome = test_data / "genomes" / f"{source_id}.fna"
    shutil.copyfile(source_genome, test_data / "genomes" / f"{OVER_CAP_ID}.fna")

    for item in sorted((test_data / "intermediate").rglob("*")):
        if item.is_file() and source_id in item.name:
            target = item.with_name(item.name.replace(source_id, OVER_CAP_ID))
            shutil.copyfile(item, target)

    # Most stage inputs are one file per sample, handled above. The rest are
    # single matrices keyed by sample_id, and those need a row appended rather
    # than a file copied - a missing row fails in a later stage for a reason
    # that has nothing to do with the cap.
    for item in sorted((test_data / "intermediate").rglob("*.tsv")):
        if source_id in item.name:
            continue
        lines = item.read_text(encoding="utf-8").splitlines()
        header = next(
            (line for line in lines if line.startswith("sample_id\t")), None
        )
        if header is None:
            continue
        template = next(
            (line for line in lines if line.split("\t")[0] == source_id), None
        )
        if template is None:
            continue
        if any(line.split("\t")[0] == OVER_CAP_ID for line in lines):
            continue
        item.write_text(
            "\n".join(lines + [template.replace(source_id, OVER_CAP_ID)]) + "\n",
            encoding="utf-8",
        )

    phenotype = test_data / "phenotype" / "imipenem_phenotype.tsv"
    text = phenotype.read_text(encoding="utf-8").splitlines()
    text.append(next(l for l in text[1:] if l.startswith(source_id + "\t")).replace(source_id, OVER_CAP_ID))
    phenotype.write_text("\n".join(text) + "\n", encoding="utf-8")

    # The phylogeny fixtures must GAIN a tip, not swap one. Replacing
    # TEST_PA_001 with TEST_PA_021 would leave the tree short a sample, and the
    # run would then fail in stage 10 for a reason that has nothing to do with
    # the cap under test.
    alignment = test_data / "phylogeny" / "core_alignment.fasta"
    lines = alignment.read_text(encoding="utf-8").splitlines()
    record = [i for i, line in enumerate(lines) if line.startswith(f">{source_id}")]
    start = record[0]
    end = next(
        (i for i in range(start + 1, len(lines)) if lines[i].startswith(">")),
        len(lines),
    )
    cloned = [f">{OVER_CAP_ID}"] + [line.replace(source_id, OVER_CAP_ID) for line in lines[start + 1:end]]
    alignment.write_text("\n".join(lines[:start] + cloned + lines[start:]) + "\n", encoding="utf-8")

    tree = test_data / "phylogeny" / "tree.nwk"
    text = tree.read_text(encoding="utf-8")
    # Graft the new tip as a child of the outermost node. The file's root is
    # labelled and does not end in ")", so anchoring on the final parenthesis
    # would corrupt the string; the outermost "(" is the one stable landmark.
    open_paren = text.index("(")
    text = (
        text[: open_paren + 1]
        + f"{OVER_CAP_ID}:0.05,"
        + text[open_paren + 1:]
    )
    # A leaf tip adds no parentheses, so the string must stay balanced. If this
    # ever fails, the graft corrupted the tree rather than extending it.
    assert text.count("(") == text.count(")"), "newick parentheses must stay balanced"
    tree.write_text(text, encoding="utf-8")

    meta = test_data / "phylogeny" / "tree_metadata.tsv"
    mlines = meta.read_text(encoding="utf-8").splitlines()
    template = next(
        line for line in mlines
        if line and not line.startswith("#") and line.split("\t")[0] == source_id
    )
    meta.write_text("\n".join(mlines + [template.replace(source_id, OVER_CAP_ID)]) + "\n", encoding="utf-8")


def _run(test_data: Path, machine: str):
    config = load_config(
        test_data.parent / "config" / "science.yaml", machine=machine
    )
    return run_pipeline(mode="TEST", config=config)


@pytest.mark.slow
def test_a_cohort_of_twenty_runs_on_the_laptop(tmp_path: Path):
    """The boundary itself must remain usable, or the cap is useless."""
    test_data = _clone_fixtures(tmp_path / "at-cap")
    result = _run(test_data, "laptop")
    assert result.n_samples == 20


@pytest.mark.slow
def test_a_cohort_of_twenty_one_stops_the_run(tmp_path: Path):
    """Hard rule 7: more than 20 samples on the laptop must not proceed."""
    test_data = _clone_fixtures(tmp_path / "over-cap")
    _add_one_sample(test_data)

    with pytest.raises(SampleCapExceeded) as excinfo:
        _run(test_data, "laptop")

    message = str(excinfo.value)
    assert "20" in message
    assert "21" in message
    assert "max_samples" in message


@pytest.mark.slow
def test_the_cap_stops_the_run_before_any_stage_does_work(tmp_path: Path):
    """A refusal after two hours of annotation is not a refusal."""
    test_data = _clone_fixtures(tmp_path / "over-cap-early")
    _add_one_sample(test_data)
    results_root = test_data.parent / "results"

    with pytest.raises(SampleCapExceeded):
        _run(test_data, "laptop")

    stage_files = list(results_root.rglob("stages/*")) if results_root.exists() else []
    assert not stage_files, (
        f"the run produced stage output before refusing: {stage_files[:5]}"
    )


@pytest.mark.slow
def test_the_analysis_machine_accepts_the_same_cohort(tmp_path: Path):
    """The cap is a laptop safety rail, not a property of the science."""
    test_data = _clone_fixtures(tmp_path / "over-cap-big")
    _add_one_sample(test_data)
    result = _run(test_data, "bigmachine")
    assert result.n_samples == 21
