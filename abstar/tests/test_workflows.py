"""Execute scheduled workflow shell contracts against bounded local inputs."""

import csv
import json
import os
from pathlib import Path
import shutil
import subprocess

import pytest
import yaml

from .corpus import load_real_bcr_cases


ROOT = Path(__file__).resolve().parents[2]
WORKFLOWS = ROOT / ".github" / "workflows"


def nightly_job():
    return yaml.load((WORKFLOWS / "nightly-corpus.yml").read_text(), Loader=yaml.BaseLoader)["jobs"]["corpus-sentinel"]


def run_step(job, name, cwd, env):
    step = next(step for step in job["steps"] if step.get("name") == name)
    return subprocess.run(["bash", "-e", "-o", "pipefail", "-c", step["run"]],
                          cwd=cwd, env=env, capture_output=True, text=True)


def nightly_env(job, tmp_path):
    # Resolve the runner context exactly as Actions does, including a space.
    context = str(tmp_path / "runner temp")
    return {**os.environ, **{key: value.replace("${{ runner.temp }}", context)
                            for key, value in job["env"].items()},
            "CORPUS_RUN_ID": "123", "CORPUS_ARTIFACT": "tiny", "PER_DATASET": "1"}


@pytest.mark.e2e
def test_nightly_exact_shell_layout_runs_real_tiny_corpus(tmp_path):
    job = nightly_job()
    checkout = tmp_path / "checkout"
    (checkout / "scripts").mkdir(parents=True)
    # _paths must see the same checkout root as the workflow shell.
    shutil.copy2(ROOT / "scripts/discover_bcr_cases.py", checkout / "scripts")
    corpus = checkout / "nightly-corpus"
    (corpus / "bcr_fastas").mkdir(parents=True)
    evidence = corpus / "cellranger" / "001"
    evidence.mkdir(parents=True)
    parent = next(case for case in load_real_bcr_cases() if case.selection_reasons == ("concordant_IGH",))
    (corpus / "bcr_fastas/001.fasta").write_text(f">{parent.sequence_id}\n{parent.sequence}\n")
    (corpus / "sample_manifest.csv").write_text("dataset,donor,flow_class\n001,donor,IgG\n")
    record = {"contig_id": parent.sequence_id, "chain": "IGH", "v_gene": "IGHV3-30", "d_gene": "None",
              "j_gene": "IGHJ4", "c_gene": "IGHM", "productive": "True",
              "cdr3": parent.expected["junction_aa"], "cdr3_nt": parent.expected["junction"], "reads": "1", "umis": "1"}
    with (evidence / "filtered_contig_annotations.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(record))
        writer.writeheader()
        writer.writerow(record)
    env = nightly_env(job, tmp_path)
    for name in ("Initialize outcome report", "Require provisioned corpus configuration",
                 "Validate corpus artifact contract", "Run bounded deterministic cohort"):
        result = run_step(job, name, checkout, env)
        assert result.returncode == 0, result.stdout + result.stderr
    upload = next(step for step in job["steps"] if step.get("name") == "Upload candidate and outcome report")
    path = upload["with"]["path"]
    for key, value in env.items():
        path = path.replace("${{ env." + key + " }}", value)
    output = Path(path)
    assert output.is_absolute() and output.is_relative_to(tmp_path / "runner temp")
    header, row = [json.loads(line) for line in (output / "candidates.jsonl").read_text().splitlines()]
    assert header["selection"]["per_dataset"] == 1
    assert row["source"]["sequence_id"] == parent.sequence_id
    assert row["status"] == "annotated"
    assert row["comparisons"]["productive"] is True
    assert row["comparisons"]["junction_aa"] is True
    assert (output / "analysis.log").is_file()


@pytest.mark.parametrize("quota", ["0", "00", "01", "101", "-1", "1+1", "9" * 100, "$(touch injected)"])
def test_nightly_shell_rejects_noncanonical_or_unbounded_quota(tmp_path, quota):
    job = nightly_job()
    env = nightly_env(job, tmp_path)
    result = run_step(job, "Initialize outcome report", tmp_path, env)
    assert result.returncode == 0, result.stderr
    result = run_step(job, "Require provisioned corpus configuration", tmp_path, {**env, "PER_DATASET": quota})
    assert result.returncode != 0
    assert "canonical decimal integer" in result.stdout
    assert not (tmp_path / "injected").exists()


def test_workflow_shell_syntax_and_expression_boundary():
    for path in WORKFLOWS.glob("*.yml"):
        workflow = yaml.load(path.read_text(), Loader=yaml.BaseLoader)
        for job in workflow["jobs"].values():
            for step in job.get("steps", []):
                if "run" in step:
                    assert "${{" not in step["run"], (path, step.get("name"))
                    result = subprocess.run(["bash", "-n"], input=step["run"], text=True, capture_output=True)
                    assert result.returncode == 0, (path, result.stderr)


def test_scheduled_linkcheck_is_independent_of_corpus_and_push_pr():
    candidates = []
    for path in WORKFLOWS.glob("*.yml"):
        workflow = yaml.load(path.read_text(), Loader=yaml.BaseLoader)
        for job in workflow["jobs"].values():
            if any("-b linkcheck" in step.get("run", "") for step in job.get("steps", [])):
                candidates.append((workflow, job))
    assert len(candidates) == 1
    workflow, job = candidates[0]
    assert "schedule" in workflow["on"]
    assert not {"push", "pull_request", "workflow_call"} & workflow["on"].keys()
    assert not job.get("needs") and not job.get("if")
    assert not any("CORPUS" in str(step) for step in job["steps"])
    assert any("docs/doc_requirements.txt" in step.get("run", "") for step in job["steps"])
