"""Validate workflow shell syntax and the independent scheduled linkcheck."""

from pathlib import Path
import json
import os
import subprocess

import yaml
import pytest


ROOT = Path(__file__).resolve().parents[2]
WORKFLOWS = ROOT / ".github" / "workflows"


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


@pytest.mark.integration
def test_corpus_workflow_propagates_invalid_corpus_and_preserves_report(tmp_path):
    """The actual workflow command must surface runner failures and retain evidence."""
    workflow = yaml.load((WORKFLOWS / "corpus.yml").read_text(), Loader=yaml.BaseLoader)
    job = workflow["jobs"]["corpus"]
    step = next(step for step in job["steps"] if step.get("name") == "Check committed BCR corpus")
    corpus = tmp_path / "invalid corpus"
    corpus.mkdir()
    sentinel = corpus / "source.txt"
    sentinel.write_text("read-only fixture input\n")
    output = tmp_path / "runner output"
    env = {**os.environ, "CORPUS_PATH": str(corpus), "CORPUS_OUTPUT": str(output),
           "MPLCONFIGDIR": str(tmp_path / "matplotlib"), "POLARS_MAX_THREADS": "2"}
    result = subprocess.run(["bash", "-e", "-o", "pipefail", "-c", step["run"]],
                            cwd=ROOT, env=env, capture_output=True, text=True)
    assert result.returncode != 0, result.stdout + result.stderr
    report = json.loads((output / "report.json").read_text())
    assert "manifest" in json.dumps(report).lower()
    assert sentinel.read_text() == "read-only fixture input\n"
