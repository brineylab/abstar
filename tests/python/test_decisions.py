import re
from pathlib import Path

DIR = Path(__file__).parents[2] / "python/abstar/decisions"
NAME = re.compile(r"^(\d{4})-[a-z0-9-]+\.md$")
HEADINGS = ["## Status", "## Context", "## Decision", "## Alternatives considered", "## Consequences"]
STATUS = re.compile(r"^(Accepted|Superseded by \d{4}) \(\d{4}-\d{2}-\d{2}\)$")


def _records():
    return sorted(p for p in DIR.glob("*.md") if p.name != "README.md")


def test_records_are_well_formed():
    recs = _records()
    assert recs, "no decision records found"
    nums = []
    for p in recs:
        m = NAME.match(p.name)
        assert m, p.name
        nums.append(int(m.group(1)))
        lines = p.read_text().splitlines()
        assert re.match(rf"^# {m.group(1)}: \S", lines[0]), p.name
        heads = [l for l in lines if l.startswith("## ")]
        assert heads == HEADINGS, (p.name, heads)
        i = lines.index("## Status")
        status = next(l for l in lines[i + 1:] if l.strip())
        assert STATUS.match(status), (p.name, status)
    assert nums == list(range(1, len(nums) + 1))


def test_seed_records_present():
    nums = {int(p.name[:4]) for p in _records()}
    assert set(range(1, 17)) <= nums
