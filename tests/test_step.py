"""Step-by-step pipeline runner: ingest, inspect prompts, run validator gates and record the ledger without a `claude` CLI."""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import pytest
import yaml

from pipeline import step
from pipeline.corpus import load_completed

FX = Path(__file__).parent / "fixtures"
needs_poppler = pytest.mark.skipif(shutil.which("pdftotext") is None, reason="needs poppler")


def _quotes(rec: Path) -> dict[int, list[str]]:
    pages: dict[int, list[str]] = {}

    def walk(x):
        if isinstance(x, dict):
            if isinstance(x.get("page"), int) and isinstance(x.get("quote"), str):
                pages.setdefault(x["page"], []).append(x["quote"])
            for v in x.values():
                walk(v)
        elif isinstance(x, list):
            for v in x:
                walk(v)

    for f in sorted(rec.rglob("*.yaml")):
        walk(yaml.safe_load(f.read_text()))
    return pages


def _pdf(path: Path, pages: dict[int, list[str]]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    path.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(path) as pdf:
        for n in range(1, max(pages) + 1):
            fig = plt.figure(figsize=(8.5, 11))
            for i, q in enumerate(pages.get(n, [])):
                fig.text(0.08, 0.9 - 0.04 * i, q, fontsize=11)
            pdf.savefig(fig)
            plt.close(fig)


@pytest.fixture
def mini_corpus(tmp_path):
    """A self-contained single-paper corpus under `tmp_path` using the toy fixture."""
    root = tmp_path / "proj"
    (root / "corpora" / "demo").mkdir(parents=True)
    (root / ".venv").mkdir()
    (root / ".venv" / "bin").symlink_to(Path(sys.executable).parent)
    q = _quotes(FX / "toy")
    _pdf(root / "corpus" / "demo" / "toy.pdf", q)
    papers = [{"id": "toy", "file": "corpus/demo/toy.pdf", "pages": [1, max(q)]}]
    (root / "corpora" / "demo" / "papers.yaml").write_text(yaml.safe_dump(papers, sort_keys=False))
    return root


@needs_poppler
def test_step_status_blocks_read_until_ingested_and_prompt_emits_constraints(mini_corpus, capsys):
    """Before ingest, `read:toy` is blocked on ingest; after ingest it becomes ready and `prompt` shows the guide, writable paths and gate."""
    rc = step.main(["--root", str(mini_corpus), "--corpus", "demo", "status"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "Needs ingest (1): toy" in out
    assert "read:toy" in out and "needs ingest (toy)" in out

    rc = step.main(["--root", str(mini_corpus), "--corpus", "demo", "ingest"])
    assert rc == 0
    capsys.readouterr()

    rc = step.main(["--root", str(mini_corpus), "--corpus", "demo", "status"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "read:toy" in out and "ready" in out
    assert "audit:toy" in out and "waits for read:toy" in out

    rc = step.main(["--root", str(mini_corpus), "--corpus", "demo", "prompt", "read:toy"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "# BoxGraph Job: read:toy (read)" in out
    assert "record_guide.md" in out
    assert "corpora/demo/records/toy/**" in out
    assert "Read the paper" in out


@needs_poppler
def test_step_gate_checks_soundness_hygiene_fresh_logs_and_records_ledger(mini_corpus, capsys):
    """`gate` refuses paper content inside a record and missing audit logs, then passes and records into `completed.yaml` when clean and sound."""
    assert step.main(["--root", str(mini_corpus), "--corpus", "demo", "ingest"]) == 0
    rec = mini_corpus / "corpora" / "demo" / "records" / "toy"
    shutil.copytree(FX / "toy", rec, dirs_exist_ok=True)
    capsys.readouterr()

    # Paper content leaked directly inside the record must fail the gate.
    stray_png = rec / "screenshot.png"
    stray_png.write_bytes(b"not-a-real-png")
    rc = step.main(["--root", str(mini_corpus), "--corpus", "demo", "gate", "read:toy", "--record"])
    assert rc == 1
    out = capsys.readouterr().out
    assert "screenshot.png" in out and "FAILED" in out
    assert not load_completed(mini_corpus / "corpora" / "demo" / "completed.yaml")

    # Removing the leaked image lets `read:toy` pass and record into the ledger.
    stray_png.unlink()
    rc = step.main(["--root", str(mini_corpus), "--corpus", "demo", "gate", "read:toy", "--record", "--reader", "test-agent"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "read:toy: PASSED" in out
    done = load_completed(mini_corpus / "corpora" / "demo" / "completed.yaml")
    assert "read:toy" in done and done["read:toy"]["reader"] == "test-agent"

    # `audit:toy` must fail without a non-empty audit.md, then pass once written.
    rc = step.main(["--root", str(mini_corpus), "--corpus", "demo", "gate", "audit:toy", "--record"])
    assert rc == 1
    out = capsys.readouterr().out
    assert "audit.md" in out and "FAILED" in out

    (rec / "audit.md").write_text("- Checked toy record values against pages.\n")
    rc = step.main(["--root", str(mini_corpus), "--corpus", "demo", "gate", "audit:toy", "--record"])
    assert rc == 0
    done = load_completed(mini_corpus / "corpora" / "demo" / "completed.yaml")
    assert "audit:toy" in done


def test_find_paper_content_ignores_text_pages_and_figwork(tmp_path):
    """Files inside `text/`, `pages/` and `figwork/` are ignored, while paper/code/data suffixes elsewhere in the record are flagged."""
    rec = tmp_path / "rec"
    for sub in ("text", "pages", "figwork", "models/m1"):
        (rec / sub).mkdir(parents=True)
    (rec / "text" / "page-001.txt").write_text("ok")
    (rec / "pages" / "page-001.png").write_bytes(b"ok")
    (rec / "figwork" / "p1.png").write_bytes(b"ok")
    (rec / "models" / "m1" / "model.py").write_text("def run(i, s): return {}\n")
    assert step.find_paper_content(rec) == []

    (rec / "copied_script.m").write_text("% matlab")
    assert step.find_paper_content(rec) == ["copied_script.m"]
