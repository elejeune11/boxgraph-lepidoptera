"""Run BoxGraph pipeline stages one job at a time without requiring the `claude` CLI.

  python -m pipeline.step --corpus actin_bundling ingest [--papers a,b]
  python -m pipeline.step --corpus actin_bundling status [--papers a,b]
  python -m pipeline.step --corpus actin_bundling prompt read:grason2007
  python -m pipeline.step --corpus actin_bundling gate read:grason2007 [--record]
  python -m pipeline.step --corpus actin_bundling build [--title "..."]

Reuses the exact job graph, prompts, gates, paper-content checks and ledger
(`corpora/<name>/completed.yaml`) of `pipeline.run_papers`, so any agent or
human reviewer can author records step by step while keeping every BoxGraph
validator and hygiene gate intact.
"""

from __future__ import annotations

import argparse
import contextlib
import datetime as dt
import os
import subprocess
import sys
from pathlib import Path

from boxgraph.ingest import ingest as ingest_pdf
from pipeline import run_papers as rp
from pipeline.corpus import (
    check_papers,
    forget_completed,
    layout,
    load_completed,
    load_papers,
    record_completed,
)
from pipeline.orchestrator import PAPER_CONTENT, RECORD_IGNORED, Job, child_env
from pipeline.pairs import candidate_pairs, split_pairs

ROOT = Path(__file__).resolve().parents[1]


@contextlib.contextmanager
def _scoped(root: Path, corpus: str | None):
    """Temporarily point `pipeline.run_papers` at `(root, corpus)` and restore afterwards."""
    saved = {k: getattr(rp, k) for k in ("ROOT", "PAPERS", "LEDGER", "OVERVIEW", "RECORDS", "BUILD", "STATE")}
    try:
        rp.ROOT = Path(root).resolve()
        rp.use_corpus(corpus)
        yield layout(rp.ROOT, corpus)
    finally:
        for k, v in saved.items():
            setattr(rp, k, v)


def _paper_ids(job_id: str) -> list[str]:
    return job_id.split(":", 1)[1].split("+")


def _ready_papers(ps: list[dict], done: dict[str, dict]) -> list[str]:
    return [e["id"] for e in ps if all(f"{s}:{e['id']}" in done for s in rp.stages_of(e))]


def _all_jobs(
    root: Path,
    corpus: str | None,
    *,
    model: str | None = None,
    extra: str = "",
) -> tuple[list[dict], list[Job], dict[tuple[str, str], list[str]], dict[tuple[str, str], list[str]], dict[str, dict]]:
    """Build every per-paper, relation and validation job currently defined for the corpus."""
    with _scoped(root, corpus) as lay:
        ps = load_papers(lay.papers)
        done = load_completed(lay.ledger)
        per_map = rp.jobs(model, 40.0, 20.0, 15.0, extra=extra, budget_digitize=30.0)
        ordered: list[Job] = [*per_map["read"], *per_map["audit"], *per_map["digitize"]]
        manual = {p["id"]: [*(p.get("relate_with") or []), *(p.get("papers") or [])] for p in ps}
        ready = _ready_papers(ps, done)
        pairs, shape_only = split_pairs(candidate_pairs(rp.ROOT / lay.records, ready, manual))
        rel = rp.relate_jobs(pairs, model, 15.0, extra=extra)
        related_pairs = {k: v for k, v in pairs.items() if f"relate:{k[0]}+{k[1]}" in done}
        val = rp.validate_jobs(related_pairs, model, 15.0, extra=extra)
        ordered.extend(rel)
        ordered.extend(val)
        return ps, ordered, pairs, shape_only, done


def _job_state(job: Job, done: dict[str, dict], needs_ingest: set[str]) -> tuple[str, str]:
    if job.id in done:
        return "done", f"done {done[job.id].get('date', '?')}"
    if job.kind == "read":
        pid = _paper_ids(job.id)[0]
        if pid in needs_ingest:
            return "blocked", f"needs ingest ({pid})"
    missing_deps = [d for d in job.depends_on if d not in done]
    if missing_deps:
        return "blocked", "waits for " + ", ".join(missing_deps)
    return "ready", "ready"


def find_paper_content(record_dir: Path) -> list[str]:
    """Any file inside `record_dir` (outside text/, pages/, figwork/) with a paper/code/data suffix."""
    bad: list[str] = []
    if not record_dir.is_dir():
        return bad
    for dirpath, dirs, files in os.walk(record_dir):
        rel = Path(dirpath).relative_to(record_dir)
        if rel == Path("."):
            dirs[:] = [d for d in dirs if d not in RECORD_IGNORED]
        for f in sorted(files):
            p = rel / f
            if p.suffix.lower() in PAPER_CONTENT:
                bad.append(p.as_posix())
    return bad


def run_ingest(root: Path, corpus: str | None, only: list[str] | None = None) -> tuple[list[str], list[str]]:
    """Ingest any new or changed PDFs in the corpus list; returns `(ingested, failed)`."""
    with _scoped(root, corpus) as lay:
        ps = load_papers(lay.papers)
        done = load_completed(lay.ledger)
        ingested: list[str] = []
        failed: list[str] = []
        for p in ps:
            if only is not None and p["id"] not in only:
                continue
            if p.get("kind") in ("software", "dataset") or not rp._ingest_needed(p):
                continue
            try:
                ingest_pdf(rp.ROOT / p["file"], p["id"], rp.ROOT / lay.records, pages=tuple(p["pages"]))
            except Exception as e:  # noqa: BLE001
                failed.append(f"{p['id']}: {e}")
                continue
            stale = rp._stale_rows(done, p["id"])
            if stale:
                forget_completed(lay.ledger, stale)
                done = load_completed(lay.ledger)
            ingested.append(p["id"])
        return ingested, failed


def run_gate(root: Path, corpus: str | None, job: Job, *, record: bool = False, reader: str | None = None) -> tuple[bool, str]:
    """Run a job's exact BoxGraph validator gate plus paper-content and fresh-log hygiene checks."""
    lay = layout(Path(root).resolve(), corpus)
    issues: list[str] = []

    # 1. Paper-content hygiene check on every target record
    for rid in _paper_ids(job.id):
        rec_dir = Path(root).resolve() / lay.records / rid
        leaked = find_paper_content(rec_dir)
        if leaked:
            issues.append(
                f"paper/code/data content found in {lay.records}/{rid}: {', '.join(leaked)} "
                "(remove or move scratch images to figwork/)"
            )

    # 2. Required log files (audit.md, figures.md) must exist and be non-empty
    for rel_path in job.fresh:
        fp = Path(job.cwd) / rel_path
        if not fp.is_file() or fp.stat().st_size == 0:
            issues.append(f"required log missing or empty: {rel_path}")

    if issues:
        return False, "\n".join(f"ERROR   {m}" for m in issues)

    # 3. Exact validator command defined by `pipeline.run_papers`
    env = child_env() | job.env
    proc = subprocess.run(
        job.validate,
        cwd=job.cwd,
        env=env,
        capture_output=True,
        text=True,
        timeout=job.validate_timeout_s,
    )
    out = (proc.stdout + proc.stderr).strip()
    ok = proc.returncode == 0
    if ok and record:
        row: dict = {
            "job": job.id,
            "date": dt.date.today().isoformat(),
            "cost_usd": 0.0,
            "turns": 0,
        }
        if reader:
            row["reader"] = reader
        record_completed(lay.ledger, [row])
    return ok, out


def _cmd_ingest(a) -> int:
    root = Path(a.root).resolve()
    lay = layout(root, a.corpus)
    ps = load_papers(lay.papers)
    errs = check_papers(ps, root)
    if errs:
        print(f"{lay.papers} has problems:", *errs, sep="\n  ")
        return 2
    only = a.papers.split(",") if a.papers else None
    ingested, failed = run_ingest(root, a.corpus, only=only)
    for pid in ingested:
        print(f"ingested {pid} -> {lay.records}/{pid}")
    if not ingested and not failed:
        print("nothing to ingest (all PDFs up to date)")
    for err in failed:
        print(f"ERROR   {err}")
    return 1 if failed else 0


def _cmd_status(a) -> int:
    root = Path(a.root).resolve()
    lay = layout(root, a.corpus)
    ps = load_papers(lay.papers)
    errs = check_papers(ps, root)
    if errs:
        print(f"{lay.papers} has problems:", *errs, sep="\n  ")
        return 2
    with _scoped(root, a.corpus):
        needs_ingest = {
            p["id"]
            for p in ps
            if p.get("kind") not in ("software", "dataset") and rp._ingest_needed(p)
        }
    only = set(a.papers.split(",")) if a.papers else None
    _, all_js, pairs, shape_only, done = _all_jobs(root, a.corpus)
    if only is not None:
        all_js = [j for j in all_js if set(_paper_ids(j.id)) & only]
    if needs_ingest:
        todo_ing = sorted(needs_ingest & only) if only is not None else sorted(needs_ingest)
        if todo_ing:
            print(f"Needs ingest ({len(todo_ing)}): {', '.join(todo_ing)}")
            print(f"  run: .venv/bin/python -m pipeline.step --corpus {a.corpus} ingest")
            print()
    counts = {"done": 0, "ready": 0, "blocked": 0}
    ready_ids: list[str] = []
    for j in all_js:
        st, detail = _job_state(j, done, needs_ingest)
        counts[st] += 1
        if st == "ready":
            ready_ids.append(j.id)
        print(f"  {j.id:36} {st:8} {detail}")
    if shape_only:
        print(f"\nShape-only pairs (listed for a person, not scheduled): {len(shape_only)}")
        for (x, y), why in sorted(shape_only.items()):
            print(f"  {x} + {y}: " + "; ".join(why))
    print(
        f"\nSummary: {counts['done']} done, {counts['ready']} ready, {counts['blocked']} blocked "
        f"(corpus: {a.corpus}, {len(ps)} source(s))"
    )
    if ready_ids:
        print(f"Next: .venv/bin/python -m pipeline.step --corpus {a.corpus} prompt {ready_ids[0]}")
    return 0


def _cmd_prompt(a) -> int:
    root = Path(a.root).resolve()
    lay = layout(root, a.corpus)
    ps = load_papers(lay.papers)
    errs = check_papers(ps, root)
    if errs:
        print(f"{lay.papers} has problems:", *errs, sep="\n  ")
        return 2
    extra = Path(a.extra).read_text() if a.extra and Path(a.extra).is_file() else (a.extra or "")
    _, all_js, _, _, done = _all_jobs(root, a.corpus, model=a.model, extra=extra)
    by_id = {j.id: j for j in all_js}
    job = by_id.get(a.job)
    if job is None:
        print(f"unknown or not-yet-unlocked job {a.job!r}; run `status` to list available jobs")
        return 2
    missing_deps = [d for d in job.depends_on if d not in done]
    print(f"# BoxGraph Job: {job.id} ({job.kind})")
    if missing_deps:
        print(f"# WARNING: dependencies not yet recorded in ledger: {', '.join(missing_deps)}")
    print(f"# Guide:    {job.append_system_prompt_file}")
    print(f"# Writable: {', '.join(job.writable)}")
    if job.fresh:
        print(f"# Fresh:    {', '.join(job.fresh)}")
    print(f"# Gate:     {' '.join(job.validate)}")
    print("-" * 72)
    print(job.prompt)
    return 0


def _cmd_gate(a) -> int:
    root = Path(a.root).resolve()
    lay = layout(root, a.corpus)
    ps = load_papers(lay.papers)
    errs = check_papers(ps, root)
    if errs:
        print(f"{lay.papers} has problems:", *errs, sep="\n  ")
        return 2
    _, all_js, _, _, done = _all_jobs(root, a.corpus)
    by_id = {j.id: j for j in all_js}
    job = by_id.get(a.job)
    if job is None:
        print(f"unknown or not-yet-unlocked job {a.job!r}; run `status` to list available jobs")
        return 2
    missing_deps = [d for d in job.depends_on if d not in done]
    if missing_deps and not a.allow_unmet_deps:
        print(f"ERROR   cannot pass {job.id}: dependencies not in ledger: {', '.join(missing_deps)}")
        return 1
    ok, out = run_gate(root, a.corpus, job, record=a.record, reader=a.reader)
    if out:
        print(out)
    status = "PASSED" if ok else "FAILED"
    recorded = f" (recorded in {lay.ledger.relative_to(root).as_posix()})" if (ok and a.record) else ""
    print(f"{job.id}: {status}{recorded}")
    return 0 if ok else 1


def _cmd_build(a) -> int:
    root = Path(a.root).resolve()
    lay = layout(root, a.corpus)
    ps = load_papers(lay.papers)
    title = a.title or f"BoxGraph Synthesis: {len(ps)} sources ({a.corpus})"
    cmd = [
        str(root / ".venv" / "bin" / "boxgraph"),
        "build",
        lay.records,
        "--out",
        lay.build,
        "--cache",
        rp.CACHE,
        "--save-readings",
        "--title",
        title,
    ]
    if lay.overview.exists():
        cmd += ["--overview", lay.overview.relative_to(root).as_posix()]
    return subprocess.run(cmd, cwd=root, check=False).returncode


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--corpus", required=True, help="the corpus, kept in corpora/<name>/")
    ap.add_argument("--root", default=str(ROOT), help=argparse.SUPPRESS)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_ing = sub.add_parser("ingest", help="ingest new or changed PDFs in the corpus list")
    p_ing.add_argument("--papers", help="comma-separated paper ids")
    p_ing.set_defaults(fn=_cmd_ingest)

    p_st = sub.add_parser("status", help="show done, ready and blocked jobs")
    p_st.add_argument("--papers", help="comma-separated paper ids")
    p_st.set_defaults(fn=_cmd_status)

    p_pr = sub.add_parser("prompt", help="print the exact BoxGraph prompt and constraints for one job")
    p_pr.add_argument("job", help="job id, e.g. read:grason2007")
    p_pr.add_argument("--model", help="reader model label for prompt metadata")
    p_pr.add_argument("--extra", help="extra instruction text or path appended to prompt")
    p_pr.set_defaults(fn=_cmd_prompt)

    p_gt = sub.add_parser("gate", help="run the exact BoxGraph validator gate for one job")
    p_gt.add_argument("job", help="job id, e.g. read:grason2007")
    p_gt.add_argument("--record", action="store_true", help="record a passing job in corpora/<name>/completed.yaml")
    p_gt.add_argument("--reader", help="optional reader identifier stored with the ledger entry")
    p_gt.add_argument("--allow-unmet-deps", action="store_true", help="run gate even if dependencies are not yet recorded")
    p_gt.set_defaults(fn=_cmd_gate)

    p_bd = sub.add_parser("build", help="build the Synthesis, report, figures and overview for the corpus")
    p_bd.add_argument("--title", help="title of the Synthesis")
    p_bd.set_defaults(fn=_cmd_build)

    a = ap.parse_args(argv)
    try:
        lay = layout(Path(a.root).resolve(), a.corpus)
    except ValueError as e:
        print(e)
        return 2
    if not lay.papers.exists():
        print(f"{lay.papers} does not exist: add a paper with python -m pipeline.add --corpus {a.corpus}")
        return 2
    return a.fn(a)


if __name__ == "__main__":
    sys.exit(main())
