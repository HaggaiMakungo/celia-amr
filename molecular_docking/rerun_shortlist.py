"""Resumable higher-exhaustiveness KatG rerun for the ranked shortlist."""

from __future__ import annotations

import argparse
import csv
import os
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from molecular_docking.batch_docking import (
    DEFAULT_WORKERS,
    _parse_existing_affinity,
    _valid_docked_result,
    prepare_katg_target,
)
from molecular_docking.run_docking import run_vina

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_LIBRARY_DIR = REPO_ROOT / "data" / "ligand_library" / "chembl_approved_v1"
SHORTLIST_RELATIVE_PATH = Path("ranking") / "katg_second_pass_shortlist.csv"
SECOND_PASS_COLUMNS = (
    "second_pass_status",
    "second_pass_pdbqt_path",
    "second_pass_best_affinity_kcal_per_mol",
    "second_pass_attempts",
    "second_pass_notes",
)
RESULTS_DIR_NAME = "second_pass_docking_results"
FAILURES_FILENAME = "second_pass_failures.csv"
SUMMARY_FILENAME = "second_pass_summary.md"
SECOND_PASS_EXHAUSTIVENESS = 8
DEFAULT_TIMEOUT_SECONDS = 30 * 60


def load_shortlist(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or "chembl_id" not in reader.fieldnames:
            raise ValueError("Second-pass shortlist must have a chembl_id column")
        rows, columns = list(reader), list(reader.fieldnames)
    for column in SECOND_PASS_COLUMNS:
        if column not in columns:
            columns.append(column)
        for row in rows:
            row.setdefault(column, "")
    return rows, columns


def write_shortlist(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    temporary = path.with_suffix(".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def _dock_one(index, row, library_dir, results_dir, target, exhaustiveness, timeout_seconds):
    chembl_id = row["chembl_id"]
    ligand = library_dir / "pdbqt" / f"{chembl_id}.pdbqt"
    output = results_dir / f"{chembl_id}.pdbqt"
    relative_output = output.relative_to(library_dir).as_posix()
    if _valid_docked_result(output):
        return index, "ok", relative_output, _parse_existing_affinity(output), "Skipped: valid result file already exists", False
    if not ligand.is_file():
        return index, "failed", "", "", f"Prepared ligand file is missing: {ligand}", True
    try:
        poses = run_vina(
            target.vina_exe, target.receptor_pdbqt, ligand, target.box_center, (20.0, 20.0, 20.0), output,
            seed=42, exhaustiveness=exhaustiveness, cpu=1, timeout_seconds=timeout_seconds,
        )
        if not poses or not _valid_docked_result(output):
            raise RuntimeError("Vina completed without a valid docked pose file")
        return index, "ok", relative_output, f"{poses[0].affinity_kcal_per_mol:.3f}", "", True
    except Exception as exc:
        if output.exists():
            output.unlink()
        return index, "failed", "", "", " ".join(str(exc).split())[:1000], True


def write_failures(path: Path, rows: list[dict[str, str]]) -> None:
    failed = [row for row in rows if row.get("second_pass_status") == "failed"]
    if not failed:
        if path.exists():
            path.unlink()
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=("chembl_id", "pref_name", "second_pass_attempts", "error_message"))
        writer.writeheader()
        for row in failed:
            writer.writerow({
                "chembl_id": row["chembl_id"], "pref_name": row.get("pref_name", ""),
                "second_pass_attempts": row.get("second_pass_attempts", ""), "error_message": row.get("second_pass_notes", ""),
            })


def rerun_shortlist(library_dir: Path = DEFAULT_LIBRARY_DIR, workers: int | None = None, exhaustiveness: int = SECOND_PASS_EXHAUSTIVENESS, timeout_seconds: float | None = DEFAULT_TIMEOUT_SECONDS) -> dict[str, int]:
    shortlist_path = library_dir / SHORTLIST_RELATIVE_PATH
    rows, columns = load_shortlist(shortlist_path)
    candidates = [i for i, row in enumerate(rows) if row.get("second_pass_status") not in {"ok", "failed"}]
    results_dir = library_dir / RESULTS_DIR_NAME
    results_dir.mkdir(parents=True, exist_ok=True)
    worker_count = workers if workers is not None else min(DEFAULT_WORKERS, max(1, (os.cpu_count() or 2) - 2))
    if worker_count < 1 or exhaustiveness < 1:
        raise ValueError("Worker count and exhaustiveness must be at least 1")
    target = prepare_katg_target()
    started = time.monotonic()
    with ThreadPoolExecutor(max_workers=worker_count) as executor:
        futures = {
            executor.submit(_dock_one, index, rows[index], library_dir, results_dir, target, exhaustiveness, timeout_seconds): index
            for index in candidates
        }
        for future in as_completed(futures):
            index, status, result_path, affinity, notes, attempted = future.result()
            row = rows[index]
            row["second_pass_status"] = status
            row["second_pass_pdbqt_path"] = result_path
            row["second_pass_best_affinity_kcal_per_mol"] = affinity
            row["second_pass_notes"] = notes
            if attempted:
                row["second_pass_attempts"] = str(int(row.get("second_pass_attempts") or 0) + 1)
            write_shortlist(shortlist_path, rows, columns)
            write_failures(library_dir / "ranking" / FAILURES_FILENAME, rows)
    counts = {status: sum(row.get("second_pass_status") == status for row in rows) for status in ("ok", "failed")}
    counts["not_attempted"] = len(rows) - counts["ok"] - counts["failed"]
    (library_dir / "ranking" / SUMMARY_FILENAME).write_text(
        "# KatG second-pass docking summary\n\n"
        f"- Shortlist rows: `{len(rows)}`\n"
        f"- Vina exhaustiveness: `{exhaustiveness}`\n"
        f"- Concurrent Vina processes: `{worker_count}`\n"
        f"- Per-ligand timeout: `{('disabled' if timeout_seconds is None else f'{timeout_seconds / 60:.0f} minutes')}`\n"
        f"- Elapsed time: `{time.monotonic() - started:.1f}` seconds\n"
        f"- `ok`: `{counts['ok']}`\n- `failed`: `{counts['failed']}`\n- `not_attempted`: `{counts['not_attempted']}`\n",
        encoding="utf-8",
    )
    return {"selected": len(candidates), **counts}


def main() -> None:
    parser = argparse.ArgumentParser(description="Rerun ranked KatG shortlist at higher Vina exhaustiveness")
    parser.add_argument("--library-dir", type=Path, default=DEFAULT_LIBRARY_DIR)
    parser.add_argument("--workers", type=int)
    parser.add_argument("--exhaustiveness", type=int, default=SECOND_PASS_EXHAUSTIVENESS)
    parser.add_argument("--ligand-timeout-minutes", type=float, default=DEFAULT_TIMEOUT_SECONDS / 60)
    args = parser.parse_args()
    print(rerun_shortlist(args.library_dir, args.workers, args.exhaustiveness, None if args.ligand_timeout_minutes == 0 else args.ligand_timeout_minutes * 60))


if __name__ == "__main__":
    main()
