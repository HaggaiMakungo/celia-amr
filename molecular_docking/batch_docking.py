"""Resumable, manifest-driven first-pass batch docking against a TB target.

Gene-parameterized: pass ``--gene katG`` (default) or ``--gene rpoB`` to pick
the receptor/target from ``run_docking.TARGETS``. Output filenames and
receptor prep are derived from the gene so a KatG run and an rpoB run never
collide in the same library directory.

The default is a safe 24-ligand timing/health check. Use ``--full`` only
after reviewing that test batch. Each Vina process gets one CPU; concurrency
is controlled by ``--workers`` so the machine remains responsive.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

from molecular_docking.run_docking import (
    TARGETS,
    DockingTarget,
    clean_receptor_pdb,
    compute_box_center,
    download_pdb_structure,
    ensure_vina_binary,
    prepare_receptor_pdbqt,
    run_vina,
)

MODULE_DIR = Path(__file__).resolve().parent
REPO_ROOT = MODULE_DIR.parent
DEFAULT_LIBRARY_DIR = REPO_ROOT / "data" / "ligand_library" / "chembl_approved_v1"
DEFAULT_GENE = "katG"
FIRST_PASS_EXHAUSTIVENESS = 4


def _gene_key(gene: str) -> str:
    """Lowercase filename-safe key, e.g. 'katG' -> 'katg', 'rpoB' -> 'rpob'."""
    return gene.lower()


def results_dir_name(gene: str) -> str:
    # KatG keeps its original unprefixed name so the already-completed KatG
    # run's output directory is still found by a plain re-run; every other
    # gene gets a prefixed name so runs never collide in one library dir.
    return "docking_results" if _gene_key(gene) == "katg" else f"{_gene_key(gene)}_docking_results"


def failures_filename(gene: str) -> str:
    return "docking_failures.csv" if _gene_key(gene) == "katg" else f"{_gene_key(gene)}_docking_failures.csv"


def summary_filename(gene: str) -> str:
    return "docking_summary.md" if _gene_key(gene) == "katg" else f"{_gene_key(gene)}_docking_summary.md"


DEFAULT_TEST_BATCH_SIZE = 24
DEFAULT_WORKERS = 6  # i7-9750H physical-core count; avoids SMT oversubscription.
DEFAULT_LIGAND_TIMEOUT_SECONDS = 15 * 60
DOCKING_COLUMNS = ("docking_status", "docked_pdbqt_path", "best_affinity_kcal_per_mol", "docking_attempts", "docking_notes")
VINA_RESULT = re.compile(r"^REMARK VINA RESULT:\s+(-?\d+(?:\.\d+)?)", re.MULTILINE)


@dataclass(frozen=True)
class PreparedTarget:
    target: DockingTarget
    vina_exe: Path
    receptor_pdbqt: Path
    box_center: tuple[float, float, float]


@dataclass(frozen=True)
class DockingOutcome:
    row_index: int
    status: str
    result_path: str = ""
    affinity: str = ""
    notes: str = ""
    attempted: bool = False
    elapsed_seconds: float = 0.0


def _parse_existing_affinity(path: Path) -> str:
    match = VINA_RESULT.search(path.read_text(errors="replace"))
    return match.group(1) if match else ""


def _valid_docked_result(path: Path) -> bool:
    return path.exists() and path.stat().st_size > 0 and "MODEL" in path.read_text(errors="replace")


def load_manifest(manifest_path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with manifest_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"Manifest has no header: {manifest_path}")
        missing = {"chembl_id", "prep_status", "pdbqt_path"}.difference(reader.fieldnames)
        if missing:
            raise ValueError(f"Manifest is missing required columns: {', '.join(sorted(missing))}")
        rows = [dict(row) for row in reader]
        columns = list(reader.fieldnames)
    for column in DOCKING_COLUMNS:
        if column not in columns:
            columns.append(column)
        for row in rows:
            row.setdefault(column, "")
    return rows, columns


def write_manifest(manifest_path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    temporary_path = manifest_path.with_suffix(".tmp")
    with temporary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    temporary_path.replace(manifest_path)


def eligible_row_indices(rows: list[dict[str, str]]) -> list[int]:
    """Rows suitable for a fresh docking attempt (never retry failures implicitly)."""
    return [
        index
        for index, row in enumerate(rows)
        if row.get("prep_status") == "ok" and row.get("docking_status", "not_attempted") not in {"ok", "failed"}
    ]


def select_test_batch(indices: list[int], rows: list[dict[str, str]], size: int) -> list[int]:
    """Select a deterministic spread across molecular weights for the test run."""
    if size < 1:
        raise ValueError("Test batch size must be at least 1")
    ordered = sorted(indices, key=lambda index: float(rows[index].get("molecular_weight") or 0.0))
    if len(ordered) <= size:
        return ordered
    positions = {round(position * (len(ordered) - 1) / (size - 1)) for position in range(size)} if size > 1 else {0}
    return [ordered[position] for position in sorted(positions)]


def prepare_target(gene: str) -> PreparedTarget:
    """Create the receptor/grid for ``gene``, the same way as the validated one-off run."""
    target = TARGETS[gene]
    if not target.receptor_pdb_path.exists():
        download_pdb_structure(target.pdb_id, target.receptor_pdb_path)
    vina_exe = ensure_vina_binary()
    box_center = compute_box_center(target.receptor_pdb_path, target.box_ligand_resname, target.box_chain)
    work_dir = MODULE_DIR / "work"
    cleaned_pdb = work_dir / f"{target.pdb_id}_cleaned.pdb"
    receptor_pdbqt = work_dir / f"{target.pdb_id}_receptor.pdbqt"
    # Receptor preparation is intentionally the same as the validated one-off path.
    clean_receptor_pdb(target.receptor_pdb_path, target.strip_resnames, cleaned_pdb)
    prepare_receptor_pdbqt(cleaned_pdb, target.bonds_to_delete, receptor_pdbqt)
    return PreparedTarget(target, vina_exe, receptor_pdbqt, box_center)


# Backward-compatible alias: existing call sites/imports (e.g. rerun_shortlist.py)
# may still reference the original KatG-only name.
def prepare_katg_target() -> PreparedTarget:
    return prepare_target("katG")


def _dock_one(
    row_index: int,
    row: dict[str, str],
    library_dir: Path,
    results_dir: Path,
    prepared_target: PreparedTarget,
    exhaustiveness: int,
    seed: int,
    timeout_seconds: float | None,
) -> DockingOutcome:
    chembl_id = row["chembl_id"]
    ligand_path = library_dir / row["pdbqt_path"]
    result_path = results_dir / f"{chembl_id}.pdbqt"
    relative_result_path = result_path.relative_to(library_dir).as_posix()
    if _valid_docked_result(result_path):
        return DockingOutcome(
            row_index=row_index,
            status="ok",
            result_path=relative_result_path,
            affinity=_parse_existing_affinity(result_path),
            notes="Skipped: valid result file already exists",
        )
    if not ligand_path.is_file():
        return DockingOutcome(row_index, "failed", notes=f"Prepared ligand file is missing: {ligand_path}", attempted=True)

    started = time.monotonic()
    try:
        poses = run_vina(
            prepared_target.vina_exe,
            prepared_target.receptor_pdbqt,
            ligand_path,
            prepared_target.box_center,
            (20.0, 20.0, 20.0),
            result_path,
            seed=seed,
            exhaustiveness=exhaustiveness,
            cpu=1,
            timeout_seconds=timeout_seconds,
        )
        if not poses or not _valid_docked_result(result_path):
            raise RuntimeError("Vina completed without a valid docked pose file")
        return DockingOutcome(
            row_index=row_index,
            status="ok",
            result_path=relative_result_path,
            affinity=f"{poses[0].affinity_kcal_per_mol:.3f}",
            attempted=True,
            elapsed_seconds=time.monotonic() - started,
        )
    except Exception as exc:
        if result_path.exists():
            result_path.unlink()
        message = " ".join(str(exc).split())[:1000]
        return DockingOutcome(
            row_index=row_index,
            status="failed",
            notes=message,
            attempted=True,
            elapsed_seconds=time.monotonic() - started,
        )


def write_failures(path: Path, rows: list[dict[str, str]]) -> None:
    failures = [row for row in rows if row.get("prep_status") == "ok" and row.get("docking_status") == "failed"]
    if not failures:
        if path.exists():
            path.unlink()
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=("chembl_id", "pref_name", "docking_attempts", "error_message"))
        writer.writeheader()
        for row in failures:
            writer.writerow(
                {
                    "chembl_id": row["chembl_id"],
                    "pref_name": row.get("pref_name", ""),
                    "docking_attempts": row.get("docking_attempts", ""),
                    "error_message": row.get("docking_notes", ""),
                }
            )


def write_summary(path: Path, rows: list[dict[str, str]], mode: str, elapsed_seconds: float, workers: int, gene: str = "katG") -> None:
    ready_rows = [row for row in rows if row.get("prep_status") == "ok"]
    status_counts: dict[str, int] = {}
    for row in ready_rows:
        status = row.get("docking_status") or "not_attempted"
        status_counts[status] = status_counts.get(status, 0) + 1
    failure_patterns: dict[str, int] = {}
    for row in ready_rows:
        if row.get("docking_status") == "failed":
            message = row.get("docking_notes", "")
            failure_patterns[message] = failure_patterns.get(message, 0) + 1
    lines = [
        f"# {gene} batch docking summary",
        "",
        f"- Mode: `{mode}`",
        f"- First-pass Vina exhaustiveness: `{FIRST_PASS_EXHAUSTIVENESS}`",
        f"- Concurrent Vina processes: `{workers}` (one CPU each)",
        f"- Elapsed time: `{elapsed_seconds:.1f}` seconds",
        f"- Ready library rows: `{len(ready_rows)}`",
    ]
    for status in ("ok", "failed", "not_attempted"):
        lines.append(f"- `{status}`: `{status_counts.get(status, 0)}`")
    if failure_patterns:
        lines.extend(["", "## Failure patterns", ""])
        for message, count in sorted(failure_patterns.items(), key=lambda item: (-item[1], item[0])):
            lines.append(f"- `{count}`: {message}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_batch(
    library_dir: Path,
    gene: str = DEFAULT_GENE,
    test_batch_size: int | None = DEFAULT_TEST_BATCH_SIZE,
    workers: int | None = None,
    exhaustiveness: int = FIRST_PASS_EXHAUSTIVENESS,
    seed: int = 42,
    ligand_timeout_seconds: float | None = DEFAULT_LIGAND_TIMEOUT_SECONDS,
) -> dict[str, int]:
    if exhaustiveness < 1:
        raise ValueError("Exhaustiveness must be at least 1")
    if gene not in TARGETS:
        raise ValueError(f"Unknown gene {gene!r}; choices are {sorted(TARGETS)}")
    manifest_path = library_dir / "ligand_library_manifest.csv"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Prepared-library manifest not found: {manifest_path}")
    rows, columns = load_manifest(manifest_path)
    candidates = eligible_row_indices(rows)
    if test_batch_size is not None:
        candidates = select_test_batch(candidates, rows, test_batch_size)
        mode = f"test batch ({len(candidates)} ligands)"
    else:
        mode = "full library"
    results_dir = library_dir / results_dir_name(gene)
    results_dir.mkdir(parents=True, exist_ok=True)
    worker_count = workers if workers is not None else min(DEFAULT_WORKERS, max(1, (os.cpu_count() or 2) - 2))
    if worker_count < 1:
        raise ValueError("Worker count must be at least 1")
    prepared_target = prepare_target(gene)
    started = time.monotonic()
    with ThreadPoolExecutor(max_workers=worker_count) as executor:
        futures = {
            executor.submit(
                _dock_one,
                index,
                rows[index],
                library_dir,
                results_dir,
                prepared_target,
                exhaustiveness,
                seed,
                ligand_timeout_seconds,
            ): index
            for index in candidates
        }
        for future in as_completed(futures):
            outcome = future.result()
            row = rows[outcome.row_index]
            row["docking_status"] = outcome.status
            row["docked_pdbqt_path"] = outcome.result_path
            row["best_affinity_kcal_per_mol"] = outcome.affinity
            row["docking_notes"] = outcome.notes
            if outcome.attempted:
                row["docking_attempts"] = str(int(row.get("docking_attempts") or 0) + 1)
            # Persist after each result: a restart sees completed output and does not redo it.
            write_manifest(manifest_path, rows, columns)
            write_failures(library_dir / failures_filename(gene), rows)
    write_summary(library_dir / summary_filename(gene), rows, mode, time.monotonic() - started, worker_count, gene=gene)
    return {
        "selected": len(candidates),
        "ok": sum(row.get("docking_status") == "ok" for row in rows if row.get("prep_status") == "ok"),
        "failed": sum(row.get("docking_status") == "failed" for row in rows if row.get("prep_status") == "ok"),
        "not_attempted": sum(not row.get("docking_status") or row.get("docking_status") == "not_attempted" for row in rows if row.get("prep_status") == "ok"),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Resumable first-pass batch docking against a TB target")
    parser.add_argument("--gene", default=DEFAULT_GENE, choices=sorted(TARGETS.keys()))
    parser.add_argument("--library-dir", type=Path, default=DEFAULT_LIBRARY_DIR)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--full", action="store_true", help="Dock every not-yet-attempted prepared ligand")
    mode.add_argument("--test-batch", type=int, nargs="?", const=DEFAULT_TEST_BATCH_SIZE, default=DEFAULT_TEST_BATCH_SIZE)
    parser.add_argument("--workers", type=int, help="Concurrent Vina processes (default: six physical cores on this machine)")
    parser.add_argument("--exhaustiveness", type=int, default=FIRST_PASS_EXHAUSTIVENESS)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--ligand-timeout-minutes",
        type=float,
        default=DEFAULT_LIGAND_TIMEOUT_SECONDS / 60,
        help="Per-ligand Vina limit; use 0 to disable (default: 15)",
    )
    args = parser.parse_args()
    summary = run_batch(
        args.library_dir,
        gene=args.gene,
        test_batch_size=None if args.full else args.test_batch,
        workers=args.workers,
        exhaustiveness=args.exhaustiveness,
        seed=args.seed,
        ligand_timeout_seconds=None if args.ligand_timeout_minutes == 0 else args.ligand_timeout_minutes * 60,
    )
    print(summary)


if __name__ == "__main__":
    main()
