"""Fold confirmed second-pass KatG scores into the final review ranking."""

from __future__ import annotations

import argparse
import csv
import statistics
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_LIBRARY_DIR = REPO_ROOT / "data" / "ligand_library" / "chembl_approved_v1"
FINAL_COLUMNS = (
    "final_rank", "chembl_id", "pref_name", "final_affinity_kcal_per_mol", "affinity_source",
    "first_pass_rank", "first_pass_affinity_kcal_per_mol", "second_pass_affinity_kcal_per_mol",
    "affinity_change_kcal_per_mol", "docking_percentile", "molecular_weight", "rdkit_logp", "tpsa",
    "hbd", "hba", "rotatable_bonds", "formal_charge", "review_flags", "canonical_smiles",
    "prepared_smiles", "docked_pdbqt_path", "second_pass_pdbqt_path",
)


def final_rows(first_pass: list[dict[str, str]], second_pass: list[dict[str, str]]) -> list[dict[str, str]]:
    confirmed = {
        row["chembl_id"]: row
        for row in second_pass
        if row.get("second_pass_status") == "ok" and row.get("second_pass_best_affinity_kcal_per_mol")
    }
    combined = []
    for row in first_pass:
        first_affinity = float(row["best_affinity_kcal_per_mol"])
        rerun = confirmed.get(row["chembl_id"])
        final_affinity = float(rerun["second_pass_best_affinity_kcal_per_mol"]) if rerun else first_affinity
        combined.append((row, rerun, first_affinity, final_affinity))
    combined.sort(key=lambda item: (item[3], item[0]["chembl_id"]))
    total = len(combined)
    output = []
    for rank, (row, rerun, first_affinity, final_affinity) in enumerate(combined, start=1):
        output.append({
            "final_rank": str(rank), "chembl_id": row["chembl_id"], "pref_name": row.get("pref_name", ""),
            "final_affinity_kcal_per_mol": f"{final_affinity:.3f}",
            "affinity_source": "second_pass_exhaustiveness_8" if rerun else "first_pass_exhaustiveness_4",
            "first_pass_rank": row["rank"], "first_pass_affinity_kcal_per_mol": f"{first_affinity:.3f}",
            "second_pass_affinity_kcal_per_mol": rerun.get("second_pass_best_affinity_kcal_per_mol", "") if rerun else "",
            "affinity_change_kcal_per_mol": f"{final_affinity - first_affinity:.3f}" if rerun else "",
            "docking_percentile": f"{100 * (total - rank + 1) / total:.2f}",
            "molecular_weight": row.get("molecular_weight", ""), "rdkit_logp": row.get("rdkit_logp", ""),
            "tpsa": row.get("tpsa", ""), "hbd": row.get("hbd", ""), "hba": row.get("hba", ""),
            "rotatable_bonds": row.get("rotatable_bonds", ""), "formal_charge": row.get("formal_charge", ""),
            "review_flags": row.get("review_flags", ""), "canonical_smiles": row.get("canonical_smiles", ""),
            "prepared_smiles": row.get("prepared_smiles", ""), "docked_pdbqt_path": row.get("docked_pdbqt_path", ""),
            "second_pass_pdbqt_path": rerun.get("second_pass_pdbqt_path", "") if rerun else "",
        })
    return output


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FINAL_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def finalize(library_dir: Path = DEFAULT_LIBRARY_DIR) -> dict[str, Path | int]:
    ranking_dir = library_dir / "ranking"
    with (ranking_dir / "katg_first_pass_ranked.csv").open(newline="", encoding="utf-8") as handle:
        first_pass = list(csv.DictReader(handle))
    with (ranking_dir / "katg_second_pass_shortlist.csv").open(newline="", encoding="utf-8") as handle:
        second_pass = list(csv.DictReader(handle))
    ranked = final_rows(first_pass, second_pass)
    output_path = ranking_dir / "katg_final_review_ranked.csv"
    summary_path = ranking_dir / "final_review_summary.md"
    write_csv(output_path, ranked)
    rerun_count = sum(row["affinity_source"].startswith("second_pass") for row in ranked)
    rerun_changes = [float(row["affinity_change_kcal_per_mol"]) for row in ranked if row["affinity_change_kcal_per_mol"]]
    summary_path.write_text(
        "# KatG final review ranking summary\n\n"
        f"- Successfully docked ligands ranked: `{len(ranked)}`\n"
        f"- High-exhaustiveness scores incorporated: `{rerun_count}`\n"
        f"- Median second-pass score change: `{statistics.median(rerun_changes):.3f}` kcal/mol\n"
        f"- Best final affinity: `{float(ranked[0]['final_affinity_kcal_per_mol']):.3f}` kcal/mol\n\n"
        "Second-pass scores replace first-pass scores only for the confirmed shortlist. "
        "All other rows retain their first-pass score and are explicitly labelled as such. "
        "Review flags remain context, not exclusion rules or evidence of therapeutic efficacy.\n",
        encoding="utf-8",
    )
    return {"ranked": output_path, "summary": summary_path, "count": len(ranked), "second_pass_incorporated": rerun_count}


def main() -> None:
    parser = argparse.ArgumentParser(description="Create final KatG review ranking")
    parser.add_argument("--library-dir", type=Path, default=DEFAULT_LIBRARY_DIR)
    args = parser.parse_args()
    print(finalize(args.library_dir))


if __name__ == "__main__":
    main()
