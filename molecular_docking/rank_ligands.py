"""Create an auditable first-pass ranking from completed docking results.

Gene-parameterized: pass ``--gene katG`` (default) or ``--gene rpoB``. Reads
docking columns written by ``batch_docking.py`` for that gene and writes
gene-prefixed ranking files, so KatG's original filenames are unchanged.

Rows are ordered by Vina affinity alone. RDKit-derived properties are review
flags, not an unvalidated composite efficacy score.
"""

from __future__ import annotations

import argparse
import csv
import statistics
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_LIBRARY_DIR = REPO_ROOT / "data" / "ligand_library" / "chembl_approved_v1"
DEFAULT_GENE = "katG"


def _gene_key(gene: str) -> str:
    return gene.lower()


def ranked_filename(gene: str) -> str:
    return f"{_gene_key(gene)}_first_pass_ranked.csv"


def shortlist_filename(gene: str) -> str:
    return f"{_gene_key(gene)}_second_pass_shortlist.csv"


def summary_filename(gene: str) -> str:
    # KatG keeps its original unprefixed name to match the already-published run.
    return "ranking_summary.md" if _gene_key(gene) == "katg" else f"{_gene_key(gene)}_ranking_summary.md"


RANKING_COLUMNS = (
    "rank", "chembl_id", "pref_name", "best_affinity_kcal_per_mol", "docking_percentile",
    "molecular_weight", "rdkit_logp", "tpsa", "hbd", "hba", "rotatable_bonds", "formal_charge",
    "review_flags", "second_pass_recommendation", "canonical_smiles", "prepared_smiles", "docked_pdbqt_path",
)


def property_profile(smiles: str) -> dict[str, float | int]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("RDKit could not parse prepared SMILES")
    return {
        "molecular_weight": Descriptors.MolWt(mol),
        "rdkit_logp": Crippen.MolLogP(mol),
        "tpsa": rdMolDescriptors.CalcTPSA(mol),
        "hbd": Lipinski.NumHDonors(mol),
        "hba": Lipinski.NumHAcceptors(mol),
        "rotatable_bonds": Lipinski.NumRotatableBonds(mol),
        "formal_charge": Chem.GetFormalCharge(mol),
    }


def review_flags(profile: dict[str, float | int]) -> list[str]:
    flags = []
    if profile["molecular_weight"] > 600:
        flags.append("high_molecular_weight")
    if profile["rotatable_bonds"] > 12:
        flags.append("high_flexibility")
    if abs(profile["formal_charge"]) >= 2:
        flags.append("high_formal_charge")
    if profile["rdkit_logp"] > 5:
        flags.append("high_logp")
    if profile["tpsa"] > 140:
        flags.append("high_polar_surface_area")
    return flags


def ranked_rows(manifest_rows: list[dict[str, str]], shortlist_size: int) -> list[dict[str, str]]:
    candidates = []
    for row in manifest_rows:
        if row.get("prep_status") != "ok" or row.get("docking_status") != "ok":
            continue
        try:
            affinity = float(row["best_affinity_kcal_per_mol"])
            profile = property_profile(row["prepared_smiles"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"Cannot rank {row.get('chembl_id', '<unknown>')}: {exc}") from exc
        candidates.append((row, affinity, profile, review_flags(profile)))
    if not candidates:
        raise ValueError("No successfully docked manifest rows are available for ranking")
    candidates.sort(key=lambda item: (item[1], item[0]["chembl_id"]))
    total = len(candidates)
    ranked = []
    for rank, (row, affinity, profile, flags) in enumerate(candidates, start=1):
        shortlisted = rank <= shortlist_size
        ranked.append({
            "rank": str(rank),
            "chembl_id": row["chembl_id"],
            "pref_name": row.get("pref_name", ""),
            "best_affinity_kcal_per_mol": f"{affinity:.3f}",
            "docking_percentile": f"{100 * (total - rank + 1) / total:.2f}",
            "molecular_weight": f"{profile['molecular_weight']:.2f}",
            "rdkit_logp": f"{profile['rdkit_logp']:.2f}",
            "tpsa": f"{profile['tpsa']:.2f}",
            "hbd": str(profile["hbd"]), "hba": str(profile["hba"]),
            "rotatable_bonds": str(profile["rotatable_bonds"]), "formal_charge": str(profile["formal_charge"]),
            "review_flags": ";".join(flags),
            "second_pass_recommendation": (
                "rerun_high_exhaustiveness" if shortlisted and not flags else
                "rerun_high_exhaustiveness_with_flags" if shortlisted else "not_shortlisted"
            ),
            "canonical_smiles": row.get("canonical_smiles", ""),
            "prepared_smiles": row.get("prepared_smiles", ""),
            "docked_pdbqt_path": row.get("docked_pdbqt_path", ""),
        })
    return ranked


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=RANKING_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def rank_library(library_dir: Path = DEFAULT_LIBRARY_DIR, shortlist_size: int = 50) -> dict[str, Path | int]:
    if shortlist_size < 1:
        raise ValueError("Shortlist size must be at least 1")
    manifest_path = library_dir / "ligand_library_manifest.csv"
    with manifest_path.open(newline="", encoding="utf-8") as handle:
        ranked = ranked_rows(list(csv.DictReader(handle)), shortlist_size)
    ranking_dir = library_dir / "ranking"
    ranking_dir.mkdir(parents=True, exist_ok=True)
    ranked_path = ranking_dir / "katg_first_pass_ranked.csv"
    shortlist_path = ranking_dir / "katg_second_pass_shortlist.csv"
    summary_path = ranking_dir / "ranking_summary.md"
    write_csv(ranked_path, ranked)
    shortlist = [row for row in ranked if int(row["rank"]) <= shortlist_size]
    write_csv(shortlist_path, shortlist)
    affinities = [float(row["best_affinity_kcal_per_mol"]) for row in ranked]
    flagged = sum(bool(row["review_flags"]) for row in ranked)
    summary_path.write_text(
        "# KatG first-pass ranking summary\n\n"
        f"- Successfully docked ligands ranked: `{len(ranked)}`\n"
        f"- Second-pass shortlist size: `{len(shortlist)}`\n"
        f"- Rows with one or more review flags: `{flagged}`\n"
        f"- Best first-pass affinity: `{min(affinities):.3f}` kcal/mol\n"
        f"- Median first-pass affinity: `{statistics.median(affinities):.3f}` kcal/mol\n"
        f"- Least-favorable first-pass affinity: `{max(affinities):.3f}` kcal/mol\n\n"
        "## Interpretation\n\n"
        "Rows are ordered only by first-pass Vina affinity (more negative ranks higher). "
        "Property flags are review context, not exclusions or a validated efficacy score. "
        "Shortlisted rows should be rerun at higher exhaustiveness before further interpretation.\n",
        encoding="utf-8",
    )
    return {"ranked": ranked_path, "shortlist": shortlist_path, "summary": summary_path, "count": len(ranked)}


def main() -> None:
    parser = argparse.ArgumentParser(description="Rank completed KatG first-pass docking results")
    parser.add_argument("--library-dir", type=Path, default=DEFAULT_LIBRARY_DIR)
    parser.add_argument("--shortlist-size", type=int, default=50)
    args = parser.parse_args()
    print(rank_library(args.library_dir, args.shortlist_size))


if __name__ == "__main__":
    main()
