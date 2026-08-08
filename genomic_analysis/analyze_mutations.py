"""Branch 1: Genomic Analysis.

Looks up amino-acid substitutions in a query protein sequence against a
curated WHO TB mutation catalogue, for a given gene/drug pair (katG/Isoniazid
or rpoB/Rifampicin).

Standalone script — no dependency on the other three branches. Run directly:

    python analyze_mutations.py --gene katG --fasta path/to/katG_query.fasta
"""

import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path

CATALOGUE_PATH = (
    Path(__file__).resolve().parent.parent
    / "data"
    / "who_tb_mutation_catalogue"
    / "tb_resistance_mutations.csv"
)

MUTATION_PATTERN = re.compile(r"^([A-Z])(\d+)([A-Z])$")

# Standard single-letter amino acid codes (protein FASTA alphabet).
VALID_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


@dataclass(frozen=True)
class CatalogueEntry:
    gene: str
    drug: str
    mutation: str
    wild_type_aa: str
    position: int
    mutant_aa: str
    who_confidence_grade: int
    confidence_label: str
    notes: str


def parse_mutation(mutation: str) -> tuple[str, int, str]:
    """Parse 'S315T' into ('S', 315, 'T')."""
    match = MUTATION_PATTERN.match(mutation.strip())
    if not match:
        raise ValueError(f"Malformed mutation notation: {mutation!r}")
    wild_type_aa, position, mutant_aa = match.groups()
    return wild_type_aa, int(position), mutant_aa


def load_catalogue(csv_path: Path = CATALOGUE_PATH) -> list[CatalogueEntry]:
    """Load the curated WHO TB mutation catalogue from disk."""
    entries = []
    with open(csv_path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            wild_type_aa, position, mutant_aa = parse_mutation(row["mutation"])
            entries.append(
                CatalogueEntry(
                    gene=row["gene"],
                    drug=row["drug"],
                    mutation=row["mutation"],
                    wild_type_aa=wild_type_aa,
                    position=position,
                    mutant_aa=mutant_aa,
                    who_confidence_grade=int(row["who_confidence_grade"]),
                    confidence_label=row["confidence_label"],
                    notes=row["notes"],
                )
            )
    return entries


def read_fasta_sequence(fasta_path: Path) -> str:
    """Read a single-record protein FASTA file and return the sequence."""
    lines = Path(fasta_path).read_text(encoding="utf-8-sig").splitlines()
    seq_lines = [line.strip() for line in lines if line and not line.startswith(">")]
    sequence = "".join(seq_lines).upper()
    invalid = set(sequence) - VALID_AMINO_ACIDS
    if invalid:
        raise ValueError(f"Sequence contains non-amino-acid characters: {sorted(invalid)}")
    return sequence


@dataclass(frozen=True)
class MutationCall:
    entry: CatalogueEntry
    observed_aa: str
    status: str  # "resistant", "wild_type", "novel_variant"


def call_mutations(gene: str, sequence: str, catalogue: list[CatalogueEntry]) -> list[MutationCall]:
    """Compare a query sequence against each cataloged position for `gene`.

    For each cataloged resistance site, the residue observed in `sequence` at
    that 1-indexed position is classified as:
      - "resistant"     if it matches the cataloged mutant residue
      - "wild_type"      if it matches the cataloged wild-type residue
      - "novel_variant"  if it's neither (a residue not in the catalogue)
    """
    calls = []
    for entry in catalogue:
        if entry.gene != gene:
            continue
        if entry.position > len(sequence):
            raise ValueError(
                f"Sequence too short ({len(sequence)} aa) to cover cataloged "
                f"position {entry.position} ({entry.mutation})"
            )
        observed_aa = sequence[entry.position - 1]
        if observed_aa == entry.mutant_aa:
            status = "resistant"
        elif observed_aa == entry.wild_type_aa:
            status = "wild_type"
        else:
            status = "novel_variant"
        calls.append(MutationCall(entry=entry, observed_aa=observed_aa, status=status))
    return calls


RESISTANCE_ASSOCIATED_GRADES = {1, 2}  # WHO grades "Assoc w R" / "Assoc w R - Interim"


def summarize(calls: list[MutationCall]) -> dict:
    """Roll per-position calls up into a gene-level resistance summary.

    Only matches against WHO grade 1/2 (resistance-associated) catalogue
    entries drive the predicted phenotype — a position matching a grade 4/5
    ("not associated with resistance") variant is reported in the raw calls
    but must not flip the phenotype to resistant.
    """
    resistant_calls = [
        c
        for c in calls
        if c.status == "resistant" and c.entry.who_confidence_grade in RESISTANCE_ASSOCIATED_GRADES
    ]
    if resistant_calls:
        top = max(resistant_calls, key=lambda c: -c.entry.who_confidence_grade)
        predicted = "resistant"
    else:
        top = None
        predicted = "susceptible"
    return {
        "predicted_phenotype": predicted,
        "driving_mutation": top.entry.mutation if top else None,
        "who_confidence_grade": top.entry.who_confidence_grade if top else None,
        "confidence_label": top.entry.confidence_label if top else None,
        "resistant_mutations_found": [c.entry.mutation for c in resistant_calls],
        "positions_checked": len(calls),
    }


def analyze(gene: str, fasta_path: Path, catalogue_path: Path = CATALOGUE_PATH) -> dict:
    catalogue = load_catalogue(catalogue_path)
    if not any(entry.gene == gene for entry in catalogue):
        raise ValueError(f"No catalogue entries for gene {gene!r}")
    sequence = read_fasta_sequence(fasta_path)
    calls = call_mutations(gene, sequence, catalogue)
    summary = summarize(calls)
    summary["gene"] = gene
    summary["drug"] = catalogue[0].drug if catalogue[0].gene == gene else next(
        e.drug for e in catalogue if e.gene == gene
    )
    summary["calls"] = [
        {
            "mutation": c.entry.mutation,
            "observed_aa": c.observed_aa,
            "status": c.status,
            "who_confidence_grade": c.entry.who_confidence_grade,
            "confidence_label": c.entry.confidence_label,
        }
        for c in calls
    ]
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="TB resistance mutation lookup")
    parser.add_argument("--gene", required=True, choices=["katG", "rpoB"])
    parser.add_argument("--fasta", required=True, type=Path, help="Protein FASTA file for the gene")
    parser.add_argument("--catalogue", type=Path, default=CATALOGUE_PATH)
    args = parser.parse_args()

    result = analyze(args.gene, args.fasta, args.catalogue)

    print(f"Gene: {result['gene']}  Drug: {result['drug']}")
    print(f"Predicted phenotype: {result['predicted_phenotype']}")
    if result["driving_mutation"]:
        print(
            f"Driving mutation: {result['driving_mutation']} "
            f"(WHO grade {result['who_confidence_grade']}: {result['confidence_label']})"
        )
    print(f"Positions checked: {result['positions_checked']}")
    for call in result["calls"]:
        print(f"  {call['mutation']}: observed {call['observed_aa']} -> {call['status']}")


if __name__ == "__main__":
    main()
