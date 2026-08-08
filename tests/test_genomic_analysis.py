import pytest

from genomic_analysis.analyze_mutations import (
    CATALOGUE_PATH,
    call_mutations,
    load_catalogue,
    parse_mutation,
    read_fasta_sequence,
    summarize,
)


def make_sequence(length: int, overrides: dict[int, str]) -> str:
    """Build a filler protein sequence with specific 1-indexed positions set."""
    chars = ["A"] * length
    for position, aa in overrides.items():
        chars[position - 1] = aa
    return "".join(chars)


def test_parse_mutation():
    assert parse_mutation("S315T") == ("S", 315, "T")


def test_parse_mutation_rejects_malformed_input():
    with pytest.raises(ValueError):
        parse_mutation("not-a-mutation")


def test_load_catalogue_covers_both_genes():
    catalogue = load_catalogue()
    genes = {entry.gene for entry in catalogue}
    assert genes == {"katG", "rpoB"}


def test_call_mutations_detects_resistant_katg_s315t():
    catalogue = load_catalogue()
    sequence = make_sequence(420, {315: "T"})  # katG S315T
    calls = call_mutations("katG", sequence, catalogue)
    hit = next(c for c in calls if c.entry.mutation == "S315T")
    assert hit.status == "resistant"


def test_call_mutations_detects_wild_type():
    catalogue = load_catalogue()
    sequence = make_sequence(420, {315: "S"})  # unmutated katG
    calls = call_mutations("katG", sequence, catalogue)
    hit = next(c for c in calls if c.entry.mutation == "S315T")
    assert hit.status == "wild_type"


def test_call_mutations_flags_novel_variant():
    catalogue = load_catalogue()
    sequence = make_sequence(420, {315: "Q"})  # neither wild-type nor cataloged mutant
    calls = call_mutations("katG", sequence, catalogue)
    hit = next(c for c in calls if c.entry.mutation == "S315T")
    assert hit.status == "novel_variant"


def test_call_mutations_rejects_short_sequence():
    catalogue = load_catalogue()
    sequence = make_sequence(10, {})
    with pytest.raises(ValueError):
        call_mutations("katG", sequence, catalogue)


def test_summarize_predicts_resistant_phenotype():
    catalogue = load_catalogue()
    sequence = make_sequence(420, {315: "T"})
    calls = call_mutations("katG", sequence, catalogue)
    summary = summarize(calls)
    assert summary["predicted_phenotype"] == "resistant"
    assert summary["driving_mutation"] == "S315T"


def test_summarize_predicts_susceptible_phenotype():
    catalogue = load_catalogue()
    sequence = make_sequence(420, {315: "S", 328: "W", 419: "D"})
    calls = call_mutations("katG", sequence, catalogue)
    summary = summarize(calls)
    assert summary["predicted_phenotype"] == "susceptible"
    assert summary["driving_mutation"] is None


def test_summarize_ignores_non_resistance_associated_grade_matches():
    # D419A is WHO grade 4 ("Not assoc w R - Interim"): matching the mutant
    # residue there must not drive a resistant phenotype call by itself.
    catalogue = load_catalogue()
    sequence = make_sequence(420, {315: "S", 328: "W", 419: "A"})
    calls = call_mutations("katG", sequence, catalogue)
    summary = summarize(calls)
    assert summary["predicted_phenotype"] == "susceptible"
    assert summary["driving_mutation"] is None
    assert "D419A" not in summary["resistant_mutations_found"]


def test_rpob_s450l_detected_as_resistant():
    catalogue = load_catalogue()
    sequence = make_sequence(500, {450: "L"})  # rpoB S450L
    calls = call_mutations("rpoB", sequence, catalogue)
    hit = next(c for c in calls if c.entry.mutation == "S450L")
    assert hit.status == "resistant"


def test_read_fasta_sequence(tmp_path):
    fasta = tmp_path / "sample.fasta"
    fasta.write_text(">katG sample\nMTAG\nCDEF\n")
    assert read_fasta_sequence(fasta) == "MTAGCDEF"


def test_read_fasta_sequence_rejects_invalid_characters(tmp_path):
    fasta = tmp_path / "bad.fasta"
    fasta.write_text(">bad\nMTAGZZZ\n")
    with pytest.raises(ValueError):
        read_fasta_sequence(fasta)


def test_catalogue_file_exists_on_disk():
    assert CATALOGUE_PATH.exists()
