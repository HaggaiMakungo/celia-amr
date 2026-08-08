# WHO TB Mutation Catalogue (curated seed subset)

`tb_resistance_mutations.csv` is a small, hand-curated seed of well-established
katG/Isoniazid and rpoB/Rifampicin resistance mutations, structured for direct
lookup by `genomic_analysis/analyze_mutations.py`.

This is **not** the full WHO catalogue — the real catalogue is published as a
spreadsheet/PDF (not an API) and covers many more genes, drugs, and mutation
types (SNPs, indels, promoter variants) with per-mutation evidence grades.
Expanding this CSV by parsing the full published catalogue into this same
schema is open work (see the project brief's "open questions").

## Schema

| column                | meaning                                                              |
|------------------------|-----------------------------------------------------------------------|
| `gene`                 | Gene symbol (`katG`, `rpoB`)                                         |
| `drug`                 | Drug the mutation confers resistance to                              |
| `mutation`              | Amino-acid substitution in `<wild-type><position><mutant>` notation  |
| `who_confidence_grade` | WHO catalogue grade: 1=Assoc w R, 2=Assoc w R (Interim), 3=Uncertain, 4=Not assoc w R (Interim), 5=Not assoc w R |
| `confidence_label`      | Human-readable form of the grade                                     |
| `notes`                 | Free-text context                                                     |

Only grades 1-2 drive a "resistant" phenotype prediction in
`analyze_mutations.summarize()` — grade 3-5 matches are reported per-position
but do not by themselves flip the call.
