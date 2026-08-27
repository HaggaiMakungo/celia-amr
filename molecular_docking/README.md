# Molecular Docking

Real AutoDock Vina docking runs for the TB drug/target pairs used elsewhere
in this project. See `run_docking.py`'s module docstring for the full
rationale (why no conda/OpenBabel/PyPI `vina` here, what each target's
docking box is centered on, and the katG covalent-crosslink handling).

## Status

| Pair                        | Status                                                    |
|------------------------------|------------------------------------------------------------|
| katG / Isoniazid (PDB 2CCA)  | Verified end-to-end. Best pose ~ -5.8 kcal/mol.            |
| rpoB / Rifampicin (PDB 5UH6) | Wired into the target registry, **not yet run**. 5UH6 is a large multi-subunit cryo-EM complex; expect a longer/likely-rockier receptor-prep pass than katG's. Run it and fix whatever meeko/Vina complains about before trusting the result. |

## Setup (this machine: no conda, Python 3.14 default)

RDKit/meeko/etc. don't ship wheels for Python 3.14 yet, and the PyPI `vina`
package has no Windows wheels at all (its sdist needs Boost). Rather than
require installing conda, this project uses a plain Python 3.11 virtualenv
plus a directly-downloaded Vina binary:

```bash
py -3.11 -m venv .venv311
.venv311/Scripts/python.exe -m pip install rdkit meeko gemmi truststore pubchempy pandas numpy scipy pytest biopython requests
python molecular_docking/run_docking.py --gene katG
```

`run_docking.py` auto-downloads the real PDB structure, the real PubChem
compound, and (on first run) the official `vina.exe` binary into
`molecular_docking/bin/` — none of these are committed to the repo (see
`.gitignore`).

If your machine's Python has a working OpenSSL trust store and/or conda
available, you likely don't need any of this — `environment.yml`'s conda
env should work directly.

## Running

```bash
python molecular_docking/run_docking.py --gene katG
python molecular_docking/run_docking.py --gene rpoB   # not yet verified, see Status above
```

Intermediate/output files (cleaned receptor, prepared PDBQTs, docked poses)
are written to `molecular_docking/work/` (not committed — regenerate by
re-running).

## Approved-drug ligand library

Build the ChEMBL phase-4, small-molecule library before adding batch docking:

```bash
python molecular_docking/build_ligand_library.py
```

It writes `data/ligand_library/chembl_approved_v1/`, containing one PDBQT per
successfully prepared ChEMBL molecule, `ligand_library_manifest.csv`, and
`SOURCE.md`. The generator calls the exact RDKit/Meeko functions used by the
verified isoniazid run: deterministic ETKDG embedding (seed 42), MMFF
optimisation, then `mk_prepare_ligand`. The validated baseline has no explicit
pH adjustment, so this library intentionally does not introduce one.

Before preparation, the builder rejects carbon-free structures as
`rejected_inorganic`. It applies RDKit salt removal to every entry and, when
an organic counterion remains, uses the ChEMBL parent record already in the
pull. Entries with no usable single-fragment parent are retained in the
manifest as `rejected_unresolved_multifragment` for manual review.

The command stops before structure preparation if the post-cleanup ready set
contains fewer than 1,500 or more than 2,500 ligands. Inspect the cleanup
result before overriding this check with `--allow-unexpected-record-count`.

## Batch docking

After ligand preparation completes, begin with the resumable 24-ligand KatG
test batch:

```bash
python molecular_docking/batch_docking.py
```

This uses the prepared-library manifest as its only ligand input, runs Vina at
first-pass exhaustiveness 4, and uses one CPU per Vina process while running
six molecules concurrently (the six physical cores on this machine). It adds docking status/score columns to the
manifest, writes one result PDBQT per ChEMBL ID under `docking_results/`, and
persists progress after every ligand. Individual Vina runs time out after 15
minutes by default and are recorded as failures; override with
`--ligand-timeout-minutes`. Run `--full` only after reviewing the test batch's
`docking_summary.md`.

## First-pass ranking

After docking completes, create the first-pass KatG ranking and a review
shortlist:

```bash
python molecular_docking/rank_ligands.py
```

This writes ranked and shortlist CSV files plus `ranking_summary.md` under
`data/ligand_library/chembl_approved_v1/ranking/`. Candidates are ordered only
by first-pass Vina affinity; RDKit property flags provide review context rather
than an unvalidated efficacy score. The top 50 are marked for a later
high-exhaustiveness rerun.

## Second-pass shortlist rerun

Rerun the ranked shortlist with the validated baseline exhaustiveness of 8:

```bash
python molecular_docking/rerun_shortlist.py
```

Second-pass output is kept separate from first-pass docking, including its own
result PDBQTs, failure log, and summary. It is resumable and uses a 30-minute
per-ligand timeout because these are deliberately more thorough searches.

## Final review ranking

Fold confirmed second-pass scores into an auditable final review table:

```bash
python molecular_docking/finalize_ranking.py
```

The final table replaces first-pass scores only for successfully rerun shortlist
ligands and labels every score source explicitly.
