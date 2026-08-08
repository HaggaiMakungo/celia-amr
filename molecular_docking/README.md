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
