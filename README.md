# Celia AMR Prototype

First prototype for Celia Drug Discovery AI's antimicrobial resistance (AMR)
prediction pipeline, scoped to a single vertical slice: **Tuberculosis (TB)**.

Real computation on real data, no UI yet:

- **Genomic Analysis** — mutation lookups against a curated WHO TB resistance
  catalogue (`genomic_analysis/`)
- **Spread Modelling** — SEIR model parameterized with real WHO TB burden data
  (`spread_modelling/`)
- **Molecular Docking** — AutoDock Vina docking runs (`molecular_docking/`)
- **Risk Fusion Engine** — combines the three branches into a resistance
  score, mutation report, and alternative drug suggestions (`risk_fusion/`)

## Build status

Each branch is built and tested as an independent standalone script first,
sequenced cheapest-to-debug first. They are wired into a shared CLI only once
each branch is stable on its own.

| Branch            | Status      |
|--------------------|-------------|
| Genomic Analysis   | Done (katG + rpoB) |
| Spread Modelling   | Done (Zambia)      |
| Molecular Docking  | katG/Isoniazid verified end-to-end; rpoB/Rifampicin wired but not yet run (see `molecular_docking/README.md`) |
| Risk Fusion Engine | Done — fuses genomic + (optional) docking evidence into a resistance score, mutation report, and alternative-drug suggestions; WHO burden data is surfaced as separate public-health context, not blended into the score (see `risk_fusion/fuse_scores.py` docstring) |

All four branches are independently runnable and tested; they are not yet
wired into one shared `celia analyze` CLI entrypoint (still a TODO per the
original build order — see the brief).

## Environment setup

```bash
conda env create -f environment.yml
conda activate celia-amr
```

## Running tests

```bash
pytest
```

## First drug/target pairs

Both are covered by the WHO TB mutation catalogue and built side by side:

- Isoniazid vs. *katG*
- Rifampicin vs. *rpoB*

## Country for spread modelling

Zambia — see `data/tb_burden/README.md` for the WHO data sources and
`spread_modelling/seir_model.py`'s module docstring for how each SEIR
parameter is derived from real WHO figures vs. literature assumptions.

## Data sources

- **Genomic Analysis**: NCBI (genome data) + WHO TB mutation catalogue
  (resistance-conferring mutations)
- **Molecular Docking**: RCSB PDB (protein structures) + PubChem (compound
  structures)
- **Spread Modelling**: WHO Global TB burden estimates and treatment outcome
  data for Zambia (`https://extranet.who.int/tme/generateCSV.asp?ds=estimates`,
  `...?ds=outcomes`)
