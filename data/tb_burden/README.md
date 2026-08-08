# WHO TB Burden Data

Raw data downloaded directly from WHO's Global Tuberculosis Programme data
endpoints (2026-08-08):

- `TB_burden_estimates_full.csv` — `https://extranet.who.int/tme/generateCSV.asp?ds=estimates`
  Country-year estimates of incidence, mortality, and case fatality ratio,
  with uncertainty bounds. This is the same dataset behind the WHO Global TB
  Report country profiles.
- `TB_outcomes_full.csv` — `https://extranet.who.int/tme/generateCSV.asp?ds=outcomes`
  Country-year treatment cohort outcomes (treatment success rate, failure,
  death, loss to follow-up), by cohort year.

Both are the full multi-country files as published (~5,300-6,400 rows). For
convenience, `zambia_tb_burden.csv` and `zambia_tb_outcomes.csv` are the
Zambia-only rows extracted from each, which is what
`spread_modelling/seir_model.py` reads.

## Columns used by the SEIR model

From `zambia_tb_burden.csv`:
- `e_pop_num` — population estimate
- `e_inc_num` — estimated incident TB cases (all forms) in the year
- `cfr` — case fatality ratio (proportion of incident cases that die)

From `zambia_tb_outcomes.csv`:
- `c_new_tsr` — treatment success rate (%) among new/relapse cases, cited in
  the model's documentation for context but not a direct input (the fitted
  `cfr` already captures the population-level death outcome across both
  treated and untreated cases).

## Refreshing this data

Re-run the download and re-extract Zambia's rows if newer WHO report data is
published:

```bash
curl -sL "https://extranet.who.int/tme/generateCSV.asp?ds=estimates" -o TB_burden_estimates_full.csv
curl -sL "https://extranet.who.int/tme/generateCSV.asp?ds=outcomes" -o TB_outcomes_full.csv
```
