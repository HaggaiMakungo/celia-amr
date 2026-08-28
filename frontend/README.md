# Celia browser frontend

This is the first browser-facing slice of Celia AMR. It keeps the approved
workspace design while calling the existing local genomic-analysis and
risk-fusion code for bundled `katG` and `rpoB` examples.

Run it from the repository root:

```bash
python frontend/server.py
```

Then open `http://127.0.0.1:8000`.
