# Celia browser frontend

This browser-facing slice has two deliberately separate modes:

- **Static demo (default):** used by GitHub Pages. It runs entirely in the
  browser with bundled illustrative data, so no server is exposed publicly.
- **Live local mode:** calls the existing local analysis pipeline through the
  Python server, for technical verification.

Run it from the repository root:

```bash
python frontend/server.py
```

Then open `http://127.0.0.1:8000/?live=1` to use the local backend. Opening
the URL without `?live=1` runs the static demo instead.

The repository's **Backend verification** GitHub Actions workflow runs the
Python test suite on pushes and pull requests. This gives reviewers an
independent, repeatable check of the backend even though GitHub Pages serves
only the static demonstration.
