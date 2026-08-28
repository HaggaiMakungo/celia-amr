"""Local browser frontend for the Celia AMR prototype.

The server deliberately uses the standard library so the first UI slice can run
without changing the scientific environment.  The analysis endpoint calls the
existing Celia pipeline rather than duplicating its logic in the browser.
"""

from __future__ import annotations

import argparse
import json
import sys
from http import HTTPStatus
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlparse


FRONTEND_ROOT = Path(__file__).resolve().parent
PROJECT_ROOT = FRONTEND_ROOT.parent
sys.path.insert(0, str(PROJECT_ROOT))

from celia import report_as_dict, run_analysis  # noqa: E402


SAMPLES = {
    "katG": PROJECT_ROOT / "genomic_analysis/examples/katG_resistant_example.fasta",
    "rpoB": PROJECT_ROOT / "genomic_analysis/examples/rpoB_resistant_example.fasta",
}


def project_payload() -> dict:
    """Return stable project context while the UI waits for an assessment."""
    return {
        "project": {
            "code": "PRJ-042",
            "name": "Drug-resistant TB investigation",
            "status": "Ready to explore",
            "question": "Which candidates warrant targeted validation against the resistant profile?",
        },
        "candidates": [
            {"name": "ZMB-041", "score": "−9.2", "fit": "High", "status": "Lead candidate"},
            {"name": "ZMB-117", "score": "−8.7", "fit": "Promising", "status": "Evidence review"},
            {"name": "ZMB-203", "score": "−8.4", "fit": "Promising", "status": "Evidence review"},
        ],
        "scenarios": {"targeted": "−12%", "baseline": "+18%", "preferred": "Targeted response"},
    }


def analysis_payload(gene: str) -> dict:
    """Run a bundled example through the real genomic and fusion workflow."""
    if gene not in SAMPLES:
        raise ValueError("Choose either katG or rpoB.")

    assessment, genomic_result, docking_result = run_analysis(
        gene,
        SAMPLES[gene],
        with_burden_context=True,
    )
    return report_as_dict(assessment, genomic_result, docking_result)


class CeliaHandler(SimpleHTTPRequestHandler):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=str(FRONTEND_ROOT), **kwargs)

    def do_GET(self) -> None:  # noqa: N802
        route = urlparse(self.path).path
        if route == "/api/health":
            self.send_json({"status": "ok"})
            return
        if route == "/api/project":
            self.send_json(project_payload())
            return
        if route == "/api/analysis/default":
            self.run_analysis("katG")
            return
        super().do_GET()

    def do_POST(self) -> None:  # noqa: N802
        route = urlparse(self.path).path
        if route != "/api/analysis":
            self.send_error(HTTPStatus.NOT_FOUND, "Unknown API route")
            return

        try:
            content_length = int(self.headers.get("Content-Length", "0"))
            if not 0 < content_length <= 2048:
                raise ValueError("Request body must be between 1 and 2048 bytes.")
            request = json.loads(self.rfile.read(content_length))
            self.run_analysis(request.get("gene", ""))
        except (json.JSONDecodeError, ValueError) as error:
            self.send_json({"error": str(error)}, status=HTTPStatus.BAD_REQUEST)
        except Exception:  # The detailed error belongs in the server log, not the UI.
            self.send_json({"error": "The assessment could not be completed."}, status=HTTPStatus.INTERNAL_SERVER_ERROR)

    def run_analysis(self, gene: str) -> None:
        try:
            self.send_json(analysis_payload(gene))
        except ValueError as error:
            self.send_json({"error": str(error)}, status=HTTPStatus.BAD_REQUEST)
        except Exception:
            self.send_json({"error": "The assessment could not be completed."}, status=HTTPStatus.INTERNAL_SERVER_ERROR)

    def send_json(self, payload: dict, status: HTTPStatus = HTTPStatus.OK) -> None:
        encoded = json.dumps(payload).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(encoded)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(encoded)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the local Celia browser frontend.")
    parser.add_argument("--port", type=int, default=8000)
    args = parser.parse_args()
    server = ThreadingHTTPServer(("127.0.0.1", args.port), CeliaHandler)
    print(f"Celia frontend running at http://127.0.0.1:{args.port}")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
