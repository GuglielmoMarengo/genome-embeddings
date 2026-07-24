"""Launch the Genome Embeddings web application.

After a production frontend build, ``python app.py`` serves FastAPI and the
React application from http://127.0.0.1:8000.

When ``frontend/dist`` is absent but frontend dependencies are installed, the
same command falls back to development mode and starts Vite with hot reload.
"""

from __future__ import annotations

import sys
from pathlib import Path

import uvicorn


PROJECT_ROOT = Path(__file__).resolve().parent
FRONTEND_ROOT = PROJECT_ROOT / "frontend"
FRONTEND_DIST = FRONTEND_ROOT / "dist" / "index.html"
NODE_MODULES = FRONTEND_ROOT / "node_modules"


def main() -> None:
    if FRONTEND_DIST.is_file():
        uvicorn.run(
            "backend.api:app",
            host="127.0.0.1",
            port=8000,
            reload=False,
            log_level="info",
        )
        return

    if NODE_MODULES.is_dir():
        print(
            "Frontend production build not found; starting development mode "
            "with Vite at http://127.0.0.1:5173."
        )
        from dev import main as development_main

        development_main()
        return

    raise SystemExit(
        "Frontend dependencies are not installed.\n\n"
        "Run:\n"
        "  cd frontend\n"
        "  npm install\n"
        "  cd ..\n"
        "  python app.py\n\n"
        "For a production build, run `npm run build` inside frontend first."
    )


if __name__ == "__main__":
    main()
