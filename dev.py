"""Run the FastAPI backend and Vite frontend development servers together."""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parent
FRONTEND_ROOT = PROJECT_ROOT / "frontend"


def main() -> None:
    npm = shutil.which("npm")
    if npm is None:
        raise SystemExit("npm was not found. Install Node.js 22 or newer.")

    backend = subprocess.Popen(
        [
            sys.executable,
            "-m",
            "uvicorn",
            "backend.api:app",
            "--host",
            "127.0.0.1",
            "--port",
            "8000",
            "--reload",
        ],
        cwd=PROJECT_ROOT,
    )

    try:
        subprocess.run(
            [npm, "run", "dev"],
            cwd=FRONTEND_ROOT,
            check=True,
        )
    except KeyboardInterrupt:
        pass
    finally:
        backend.terminate()
        try:
            backend.wait(timeout=5)
        except subprocess.TimeoutExpired:
            backend.kill()


if __name__ == "__main__":
    main()
