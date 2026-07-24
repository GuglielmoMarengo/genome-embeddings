import json
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
FRONTEND_ROOT = PROJECT_ROOT / "frontend"


def test_frontend_uses_react_vite_tailwind_and_motion():
    package = json.loads((FRONTEND_ROOT / "package.json").read_text())

    assert package["dependencies"]["react"].startswith("19.")
    assert "motion" in package["dependencies"]
    assert "plotly.js-dist-min" in package["dependencies"]
    assert "@tailwindcss/vite" in package["devDependencies"]
    assert "vite" in package["devDependencies"]


def test_theme_provider_persists_and_applies_dark_class():
    source = (FRONTEND_ROOT / "src" / "components" / "ThemeProvider.tsx").read_text()

    assert "genome-embeddings-theme" in source
    assert "classList.toggle('dark'" in source
    assert "prefers-color-scheme: dark" in source


def test_vite_proxies_api_to_fastapi_backend():
    source = (FRONTEND_ROOT / "vite.config.ts").read_text()

    assert "'/api': 'http://127.0.0.1:8000'" in source


def test_frontend_has_professional_navigation_pages():
    source = (FRONTEND_ROOT / "src" / "components" / "AppShell.tsx").read_text()

    for label in (
        "Overview",
        "Dataset",
        "Descriptors V2",
        "Matrices",
        "Multiscale",
        "Exports",
        "Methodology",
    ):
        assert label in source


def test_frontend_applies_theme_before_react_hydration():
    source = (FRONTEND_ROOT / "index.html").read_text()

    assert "genome-embeddings-theme" in source
    assert "document.documentElement.classList.toggle('dark'" in source
