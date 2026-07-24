# Frontend and UX Release

This release replaces the primary NiceGUI interface with a separated FastAPI and React architecture.

## What changed

- FastAPI stateless API under `backend/`
- React and TypeScript application under `frontend/`
- Vite development server and production build
- Tailwind CSS design tokens and responsive layouts
- Motion page transitions, panels, loading states and micro-interactions
- Plotly.js heatmaps and charts with responsive resize, zoom and image export
- Persistent light and dark themes applied before React starts
- Collapsible desktop navigation and mobile drawer
- Guided off-canvas analysis setup
- Drag-and-drop FASTA upload with backend validation
- Searchable dataset inventory and sequence previews
- Descriptor V2 entropy, coverage, effective-count and dinucleotide charts
- Fullscreen matrix and multiscale visualizations
- JSON and CSV exports through stateless API endpoints

## Theme behavior

The theme is stored under:

```text
genome-embeddings-theme
```

The initial theme is applied by an inline script in `frontend/index.html` before the application bundle loads. This avoids the previous non-functional toggle and reduces light/dark flashing during startup.

## Development

```powershell
python -m pip install -r requirements.txt
cd frontend
npm install
cd ..
python dev.py
```

Open `http://127.0.0.1:5173`.

## Production

```powershell
cd frontend
npm run typecheck
npm run build
cd ..
python app.py
```

Open `http://127.0.0.1:8000`.

## Legacy interface

The previous NiceGUI application remains available as `legacy_app.py`. Install its optional dependencies with:

```powershell
python -m pip install -r requirements-legacy-ui.txt
python legacy_app.py
```
