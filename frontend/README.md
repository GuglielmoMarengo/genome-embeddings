# Genome Embeddings Frontend

Professional React and TypeScript client for the Genome Embeddings scientific API.

## Stack

- React
- TypeScript
- Vite
- Tailwind CSS
- Motion
- Plotly.js
- Lucide icons

## Development

From the repository root, install the Python dependencies first:

```powershell
python -m pip install -r requirements.txt
```

Then install the frontend dependencies:

```powershell
cd frontend
npm install
```

Run both development servers from the repository root:

```powershell
python dev.py
```

The frontend is available at `http://127.0.0.1:5173` and proxies `/api` to the FastAPI backend at `http://127.0.0.1:8000`.

## Production build

```powershell
cd frontend
npm run build
cd ..
python app.py
```

FastAPI serves the generated `frontend/dist` application and the API from the same origin.
