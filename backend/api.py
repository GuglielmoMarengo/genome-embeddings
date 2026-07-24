"""FastAPI application serving the scientific API and production frontend."""

from __future__ import annotations

from pathlib import Path

from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse, Response
from fastapi.staticfiles import StaticFiles

from backend.schemas import AnalysisRequest, ApiMetadata, DemoResponse, SequencePayload
from backend.service import demo_request, parse_upload, run_request


PROJECT_ROOT = Path(__file__).resolve().parents[1]
FRONTEND_DIST = PROJECT_ROOT / "frontend" / "dist"

app = FastAPI(
    title="Genome Embeddings API",
    description=(
        "Stateless API for interpretable genomic descriptors, multiscale "
        "embeddings, and k-mer distribution comparison."
    ),
    version="0.4.0",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "http://127.0.0.1:5173",
        "http://localhost:5173",
        "http://127.0.0.1:8000",
        "http://localhost:8000",
    ],
    allow_credentials=False,
    allow_methods=["GET", "POST"],
    allow_headers=["*"],
)


@app.get("/api/health")
def health() -> dict[str, str]:
    return {"status": "ok", "service": "genome-embeddings"}


@app.get("/api/metadata", response_model=ApiMetadata)
def metadata() -> ApiMetadata:
    return ApiMetadata(
        name="Genome Embeddings",
        tagline="Turning genomes into mathematics.",
        version=app.version,
        supported_k_min=1,
        supported_k_max=8,
        upload_extensions=[".fasta", ".fa", ".fna", ".txt"],
    )


@app.get("/api/demo", response_model=DemoResponse)
def demo() -> DemoResponse:
    request = demo_request()
    analysis = run_request(request)
    return DemoResponse(records=request.records, analysis=analysis.to_dict())


@app.post("/api/records/parse", response_model=SequencePayload)
async def parse_record(
    file: UploadFile = File(...),
    label: str | None = Form(default=None),
) -> SequencePayload:
    try:
        data = await file.read()
        return parse_upload(file.filename or "uploaded.fasta", data, label)
    except (TypeError, ValueError) as error:
        raise HTTPException(status_code=400, detail=str(error)) from error


@app.post("/api/analyze")
def analyze(request: AnalysisRequest) -> dict[str, object]:
    try:
        return run_request(request).to_dict()
    except (TypeError, ValueError) as error:
        raise HTTPException(status_code=400, detail=str(error)) from error


@app.post("/api/export/json")
def export_json(request: AnalysisRequest) -> Response:
    try:
        content = run_request(request).to_json(indent=2)
    except (TypeError, ValueError) as error:
        raise HTTPException(status_code=400, detail=str(error)) from error
    return Response(
        content=content,
        media_type="application/json",
        headers={
            "Content-Disposition": (
                'attachment; filename="genome_embeddings_analysis.json"'
            )
        },
    )


@app.post("/api/export/csv")
def export_csv(request: AnalysisRequest) -> Response:
    try:
        content = run_request(request).summary_csv()
    except (TypeError, ValueError) as error:
        raise HTTPException(status_code=400, detail=str(error)) from error
    return Response(
        content="\ufeff" + content,
        media_type="text/csv; charset=utf-8",
        headers={
            "Content-Disposition": (
                'attachment; filename="genome_embeddings_summary.csv"'
            )
        },
    )


@app.exception_handler(Exception)
async def unhandled_exception_handler(_, error: Exception) -> JSONResponse:
    return JSONResponse(
        status_code=500,
        content={"detail": "Unexpected server error.", "type": type(error).__name__},
    )


if FRONTEND_DIST.is_dir():
    assets_dir = FRONTEND_DIST / "assets"
    if assets_dir.is_dir():
        app.mount("/assets", StaticFiles(directory=assets_dir), name="assets")

    @app.get("/{path:path}", include_in_schema=False)
    def frontend(path: str) -> FileResponse:
        candidate = FRONTEND_DIST / path
        if path and candidate.is_file():
            return FileResponse(candidate)
        return FileResponse(FRONTEND_DIST / "index.html")
else:
    @app.get("/", include_in_schema=False)
    def frontend_not_built() -> JSONResponse:
        return JSONResponse(
            {
                "message": "Frontend build not found.",
                "development": "cd frontend && npm install && npm run dev",
                "production": "cd frontend && npm install && npm run build",
                "api_docs": "/docs",
            }
        )
