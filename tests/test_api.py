from fastapi.testclient import TestClient

from backend.api import app


client = TestClient(app)


def test_health_endpoint():
    response = client.get("/api/health")

    assert response.status_code == 200
    assert response.json() == {
        "status": "ok",
        "service": "genome-embeddings",
    }


def test_metadata_endpoint_describes_supported_uploads():
    response = client.get("/api/metadata")

    assert response.status_code == 200
    body = response.json()
    assert body["name"] == "Genome Embeddings"
    assert ".fasta" in body["upload_extensions"]
    assert body["supported_k_max"] == 8


def test_demo_endpoint_returns_records_and_analysis():
    response = client.get("/api/demo")

    assert response.status_code == 200
    body = response.json()
    assert len(body["records"]) == 6
    assert body["records"][0]["label"] == "Aequorea GFP"
    assert body["analysis"]["summary"]["dataset_size"] == 6
    assert "jensen_shannon" in body["analysis"]["matrices"]


def test_parse_record_endpoint_validates_fasta():
    response = client.post(
        "/api/records/parse",
        files={
            "file": (
                "example_sequence.fasta",
                b">Example\nACGTACGT\n",
                "text/plain",
            )
        },
    )

    assert response.status_code == 200
    assert response.json()["label"] == "example sequence"
    assert response.json()["sequence"] == "ACGTACGT"


def test_analyze_endpoint_is_stateless():
    request = {
        "records": [
            {
                "label": "Reference",
                "sequence": "ACGTACGT",
                "header": ">Reference",
                "source": "reference.fasta",
            },
            {
                "label": "Comparison",
                "sequence": "ACGTACGA",
                "header": ">Comparison",
                "source": "comparison.fasta",
            },
        ],
        "config": {
            "k_values": [1, 2, 3],
            "selected_k": 2,
            "reference_label": "Reference",
            "comparison_label": "Comparison",
        },
    }

    response = client.post("/api/analyze", json=request)

    assert response.status_code == 200
    body = response.json()
    assert body["summary"]["dataset_size"] == 2
    assert body["summary"]["selected_k"] == 2
    assert list(body["pair_trajectories"]["jensen_shannon"]) == [
        "1",
        "2",
        "3",
    ]


def test_export_endpoints_set_download_headers():
    request = {
        "records": [
            {
                "label": "Reference",
                "sequence": "ACGTACGT",
                "header": ">Reference",
                "source": "reference.fasta",
            },
            {
                "label": "Comparison",
                "sequence": "AAAACCCC",
                "header": ">Comparison",
                "source": "comparison.fasta",
            },
        ],
        "config": {
            "k_values": [1, 2],
            "selected_k": 2,
            "reference_label": "Reference",
            "comparison_label": "Comparison",
        },
    }

    json_response = client.post("/api/export/json", json=request)
    csv_response = client.post("/api/export/csv", json=request)

    assert json_response.status_code == 200
    assert "genome_embeddings_analysis.json" in json_response.headers[
        "content-disposition"
    ]
    assert csv_response.status_code == 200
    assert "genome_embeddings_summary.csv" in csv_response.headers[
        "content-disposition"
    ]


def test_analyze_endpoint_requires_two_kmer_scales():
    request = {
        "records": [
            {
                "label": "Reference",
                "sequence": "ACGTACGT",
                "header": ">Reference",
                "source": "reference.fasta",
            },
            {
                "label": "Comparison",
                "sequence": "ACGTACGA",
                "header": ">Comparison",
                "source": "comparison.fasta",
            },
        ],
        "config": {
            "k_values": [2],
            "selected_k": 2,
            "reference_label": "Reference",
            "comparison_label": "Comparison",
        },
    }

    response = client.post("/api/analyze", json=request)

    assert response.status_code == 422
