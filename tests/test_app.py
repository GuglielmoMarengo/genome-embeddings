import app

from src.genome import GenomeMatrix


def test_app_module_imports_without_starting_server():
    assert callable(app.create_dashboard)
    assert callable(app.main)


def test_matrix_figure_uses_matrix_labels_and_metric_title():
    matrix = GenomeMatrix(
        labels=["A", "B"],
        values=[[0.0, 0.3], [0.3, 0.0]],
        metric="jensen_shannon",
        kmer_length=2,
    )

    figure = app.matrix_figure(matrix)

    assert figure.layout.title.text == (
        f"{app.metric_title('jensen_shannon')} · k=2"
    )
    assert tuple(figure.data[0].x) == ("A", "B")


def test_unique_label_adds_incrementing_suffix():
    from src.dashboard import DatasetRecord
    from src.genome import Genome

    records = [
        DatasetRecord("Sample", Genome("ACGT"), "one.fasta"),
        DatasetRecord("Sample (2)", Genome("ACGT"), "two.fasta"),
    ]

    assert app.unique_label(records, "Sample") == "Sample (3)"
