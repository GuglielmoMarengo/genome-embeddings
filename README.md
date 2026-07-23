# Genome Embeddings

[![Python](https://img.shields.io/badge/Python-3.13%2B-blue)](https://www.python.org/)
[![License](https://img.shields.io/github/license/GuglielmoMarengo/genome-embeddings)](LICENSE)
[![Last Commit](https://img.shields.io/github/last-commit/GuglielmoMarengo/genome-embeddings)](https://github.com/GuglielmoMarengo/genome-embeddings/commits/main)
[![Status](https://img.shields.io/badge/Status-In%20Development-orange)](https://github.com/GuglielmoMarengo/genome-embeddings)
[![Tests](https://img.shields.io/badge/Tests-140%20passing-brightgreen)](tests)
[![Security Policy](https://img.shields.io/badge/Security-Policy-blue)](SECURITY.md)
[![Contributions Welcome](https://img.shields.io/badge/Contributions-Welcome-brightgreen)](CONTRIBUTING.md)

**Turning genomes into mathematics.**

Genome Embeddings is an open-source Python project for representing genomic sequences through interpretable mathematical descriptors.

Rather than relying exclusively on neural networks, the project explores transparent and reproducible representations inspired by information theory, statistics, graph theory and number theory.

The current implementation supports:

* sequence-level mathematical descriptors;
* pairwise genome comparison;
* multi-genome distance and similarity matrices;
* matrix serialization;
* multiscale k-mer sensitivity analysis;
* matrix-geometry trajectories;
* cross-scale trajectory-change analysis;
* pair contributions to matrix deformation;
* full-dataset versus biological-only multiscale comparison;
* squared-deformation partitioning for synthetic-control effects;
* graphical analysis through heatmaps, distributions and trajectory plots;
* structured security, contribution, issue-reporting and pull-request workflows.

> This project is currently intended for research and software-development purposes. Any future clinical or diagnostic use would require extensive biological, statistical, clinical and regulatory validation.

---

## Project Goals

Genome Embeddings aims to develop reusable representations of genomic and transcriptomic sequences that can be:

* measured;
* interpreted;
* compared;
* ranked;
* converted into numerical vectors;
* exported into reusable formats;
* visualized;
* analyzed across multiple sequence scales;
* integrated with statistical and machine-learning workflows.

The long-term objective is to investigate whether interpretable mathematical descriptors can produce robust and biologically meaningful genomic embeddings without sacrificing transparency.

---

## Project Principles

The project follows these principles:

1. **Interpretability before complexity**
2. **Tests before or alongside implementation**
3. **Explicit mathematical definitions**
4. **Small and reviewable changes**
5. **Separation between analysis and visualization**
6. **Reproducible behavior**
7. **No unsupported biological or clinical claims**
8. **One descriptor family or analytical concept at a time**

---

## Current Capabilities

### Sequence analysis

* Single-record FASTA parsing
* DNA sequence validation
* Automatic uppercase normalization
* Sequence length
* Reverse complement
* Nucleotide frequencies
* GC and AT content
* GC skew
* Purine and pyrimidine content
* Shannon entropy

### k-mer analysis

* Configurable k-mer length
* Sliding-window k-mer extraction
* k-mer frequencies
* Normalized k-mer diversity
* k-mer entropy

### Structured representations

* `Genome`
* `GenomeDescriptor`
* `GenomeComparison`
* `GenomeCollection`
* `GenomeMatrix`

### Comparison tools

* Raw descriptor vectors
* Normalized descriptor vectors
* Euclidean distance
* Cosine similarity
* Feature-level differences
* Sorted feature differences
* Euclidean distance matrices
* Cosine similarity matrices
* Label-based matrix lookup
* Distance and similarity ranking

### Matrix conversion and serialization

* Row-oriented matrix conversion
* Dictionary conversion
* JSON serialization
* Configurable JSON indentation
* CSV serialization
* Configurable CSV delimiter
* Upper-triangle vectorization
* Upper-triangle label pairs
* Structured upper-triangle rows

### Multiscale analysis

* Multi-`k` Euclidean matrices
* Multi-`k` cosine matrices
* Ordered matrix collections
* Empty and duplicate `k` validation
* Matrix-geometry trajectories
* Pair-specific trajectories
* Signed pair-trajectory step differences
* Cross-scale matrix-trajectory distances
* Pair contributions to matrix deformation
* Full-dataset versus biological-only trajectory comparison
* Relative step-distance reduction after control removal
* Squared-deformation shares for control-associated and biological-only pairs
* Euclidean and cosine sensitivity analysis across `k`

### Visualization

* Euclidean distance heatmaps
* Cosine similarity heatmaps
* Pair-trajectory line plots
* Matrix-value distributions
* Multi-`k` distribution plots
* PNG export
* Automatic output-directory creation
* Headless visualization testing

### Repository governance

* Security policy
* Private vulnerability-reporting guidance
* Contribution guidelines
* Code of conduct
* Structured bug-report form
* Structured feature-request form
* Standardized pull-request template
* Public-issue restrictions for security vulnerabilities

---

## Installation

Clone the repository:

```bash
git clone https://github.com/GuglielmoMarengo/genome-embeddings.git
cd genome-embeddings
```

Create a virtual environment if desired:

```bash
python -m venv .venv
```

Activate it on Windows PowerShell:

```powershell
.venv\Scripts\Activate.ps1
```

Activate it on Linux or macOS:

```bash
source .venv/bin/activate
```

Install the dependencies:

```bash
python -m pip install -r requirements.txt
```

Run the demonstration:

```bash
python main.py
```

Run the test suite:

```bash
python -m pytest
```

For detailed output:

```bash
python -m pytest -v
```

The current test suite contains **140 passing tests**.

---

## Dependencies

The project currently uses:

```text
matplotlib>=3.10
pytest>=8.3
```

Matplotlib is isolated in the visualization module. The core mathematical classes do not depend directly on plotting code.

---

## Quick Start

### Create a genome

```python
from src.genome import Genome

genome = Genome(
    sequence="ACGTACGT",
    header=">Example sequence",
)
```

### Load a FASTA file

```python
genome = Genome.from_fasta(
    "data/fluorescent_proteins/aequorea_victoria_gfp_cds.fasta"
)
```

### Calculate descriptors

```python
genome.length()
genome.gc_content()
genome.at_content()
genome.shannon_entropy()
genome.gc_skew()
genome.purine_content()
genome.pyrimidine_content()
genome.kmer_frequencies(k=3)
genome.kmer_diversity(k=3)
genome.kmer_entropy(k=3)
```

### Generate a structured descriptor

```python
descriptor = genome.descriptor(k=3)

descriptor.to_dict()
descriptor.to_vector()
descriptor.to_normalized_vector()
```

---

## GenomeDescriptor

`GenomeDescriptor` stores the mathematical properties calculated from a sequence.

Current fields:

```text
length
gc_content
at_content
shannon_entropy
gc_skew
purine_content
pyrimidine_content
kmer_length
kmer_diversity
kmer_entropy
```

The normalized comparison vector contains six non-redundant features:

```text
GC content
normalized Shannon entropy
normalized GC skew
purine content
k-mer diversity
normalized k-mer entropy
```

Sequence length, AT content, pyrimidine content and k-mer length are excluded from the normalized comparison vector because they are scale-dependent, redundant or configuration-related.

---

## Descriptor Definitions

### GC content

```text
GC content = (G + C) / sequence length
```

### AT content

```text
AT content = (A + T) / sequence length
```

For valid DNA:

```text
GC content + AT content = 1
```

### Shannon entropy

```text
H = -Σ p(x) log₂ p(x)
```

The maximum nucleotide entropy for DNA is:

```text
2 bits
```

Normalization:

```text
normalized Shannon entropy = Shannon entropy / 2
```

### GC skew

```text
GC skew = (G - C) / (G + C)
```

Normalization from `[-1, 1]` to `[0, 1]`:

```text
normalized GC skew = (GC skew + 1) / 2
```

### Purine content

```text
Purine content = (A + G) / sequence length
```

### Pyrimidine content

```text
Pyrimidine content = (C + T) / sequence length
```

### k-mer diversity

```text
k-mer diversity =
distinct observed k-mers
────────────────────────────────
min(total observed k-mers, 4ᵏ)
```

### k-mer entropy

```text
Hₖ = -Σ p(k-mer) log₂ p(k-mer)
```

The theoretical maximum for DNA is:

```text
log₂(4ᵏ) = 2k
```

Normalization:

```text
normalized k-mer entropy = k-mer entropy / (2k)
```

---

## Pairwise Genome Comparison

Two descriptors can be compared directly:

```python
first_descriptor = first_genome.descriptor(k=3)
second_descriptor = second_genome.descriptor(k=3)

distance = first_descriptor.euclidean_distance(
    second_descriptor
)

similarity = first_descriptor.cosine_similarity(
    second_descriptor
)
```

A structured comparison can be generated with:

```python
comparison = first_descriptor.compare(
    second_descriptor
)
```

The resulting `GenomeComparison` contains:

```text
euclidean_distance
cosine_similarity
feature_differences
```

Feature differences can be ranked from largest to smallest:

```python
for feature_name, difference in (
    comparison.sorted_feature_differences()
):
    print(f"{feature_name}: {difference:.4f}")
```

Cosine similarity refers only to similarity within the implemented descriptor space. It must not be interpreted directly as sequence identity, protein identity or evolutionary relatedness.

---

## Multi-Genome Comparison

Multiple genomes can be managed through `GenomeCollection`:

```python
from src.genome import GenomeCollection

collection = GenomeCollection(
    [
        first_genome,
        second_genome,
        third_genome,
    ]
)
```

Generate descriptors:

```python
descriptors = collection.descriptors(k=3)
```

Generate comparison matrices:

```python
labels = [
    "Genome A",
    "Genome B",
    "Genome C",
]

euclidean_matrix = collection.euclidean_distance_matrix(
    labels=labels,
    k=3,
)

cosine_matrix = collection.cosine_similarity_matrix(
    labels=labels,
    k=3,
)
```

Both methods return a structured `GenomeMatrix`.

---

## GenomeMatrix

`GenomeMatrix` stores matrix values together with their analytical context:

```text
labels
values
metric
kmer_length
```

Example:

```python
matrix = collection.euclidean_distance_matrix(
    labels=labels,
    k=3,
)

print(matrix.metric)
print(matrix.kmer_length)
print(matrix.labels)
print(matrix.values)
```

### Label-based lookup

```python
distance = matrix.get_value(
    row_label="Genome A",
    column_label="Genome B",
)
```

Unknown labels raise a `ValueError`.

### Ranking

A matrix can rank every other genome relative to a selected reference:

```python
ranking = matrix.rank_by_label(
    label="Genome A",
)
```

For Euclidean matrices, smaller values are ranked first.

For cosine matrices, larger values are ranked first.

The reference genome is excluded from its own ranking.

---

## Matrix Conversion

### Row-oriented conversion

```python
rows = matrix.to_rows()
```

Example:

```python
[
    {
        "label": "Genome A",
        "values": [0.0, 0.12, 0.45],
    },
    {
        "label": "Genome B",
        "values": [0.12, 0.0, 0.39],
    },
]
```

### Dictionary conversion

```python
matrix_dict = matrix.to_dict()
```

Example:

```python
{
    "labels": ["Genome A", "Genome B"],
    "values": [
        [0.0, 0.12],
        [0.12, 0.0],
    ],
    "metric": "euclidean",
    "kmer_length": 3,
}
```

Conversion methods return copies of the underlying lists, preventing accidental modification of the original matrix.

---

## Matrix Serialization

### JSON

```python
json_output = matrix.to_json()
```

Formatted JSON can be generated with:

```python
json_output = matrix.to_json(
    indent=2,
)
```

Example:

```json
{
  "labels": [
    "Genome A",
    "Genome B"
  ],
  "values": [
    [
      0.0,
      0.12
    ],
    [
      0.12,
      0.0
    ]
  ],
  "metric": "euclidean",
  "kmer_length": 3
}
```

### CSV

```python
csv_output = matrix.to_csv()
```

Example:

```csv
label,Genome A,Genome B
Genome A,0.0,0.12
Genome B,0.12,0.0
```

A custom delimiter can be selected:

```python
csv_output = matrix.to_csv(
    delimiter=";",
)
```

The current serialization methods return strings. Direct file-writing utilities are planned separately.

---

## Upper-Triangle Representation

Distance and similarity matrices are symmetric. Their diagonal and lower triangle therefore contain redundant information.

For a matrix with `n` genomes, the number of unique pairwise comparisons is:

```text
n × (n - 1) / 2
```

For six genomes:

```text
6 × 5 / 2 = 15
```

### Vector representation

```python
vector = matrix.to_upper_triangle_vector()
```

Example:

```python
[0.2, 0.4, 0.3]
```

### Associated label pairs

```python
pairs = matrix.upper_triangle_pairs()
```

Example:

```python
[
    ("Genome A", "Genome B"),
    ("Genome A", "Genome C"),
    ("Genome B", "Genome C"),
]
```

The order of the pairs matches the order of the vector coordinates.

### Structured rows

```python
rows = matrix.to_upper_triangle_rows()
```

Example:

```python
[
    {
        "row_label": "Genome A",
        "column_label": "Genome B",
        "value": 0.2,
    },
    {
        "row_label": "Genome A",
        "column_label": "Genome C",
        "value": 0.4,
    },
    {
        "row_label": "Genome B",
        "column_label": "Genome C",
        "value": 0.3,
    },
]
```

This representation preserves both numerical values and their biological labels.

---

## Configurable k-mer Resolution

The demonstration program uses:

```python
DEFAULT_KMER_LENGTH = 3
```

The core API accepts other values:

```python
descriptor = genome.descriptor(k=4)

matrix = collection.euclidean_distance_matrix(
    labels=labels,
    k=4,
)
```

Different values of `k` capture different sequence scales:

| k    | Main interpretation                           |
| ---- | --------------------------------------------- |
| `1`  | Nucleotide composition                        |
| `2`  | Short-range nucleotide dependencies           |
| `3`  | Triplet-scale organization                    |
| `4+` | Increasingly specific local sequence patterns |

Triplet k-mers are not automatically equivalent to codons because the current implementation uses a sliding window and does not enforce a biological reading frame.

Higher values of `k` increase sparsity, computational cost and sensitivity to sequence length.

---

## Multi-k Matrix Analysis

A `GenomeCollection` can generate one matrix for every requested k-mer length.

### Euclidean matrices

```python
matrices = collection.euclidean_distance_matrices(
    labels=labels,
    k_values=[1, 2, 3, 4],
)
```

### Cosine matrices

```python
matrices = collection.cosine_similarity_matrices(
    labels=labels,
    k_values=[1, 2, 3, 4],
)
```

The result is an ordered mapping:

```python
{
    1: GenomeMatrix(...),
    2: GenomeMatrix(...),
    3: GenomeMatrix(...),
    4: GenomeMatrix(...),
}
```

The requested order of `k_values` is preserved.

The collection rejects:

* empty `k_values`;
* duplicate k-mer lengths;
* invalid individual values of `k`.

---

## Matrix-Geometry Trajectories

Each symmetric comparison matrix can be converted into the vector of its unique pairwise values.

Repeating this process for several values of `k` produces a trajectory in matrix space:

```python
trajectory = collection.euclidean_matrix_trajectory(
    labels=labels,
    k_values=[1, 2, 3, 4],
)
```

Example structure:

```python
{
    1: [0.10, 0.20, 0.30],
    2: [0.12, 0.25, 0.28],
    3: [0.15, 0.27, 0.24],
    4: [0.18, 0.31, 0.22],
}
```

Cosine trajectories are generated with:

```python
trajectory = collection.cosine_matrix_trajectory(
    labels=labels,
    k_values=[1, 2, 3, 4],
)
```

Each vector represents the global geometry of the dataset at one k-mer scale.

---

## Pair Trajectories

A specific genome pair can be followed across multiple values of `k`.

### Euclidean pair trajectory

```python
trajectory = collection.euclidean_pair_trajectory(
    labels=labels,
    row_label="Genome A",
    column_label="Genome B",
    k_values=[1, 2, 3, 4],
)
```

### Cosine pair trajectory

```python
trajectory = collection.cosine_pair_trajectory(
    labels=labels,
    row_label="Genome A",
    column_label="Genome B",
    k_values=[1, 2, 3, 4],
)
```

Example:

```python
{
    1: 0.0635,
    2: 0.0641,
    3: 0.0656,
    4: 0.0862,
}
```

Pair trajectories make it possible to investigate whether the relationship between two sequences is stable or scale-dependent.

---

## Cross-Scale Change Analysis

The project can quantify changes between consecutive k-mer scales.

The insertion order of `k_values` defines the trajectory order, so non-monotonic sequences such as `[3, 1, 2]` remain valid and are analyzed in that exact order.

### Pair-trajectory step differences

Signed changes between consecutive pair-trajectory values are calculated with:

```python
differences = (
    GenomeCollection
    .pair_trajectory_step_differences(
        {
            1: 0.0635,
            2: 0.0641,
            3: 0.0656,
            4: 0.0862,
        }
    )
)
```

Example:

```python
{
    (1, 2): 0.0006,
    (2, 3): 0.0015,
    (3, 4): 0.0206,
}
```

Positive values indicate an increase in distance or similarity. Negative values indicate a decrease. The metric determines the interpretation of the sign.

Metric-specific convenience methods are also available:

```python
collection.euclidean_pair_trajectory_step_differences(
    labels=labels,
    row_label="Genome A",
    column_label="Genome B",
    k_values=[1, 2, 3, 4],
)

collection.cosine_pair_trajectory_step_differences(
    labels=labels,
    row_label="Genome A",
    column_label="Genome B",
    k_values=[1, 2, 3, 4],
)
```

### Matrix-trajectory step distances

The Euclidean distance between consecutive upper-triangle vectors quantifies how much the global comparison geometry changes between scales:

```python
step_distances = (
    GenomeCollection
    .matrix_trajectory_step_distances(
        trajectory
    )
)
```

For consecutive vectors `Mₖ` and `Mₖ₊₁`:

```text
step distance = sqrt(Σ (Mₖ₊₁[i] - Mₖ[i])²)
```

Metric-specific methods can calculate these values directly from a collection:

```python
collection.euclidean_matrix_trajectory_step_distances(
    labels=labels,
    k_values=[1, 2, 3, 4],
)

collection.cosine_matrix_trajectory_step_distances(
    labels=labels,
    k_values=[1, 2, 3, 4],
)
```

### Pair contributions to matrix deformation

A global matrix change can be decomposed into the contribution of each unique genome pair:

```python
contributions = (
    GenomeCollection
    .matrix_trajectory_pair_contributions(
        labels=labels,
        trajectory=trajectory,
    )
)
```

Each transition contains rows with:

```text
row_label
column_label
difference
absolute_difference
```

Rows are sorted by decreasing absolute change, making it possible to identify which genome relationships drive the deformation of the full matrix.

Convenience methods are available for both metrics:

```python
collection.euclidean_matrix_pair_contributions(
    labels=labels,
    k_values=[1, 2, 3, 4],
)

collection.cosine_matrix_pair_contributions(
    labels=labels,
    k_values=[1, 2, 3, 4],
)
```

---

## Visualization

Visualization functions are implemented in:

```text
src/visualization.py
```

The visualization layer consumes existing `GenomeMatrix` objects and trajectories without modifying the core mathematical model.

### Matrix heatmaps

```python
from src.visualization import (
    plot_matrix_heatmap,
)

figure = plot_matrix_heatmap(
    euclidean_matrix
)
```

The function supports:

* Euclidean matrices;
* cosine matrices;
* automatic metric titles;
* genome labels on both axes;
* optional numerical annotations;
* configurable decimal places.

Example without annotations:

```python
figure = plot_matrix_heatmap(
    cosine_matrix,
    annotate=False,
)
```

### Pair-trajectory plots

```python
from src.visualization import (
    plot_pair_trajectory,
)

figure = plot_pair_trajectory(
    trajectory=euclidean_pair_trajectory,
    row_label="Aequorea GFP",
    column_label="Acropora GFP",
    metric="euclidean",
)
```

The x-axis represents k-mer length. The y-axis represents the selected distance or similarity metric.

### Matrix-value distributions

```python
from src.visualization import (
    plot_matrix_distribution,
)

figure = plot_matrix_distribution(
    euclidean_matrix
)
```

Only unique pairwise values from the upper triangle are included.

The diagonal and mirrored lower triangle are excluded.

### Multi-k distributions

```python
from src.visualization import (
    plot_trajectory_distributions,
)

figure = plot_trajectory_distributions(
    trajectory=euclidean_trajectory,
    metric="euclidean",
)
```

Each box represents the distribution of unique pairwise comparisons at one value of `k`.

### Saving figures

```python
from src.visualization import (
    save_figure,
)

saved_path = save_figure(
    figure=figure,
    output_path=(
        "outputs/euclidean_heatmap.png"
    ),
    dpi=150,
    close=True,
)
```

The destination directory is created automatically when needed.

---

## Generated Visualizations

Running:

```bash
python main.py
```

generates:

```text
outputs/
├── aequorea_acropora_cosine_trajectory.png
├── aequorea_acropora_euclidean_trajectory.png
├── cosine_distribution.png
├── cosine_heatmap.png
├── cosine_multi_k_distribution.png
├── euclidean_distribution.png
├── euclidean_heatmap.png
└── euclidean_multi_k_distribution.png
```

The `outputs/` directory is excluded through `.gitignore` because the images are generated artifacts.

Selected figures intended for permanent documentation may later be stored separately under:

```text
docs/images/
```

---

## Example Dataset

The repository contains fluorescent-protein sequences, two biological negative controls and one synthetic control.

| File                                      | Organism                   |          Region |     Length or role |
| ----------------------------------------- | -------------------------- | --------------: | -----------------: |
| `aequorea_victoria_gfp_cds.fasta`         | *Aequorea victoria*        |             CDS |             717 nt |
| `aequorea_victoria_gfp_mrna.fasta`        | *Aequorea victoria*        |   Complete mRNA |             922 nt |
| `acropora_millepora_gfp_cds.fasta`        | *Acropora millepora*       |             CDS |             696 nt |
| `discosoma_fp583_cds.fasta`               | *Discosoma* species        |             CDS |             678 nt |
| `staphylococcus_aureus_cata_cds.fasta`    | *Staphylococcus aureus*    |             CDS | Biological control |
| `saccharomyces_cerevisiae_tpi1_cds.fasta` | *Saccharomyces cerevisiae* |             CDS | Biological control |
| `periodic_sequence.fasta`                 | Synthetic control          | Repeated `ACGT` |             120 nt |

The primary comparisons use equivalent coding regions:

```text
CDS versus CDS
```

The complete *Aequorea victoria* mRNA is retained for future experiments but is not used in the primary CDS comparison.

The biological controls help evaluate whether the descriptor space separates function, taxonomy, composition or other broad sequence properties.

The synthetic periodic sequence is a methodological control and is not treated as a biological reference.

---

## Current Example Results

The demonstration compares:

* *Aequorea victoria* GFP CDS
* *Acropora millepora* GFP CDS
* *Discosoma* FP583 CDS
* *Staphylococcus aureus* `catA` CDS
* *Saccharomyces cerevisiae* `TPI1` CDS
* a synthetic periodic sequence

### Pairwise comparison at `k = 3`

For the two GFP coding sequences:

```text
Euclidean distance: 0.0656
Cosine similarity: 0.9996
```

Largest feature differences:

```text
gc_content: 0.0590
normalized_kmer_entropy: 0.0164
normalized_shannon_entropy: 0.0160
kmer_diversity: 0.0156
purine_content: 0.0052
normalized_gc_skew: 0.0047
```

### Ranking from *Aequorea victoria* GFP

Euclidean ranking:

```text
Acropora GFP: 0.0656
S. cerevisiae TPI1: 0.0842
Discosoma FP583: 0.1102
S. aureus catA: 0.2074
Periodic control: 1.1187
```

Cosine ranking:

```text
Acropora GFP: 0.9996
S. cerevisiae TPI1: 0.9992
Discosoma FP583: 0.9990
S. aureus catA: 0.9954
Periodic control: 0.8071
```

### Multiscale matrix trajectory

With six sequences, each matrix trajectory vector contains:

```text
15 dimensions
```

because:

```text
6 × 5 / 2 = 15
```

The first Euclidean coordinates across `k = 1, 2, 3, 4` are:

```text
k=1: [0.0635, 0.1088, 0.2011, ...]
k=2: [0.0641, 0.1097, 0.1988, ...]
k=3: [0.0656, 0.1102, 0.2074, ...]
k=4: [0.0862, 0.1216, 0.2719, ...]
```

The first cosine coordinates are:

```text
k=1: [0.9996, 0.9989, 0.9951, ...]
k=2: [0.9996, 0.9989, 0.9952, ...]
k=3: [0.9996, 0.9990, 0.9954, ...]
k=4: [0.9994, 0.9989, 0.9918, ...]
```

### *Aequorea*–*Acropora* pair trajectory

Euclidean:

```text
k=1: 0.0635
k=2: 0.0641
k=3: 0.0656
k=4: 0.0862
```

Cosine:

```text
k=1: 0.9996
k=2: 0.9996
k=3: 0.9996
k=4: 0.9994
```

### Pair-trajectory step differences

Euclidean:

```text
k=1 -> k=2: +0.000584
k=2 -> k=3: +0.001426
k=3 -> k=4: +0.020621
```

Cosine:

```text
k=1 -> k=2: -0.000004
k=2 -> k=3: +0.000038
k=3 -> k=4: -0.000237
```

### Matrix-trajectory step distances

Euclidean:

```text
k=1 -> k=2: 1.702979
k=2 -> k=3: 0.499185
k=3 -> k=4: 0.264115
```

Cosine:

```text
k=1 -> k=2: 0.233006
k=2 -> k=3: 0.185453
k=3 -> k=4: 0.035751
```

### Full dataset versus biological-only dataset

The full dataset contains six sequences and therefore 15 unique pairwise coordinates.

The biological-only dataset excludes the synthetic periodic control, contains five sequences and therefore 10 unique pairwise coordinates.

Euclidean matrix-step comparison:

```text
k=1 -> k=2:
full=1.702979
biological-only=0.003742
relative reduction=99.78%

k=2 -> k=3:
full=0.499185
biological-only=0.017912
relative reduction=96.41%

k=3 -> k=4:
full=0.264115
biological-only=0.152537
relative reduction=42.25%
```

Cosine matrix-step comparison:

```text
k=1 -> k=2:
full=0.233006
biological-only=0.000061
relative reduction=99.97%

k=2 -> k=3:
full=0.185453
biological-only=0.001222
relative reduction=99.34%

k=3 -> k=4:
full=0.035751
biological-only=0.008622
relative reduction=75.88%
```

The relative reduction is calculated as:

```text
(full step distance - biological-only step distance)
────────────────────────────────────────────────────
full step distance
```

This quantity measures the reduction in the global step distance after removing the periodic control. Because Euclidean norms are not additive, it is a comparative reduction rather than an exact contribution share.

### Squared-deformation partitioning

For an exact additive decomposition, each pairwise step difference is squared:

```text
total squared deformation = Σ pairwise difference²
```

The coordinates can then be divided into:

```text
control-associated pairs
+
biological-only pairs
=
total squared deformation
```

This partition quantifies the fraction of squared matrix deformation associated with pairs containing the synthetic periodic control.

The current demonstration calculates this partition separately for Euclidean and cosine trajectories. The calculation is exploratory and dataset-specific; it must not be interpreted as a general biological effect size.

---

## Interpretation

Within the current descriptor space:

* the synthetic periodic control remains strongly separated;
* the bacterial `catA` CDS is separated from the compact eukaryotic group;
* the eukaryotic `TPI1` control remains close to the fluorescent-protein CDS sequences;
* Euclidean distance provides stronger numerical separation than cosine similarity;
* cosine similarity remains highly compressed among biological sequences;
* the *Aequorea*–*Acropora* Euclidean relationship is stable for `k = 1, 2, 3` and changes more visibly at `k = 4`;
* the largest global matrix deformation occurs between `k = 1` and `k = 2`;
* many of the strongest cross-scale changes involve the synthetic periodic control;
* removing the periodic control reduces the first two Euclidean step distances by more than 96%;
* removing the periodic control reduces the first two cosine step distances by more than 99%;
* the biological-only Euclidean geometry changes more substantially between `k = 3` and `k = 4`;
* more specific changes between biological sequences begin to appear at higher `k`.

The visualizations make these patterns easier to inspect:

* heatmaps expose matrix structure;
* distribution plots show metric compression and separation;
* pair trajectories show scale-dependent changes;
* multi-`k` distributions show how global pairwise relationships vary with `k`.

The proximity of `TPI1` to the fluorescent-protein sequences shows that the current representation does not distinguish protein function.

The descriptors likely capture broad compositional, taxonomic or coding-sequence properties.

The cross-scale analysis distinguishes three complementary levels of change:

* pair step differences quantify how one selected relationship changes;
* matrix step distances quantify deformation of the complete dataset geometry;
* pair-contribution rankings identify which relationships drive that deformation;
* full versus biological-only comparisons show how strongly a synthetic control influences global geometry;
* squared-deformation partitioning provides an additive separation between control-associated and biological-only changes.

These results are preliminary and do not represent biological validation.

---

## Demonstration Program

`main.py` currently demonstrates:

* FASTA loading;
* genome summary generation;
* descriptor calculation;
* raw and normalized descriptor vectors;
* k-mer frequency calculation;
* pairwise genome comparison;
* Euclidean and cosine matrices;
* ranking relative to a reference genome;
* multi-`k` matrix generation;
* matrix-geometry trajectories;
* pair-specific trajectories;
* signed pair-trajectory step differences;
* cross-scale matrix distances;
* pair contributions to matrix deformation;
* full-dataset versus biological-only trajectory comparison;
* relative step-distance reduction after periodic-control removal;
* squared-deformation partitioning by pair category;
* heatmap generation;
* pair-trajectory plotting;
* matrix distributions;
* multi-`k` distributions;
* PNG export;
* JSON export preview;
* CSV export preview;
* label-based lookup;
* matrix metadata and row conversion.

The single-scale demonstration uses:

```python
DEFAULT_KMER_LENGTH = 3
```

The multiscale demonstration uses:

```python
KMER_SENSITIVITY_LENGTHS = [
    1,
    2,
    3,
    4,
]
```

---

## Architecture

```text
FASTA sequence
      │
      ▼
    Genome
      │
      ├── sequence descriptors
      ├── k-mer descriptors
      └── descriptor(k)
              │
              ▼
      GenomeDescriptor
              │
              ├── to_dict()
              ├── to_vector()
              ├── to_normalized_vector()
              ├── euclidean_distance()
              ├── cosine_similarity()
              ├── feature_differences()
              └── compare()
                       │
                       ▼
               GenomeComparison

Multiple Genome objects
      │
      ▼
GenomeCollection
      │
      ├── descriptors(k)
      ├── euclidean_distance_matrix()
      ├── cosine_similarity_matrix()
      ├── euclidean_distance_matrices()
      ├── cosine_similarity_matrices()
      ├── euclidean_matrix_trajectory()
      ├── cosine_matrix_trajectory()
      ├── euclidean_pair_trajectory()
      ├── cosine_pair_trajectory()
      ├── pair_trajectory_step_differences()
      ├── matrix_trajectory_step_distances()
      └── matrix_trajectory_pair_contributions()
                       │
                       ├── full-dataset analysis
                       ├── biological-only analysis
                       └── squared-deformation partitioning
                       │
                       ▼
                 GenomeMatrix
                       │
                       ├── get_value()
                       ├── rank_by_label()
                       ├── to_rows()
                       ├── to_dict()
                       ├── to_json()
                       ├── to_csv()
                       ├── to_upper_triangle_vector()
                       ├── upper_triangle_pairs()
                       └── to_upper_triangle_rows()
                               │
                               ▼
                       visualization.py
                               │
                               ├── plot_matrix_heatmap()
                               ├── plot_pair_trajectory()
                               ├── plot_matrix_distribution()
                               ├── plot_trajectory_distributions()
                               └── save_figure()
```

---

## Project Structure

```text
genome-embeddings/
├── .github/
│   ├── ISSUE_TEMPLATE/
│   │   ├── bug_report.yml
│   │   ├── feature_request.yml
│   │   └── config.yml
│   └── pull_request_template.md
├── data/
│   ├── fluorescent_proteins/
│   │   ├── aequorea_victoria_gfp_mrna.fasta
│   │   ├── aequorea_victoria_gfp_cds.fasta
│   │   ├── acropora_millepora_gfp_cds.fasta
│   │   └── discosoma_fp583_cds.fasta
│   └── controls/
│       ├── biological/
│       │   ├── staphylococcus_aureus_cata_cds.fasta
│       │   └── saccharomyces_cerevisiae_tpi1_cds.fasta
│       └── periodic_sequence.fasta
├── outputs/
├── src/
│   ├── genome.py
│   └── visualization.py
├── tests/
│   ├── data/
│   │   └── example.fasta
│   ├── test_genome.py
│   ├── test_multiscale_analysis.py
│   └── test_visualization.py
├── .gitignore
├── CODE_OF_CONDUCT.md
├── CONTRIBUTING.md
├── LICENSE
├── README.md
├── SECURITY.md
├── main.py
└── requirements.txt
```

`outputs/` is generated automatically and excluded from Git.

---

## Validation

The current implementation validates:

* non-empty DNA sequences;
* string input;
* accepted nucleotides: `A`, `C`, `G`, `T`;
* invalid nucleotide positions;
* FASTA headers;
* non-empty FASTA files;
* positive integer k-mer lengths;
* k-mer lengths not exceeding sequence length;
* descriptor comparison types;
* non-empty genome collections;
* collection item types;
* matrix dimensions;
* matrix label counts;
* supported matrix metrics;
* unknown matrix labels;
* empty multi-`k` requests;
* duplicate k-mer lengths;
* trajectories with fewer than two scales;
* inconsistent matrix-trajectory vector dimensions;
* matrix-vector and label-pair mismatches;
* duplicate or empty labels in contribution analysis;
* empty visualization trajectories;
* empty distribution vectors;
* invalid histogram bins;
* negative decimal places;
* invalid output DPI values.

---

## Testing Visualizations

Visualization tests use Matplotlib's non-interactive backend:

```python
matplotlib.use("Agg")
```

This allows figures to be created and saved during automated tests without opening graphical windows.

The visualization tests verify:

* figure creation;
* titles and axis labels;
* heatmap annotations;
* trajectory values;
* k-mer tick labels;
* empty-input validation;
* output-directory creation;
* PNG file creation;
* returned output paths.

---

## Contributing

Contributions are welcome.

Before contributing, read:

* [Contribution Guidelines](CONTRIBUTING.md)
* [Code of Conduct](CODE_OF_CONDUCT.md)
* [Security Policy](SECURITY.md)

Bug reports and feature proposals should use the structured GitHub issue forms.

Before submitting a pull request, run:

```bash
python -m pytest
python main.py
```

Generated files under `outputs/` must not be committed.

---

## Reporting Bugs

Use the repository's **Bug report** issue form.

Include:

* a minimal reproducible example;
* expected behavior;
* actual behavior;
* traceback or relevant output;
* Python version;
* operating system;
* sequence and k-mer characteristics.

Do not upload sensitive, proprietary, clinical or personally identifiable genomic data.

---

## Requesting Features

Use the repository's **Feature request** issue form.

Feature proposals should explain:

* the problem being solved;
* proposed behavior;
* mathematical or scientific motivation;
* interpretation and limitations;
* suggested tests;
* architecture and compatibility impact;
* relevant references.

---

## Security

Suspected security vulnerabilities must not be reported through public issues, discussions or pull requests.

Read the [Security Policy](SECURITY.md) and use GitHub Private Vulnerability Reporting when available.

Security-relevant concerns may include:

* unsafe handling of untrusted files;
* path traversal;
* unintended file access;
* arbitrary code execution;
* denial-of-service inputs;
* exposed credentials;
* vulnerable dependencies.

---

## Code of Conduct

Participation in Genome Embeddings is governed by the [Code of Conduct](CODE_OF_CONDUCT.md).

Technical and scientific disagreement is welcome when it remains respectful, evidence-based and focused on ideas rather than individuals.

---

## Roadmap

### Core representation

* [x] Core genome representation
* [x] Mathematical descriptors
* [x] Normalized descriptor vectors
* [x] Pairwise comparison
* [x] Structured comparison results
* [x] Multi-genome collections
* [x] Distance and similarity matrices
* [x] Biological negative controls

### Matrix utilities

* [x] Label-based matrix lookup
* [x] Similarity and distance ranking
* [x] Matrix row conversion
* [x] Matrix dictionary conversion
* [x] JSON serialization
* [x] CSV serialization
* [x] Upper-triangle vectorization
* [x] Upper-triangle label mapping
* [x] Structured pairwise rows
* [ ] Direct file-writing utilities

### Multiscale analysis

* [x] Multi-`k` Euclidean matrices
* [x] Multi-`k` cosine matrices
* [x] Initial k-mer sensitivity analysis
* [x] Matrix-geometry trajectories
* [x] Pair-specific trajectories
* [x] Pair-trajectory step differences
* [x] Cross-scale matrix distances
* [x] Pair contributions to matrix deformation
* [x] Full-dataset versus biological-only comparison
* [x] Relative step-distance reduction analysis
* [x] Synthetic-control squared-deformation partitioning
* [ ] Cross-scale stability metrics
* [ ] Ranking stability analysis
* [ ] Clustering stability analysis
* [ ] Sequence-length sensitivity analysis
* [ ] Sparse high-`k` representation
* [ ] Multiscale mutation signatures
* [ ] Multiscale genome embeddings
* [ ] Scale weighting
* [ ] Cross-scale normalization

### Visualization

* [x] Euclidean heatmap
* [x] Cosine heatmap
* [x] Pair-trajectory plots
* [x] Matrix-value distributions
* [x] Multi-`k` distribution plots
* [x] PNG export
* [x] Headless visualization testing
* [ ] Highlighted pairwise-value plots
* [ ] Full matrix-trajectory plots
* [ ] Cross-scale deformation plots
* [ ] Ranking-stability plots
* [ ] Clustering dendrograms
* [ ] Dimensionality-reduction plots
* [ ] Permanent README figures under `docs/images/`

### Repository governance and security

* [x] Security policy
* [x] Contribution guidelines
* [x] Code of conduct
* [x] Bug-report issue form
* [x] Feature-request issue form
* [x] Pull-request template
* [x] Private security-reporting link
* [ ] Dependabot version-update configuration
* [ ] Continuous-integration test workflow
* [ ] Automated code-quality checks

### Future descriptors and metrics

* [ ] Conditional entropy
* [ ] Block entropy
* [ ] Entropy rate
* [ ] Mutual information
* [ ] Jensen-Shannon divergence
* [ ] Dinucleotide statistics
* [ ] Full k-mer frequency vectors
* [ ] Hybrid descriptor vectors
* [ ] Codon-usage descriptors
* [ ] Reading-frame-aware descriptors
* [ ] Graph-based descriptors
* [ ] Spectral descriptors
* [ ] Compression-based descriptors
* [ ] Fractal descriptors
* [ ] Number-theoretic descriptors

### Platform development

* [ ] Multiple FASTA records
* [ ] Ambiguous nucleotide support
* [ ] RNA support
* [ ] Command-line interface
* [ ] Package distribution
* [ ] Descriptor export
* [ ] Comparison export
* [x] Matrix serialization
* [ ] Trajectory export
* [ ] Embedding export
* [ ] Parallel processing
* [ ] Biological benchmark datasets
* [ ] Statistical validation
* [ ] External validation
* [ ] Transcript and transcriptomic embeddings

---

## Multiscale Research Direction

The current implementation distinguishes between single-scale and multiscale analysis.

### Single-scale representation

A descriptor or comparison matrix is calculated for one selected value of `k`.

This provides a snapshot of sequence relationships at one resolution.

### k-mer sensitivity analysis

Several values of `k` are evaluated independently to study:

* descriptor stability;
* distance stability;
* similarity stability;
* ranking changes;
* clustering changes;
* sparsity;
* dependence on sequence length.

### Matrix-geometry trajectories

Each comparison matrix is transformed into the vector of its unique pairwise values.

As `k` changes, these vectors form a trajectory in a fixed-dimensional matrix space.

This representation may support the identification of:

* stable sequence relationships;
* scale-dependent relationships;
* abrupt geometric changes;
* sequence pairs driving matrix deformation;
* differences between distance metrics;
* exploratory multiscale mutation signatures;
* differences between full and biological-only dataset geometry;
* the fraction of squared deformation associated with synthetic-control pairs.

### Visual exploration

The visualization layer provides graphical views of:

* matrix structure;
* metric compression;
* pairwise-value distributions;
* changes across k-mer scales;
* pair-specific trajectories.

Visual evidence remains exploratory. The implemented step differences, matrix-trajectory distances, pair-contribution rankings, full-versus-biological comparisons and squared-deformation partitioning provide a quantitative layer beneath those visual patterns, but they still require statistical and biological validation.

### Future multiscale embeddings

Future work may combine features from several k-mer scales into one representation:

```text
global descriptors
        +
k = 1 features
        +
k = 2 features
        +
k = 3 features
        +
higher-order features
        =
multiscale genome embedding
```

The main challenge will be to combine scales without sacrificing interpretability.

Important research topics include:

* normalization;
* scale weighting;
* redundancy;
* sparsity;
* sequence-length dependence;
* statistical validation;
* biological validation.

---

## Research Questions

> Can interpretable mathematical descriptors produce biologically meaningful, robust and reusable embeddings of genomic and transcriptomic data?

> Can information across multiple sequence resolutions be combined into more informative representations without sacrificing interpretability?

> Can matrix-geometry trajectories reveal stable, scale-dependent or anomalous sequence relationships?

> Which sequence pairs drive changes in dataset geometry across k-mer scales?

> How much of the observed matrix deformation is associated with a synthetic control rather than relationships among biological sequences?

> Can graphical analysis expose meaningful multiscale patterns that can later be quantified through formal stability metrics?

> Can multiscale comparison geometries support an exploratory method for detecting and characterizing mutation signatures?

---

## Current Limitations

* Only DNA sequences containing `A`, `C`, `G` and `T` are supported.
* FASTA parsing currently supports a single record.
* The descriptor vector contains a limited number of aggregated features.
* Sliding-window triplets are not equivalent to reading-frame-aware codons.
* The example dataset is too small for biological validation.
* The synthetic periodic control dominates several global cross-scale changes.
* Full-versus-biological comparisons are currently demonstrated on one small dataset.
* Relative step-distance reduction is comparative and is not an additive contribution measure.
* Squared-deformation shares are exact for the current matrix coordinates but remain dataset-specific exploratory quantities.
* Cosine similarity is highly compressed in the current descriptor space.
* Histograms currently contain only 15 unique pairwise observations for the six-sequence dataset.
* Box plots summarize small exploratory samples.
* Visualization does not replace statistical testing.
* No statistical significance testing is currently implemented.
* No phylogenetic, structural or functional inference should be made from the current results.
* Multiscale trajectories are exploratory representations, not validated biological signatures.

---

## License

This project is released under the MIT License.

See the [LICENSE](LICENSE) file for details.