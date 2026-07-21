# Genome Embeddings

[![Python](https://img.shields.io/badge/Python-3.13%2B-blue)](https://www.python.org/)
[![License](https://img.shields.io/github/license/GuglielmoMarengo/genome-embeddings)](LICENSE)
[![Last Commit](https://img.shields.io/github/last-commit/GuglielmoMarengo/genome-embeddings)](https://github.com/GuglielmoMarengo/genome-embeddings/commits/main)
[![Status](https://img.shields.io/badge/Status-In%20Development-orange)](https://github.com/GuglielmoMarengo/genome-embeddings)

**Turning genomes into mathematics.**

Genome Embeddings is an open-source Python project for representing genomic sequences through interpretable mathematical descriptors.

Instead of relying exclusively on neural networks, the project explores transparent and reproducible representations inspired by information theory, statistics, graph theory and number theory.

> The project is currently intended for research and software-development purposes. Any future clinical or diagnostic application would require extensive biological, statistical, clinical and regulatory validation.

---

## Project Goals

Genome Embeddings aims to build reusable representations of genomic and transcriptomic sequences that can be:

* measured;
* interpreted;
* compared;
* converted into numerical vectors;
* analyzed across multiple sequence scales;
* integrated with statistical and machine-learning workflows.

The current implementation focuses on explainable DNA descriptors, pairwise comparison and multi-genome comparison matrices.

A major future direction is **multiscale representation**, combining information calculated across multiple k-mer lengths rather than relying on a single value of `k`.

---

## Implemented Features

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
* k-mer extraction
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
* Row-oriented matrix conversion
* Dictionary matrix conversion

---

## Installation

Clone the repository:

```bash
git clone https://github.com/GuglielmoMarengo/genome-embeddings.git
cd genome-embeddings
```

Install the dependencies:

```bash
pip install -r requirements.txt
```

Run the demonstration:

```bash
python main.py
```

Run the tests:

```bash
python -m pytest
```

For detailed test output:

```bash
python -m pytest -v
```

The current test suite contains **71 tests**.

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

### Calculate individual descriptors

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

The current fields are:

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

Sequence length, AT content, pyrimidine content and k-mer length are excluded from the normalized vector because they are scale-dependent, redundant or configuration-related.

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

The maximum nucleotide entropy for DNA is `2 bits`.

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

A complete structured comparison can be generated with:

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

Feature differences can be sorted from largest to smallest:

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

`GenomeMatrix` stores matrix values together with their context:

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

## Configurable k-mer Resolution

The demonstration program currently uses:

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

## Example Dataset

The repository contains a small exploratory dataset of fluorescent-protein sequences.

| File                               | Organism             |          Region | Length |
| ---------------------------------- | -------------------- | --------------: | -----: |
| `aequorea_victoria_gfp_cds.fasta`  | *Aequorea victoria*  |             CDS | 717 nt |
| `aequorea_victoria_gfp_mrna.fasta` | *Aequorea victoria*  |   Complete mRNA | 922 nt |
| `acropora_millepora_gfp_cds.fasta` | *Acropora millepora* |             CDS | 696 nt |
| `discosoma_fp583_cds.fasta`        | *Discosoma* species  |             CDS | 678 nt |
| `periodic_sequence.fasta`          | Synthetic control    | Repeated `ACGT` | 120 nt |

The primary comparisons use equivalent coding regions:

```text
CDS versus CDS
```

The complete *Aequorea victoria* mRNA is retained for future experiments but is not used in the primary CDS comparison.

The synthetic periodic sequence is a methodological control and is not treated as a biological reference.

---

## Current Example Results

The demonstration compares:

* *Aequorea victoria* GFP CDS
* *Acropora millepora* GFP CDS
* *Discosoma* FP583 CDS
* a synthetic periodic sequence

### Pairwise comparison

For the two GFP coding sequences at `k = 3`:

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

### Euclidean distance matrix

```text
                      Aequorea GFP   Acropora GFP   Discosoma FP583   Periodic control
Aequorea GFP                 0.0000         0.0656            0.1102             1.1187
Acropora GFP                 0.0656         0.0000            0.0795             1.1358
Discosoma FP583              0.1102         0.0795            0.0000             1.1393
Periodic control             1.1187         1.1358            1.1393             0.0000
```

### Cosine similarity matrix

```text
                      Aequorea GFP   Acropora GFP   Discosoma FP583   Periodic control
Aequorea GFP                 1.0000         0.9996            0.9990             0.8071
Acropora GFP                 0.9996         1.0000            0.9993             0.8110
Discosoma FP583              0.9990         0.9993            1.0000             0.8189
Periodic control             0.8071         0.8110            0.8189             1.0000
```

Within the current descriptor space:

* the fluorescent-protein coding sequences form a compact group;
* the periodic control is clearly separated;
* Euclidean distance provides stronger numerical separation than cosine similarity.

These observations are preliminary and do not represent biological validation.

---

## Architecture

```text
FASTA sequence
      │
      ▼
    Genome
      │
      ├── descriptor methods
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
      └── cosine_similarity_matrix()
                       │
                       ▼
                 GenomeMatrix
                       │
                       ├── get_value()
                       ├── to_rows()
                       └── to_dict()
```

---

## Project Structure

```text
genome-embeddings/
├── data/
│   ├── fluorescent_proteins/
│   │   ├── aequorea_victoria_gfp_mrna.fasta
│   │   ├── aequorea_victoria_gfp_cds.fasta
│   │   ├── acropora_millepora_gfp_cds.fasta
│   │   └── discosoma_fp583_cds.fasta
│   └── controls/
│       └── periodic_sequence.fasta
├── src/
│   └── genome.py
├── tests/
│   ├── data/
│   │   └── example.fasta
│   └── test_genome.py
├── main.py
├── README.md
├── requirements.txt
└── LICENSE
```

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
* unknown matrix lookup labels.

RNA and ambiguous nucleotide support are planned.

---

## Roadmap

### Near term

* [x] Core genome representation
* [x] Mathematical descriptors
* [x] Normalized descriptor vectors
* [x] Pairwise comparison
* [x] Structured comparison results
* [x] Multi-genome collections
* [x] Distance and similarity matrices
* [x] Label-based matrix lookup
* [x] Matrix row conversion
* [x] Matrix dictionary conversion
* [ ] Biological negative controls
* [ ] Similarity ranking
* [ ] Matrix export
* [ ] Heatmap visualization
* [ ] Clustering

### Multiscale representation

* [ ] k-mer sensitivity analysis
* [ ] Multi-k descriptor comparison
* [ ] Multi-k matrix comparison
* [ ] Cross-scale stability analysis
* [ ] Multiscale genome embeddings
* [ ] Scale weighting
* [ ] Cross-scale normalization
* [ ] Sequence-length sensitivity analysis
* [ ] Sparse high-k representation

### Future descriptors

* [ ] Conditional entropy
* [ ] Block entropy
* [ ] Entropy rate
* [ ] Mutual information
* [ ] Jensen-Shannon divergence
* [ ] Dinucleotide statistics
* [ ] Codon-usage descriptors
* [ ] Reading-frame-aware descriptors
* [ ] Graph-based descriptors
* [ ] Spectral descriptors
* [ ] Compression-based descriptors
* [ ] Fractal descriptors
* [ ] Number-theoretic descriptors

### Future platform development

* [ ] Multiple FASTA records
* [ ] Ambiguous nucleotide support
* [ ] RNA support
* [ ] Command-line interface
* [ ] Package distribution
* [ ] Descriptor and matrix export
* [ ] Visualization tools
* [ ] Parallel processing
* [ ] Biological benchmark datasets
* [ ] Statistical validation
* [ ] External validation
* [ ] Transcript and transcriptomic embeddings

---

## Multiscale Vision

The current implementation calculates one descriptor space for one selected k-mer length.

Future work will distinguish between:

### k-mer sensitivity analysis

Each value of `k` will be evaluated independently to study:

* descriptor stability;
* distance stability;
* similarity rankings;
* clustering stability;
* sparsity;
* dependence on sequence length.

### Multiscale embeddings

Features derived from multiple values of `k` will be combined into one representation:

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

The challenge will be to combine scales without sacrificing interpretability.

Key research topics will include normalization, weighting, redundancy, sparsity and biological validation.

---

## Research Question

> Can interpretable mathematical descriptors produce biologically meaningful, robust and reusable embeddings of genomic and transcriptomic data?

A complementary multiscale question is:

> Can information across multiple sequence resolutions be combined into more informative embeddings without sacrificing interpretability?

---

## License

This project is released under the MIT License.

See the [LICENSE](LICENSE) file for details.