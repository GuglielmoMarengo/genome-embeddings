# Genome Embeddings

[![Python](https://img.shields.io/badge/Python-3.13%2B-blue)](https://www.python.org/)
[![License](https://img.shields.io/github/license/GuglielmoMarengo/genome-embeddings)](LICENSE)
[![Last Commit](https://img.shields.io/github/last-commit/GuglielmoMarengo/genome-embeddings)](https://github.com/GuglielmoMarengo/genome-embeddings/commits/main)
[![Status](https://img.shields.io/badge/Status-In%20Development-orange)](https://github.com/GuglielmoMarengo/genome-embeddings)

**Turning genomes into mathematics.**

Genome Embeddings is an open-source Python project for representing genomic sequences through mathematically interpretable descriptors.

Rather than relying solely on machine learning, the project explores genomic representations inspired by information theory, number theory, graph theory and statistics.

---

# Vision

Most genomic embeddings rely on neural networks trained on large biological datasets.

Genome Embeddings explores a complementary direction: representing genomic sequences through transparent, reproducible and explainable mathematical descriptors.

The long-term goal is to build interpretable genome and transcriptome embeddings that can support downstream applications in:

* Bioinformatics
* Computational biology
* Comparative genomics
* Transcriptomics
* Genome comparison
* Clustering
* Classification
* Anomaly detection
* Biomarker research
* Machine learning
* Artificial intelligence

A central future direction of the project is multiscale sequence representation.

Rather than describing a biological sequence at a single k-mer resolution, future embeddings will combine information across multiple values of `k`.

Different k-mer scales capture different levels of sequence organization:

* `k = 1`: nucleotide composition
* `k = 2`: short-range nucleotide dependencies
* `k = 3`: triplet-scale and codon-related organization
* `k >= 4`: progressively more specific local sequence patterns

The long-term objective is to build embeddings that preserve information across complementary sequence scales while remaining mathematically interpretable.

The project is currently intended for research and software-development purposes.

Potential future clinical or diagnostic applications would require extensive biological, statistical, clinical and regulatory validation.

---

# Implemented Features

Current capabilities include:

* FASTA parsing
* DNA sequence validation
* Automatic uppercase normalization
* Sequence length calculation
* GC content
* AT content
* Reverse complement
* Nucleotide frequency analysis
* Shannon entropy
* GC skew
* Purine content
* Pyrimidine content
* k-mer extraction
* k-mer frequency analysis
* Normalized k-mer diversity
* k-mer entropy
* Configurable single k-mer length
* Genome descriptor generation
* Descriptor dictionary conversion
* Raw descriptor vector conversion
* Normalized descriptor vector conversion
* Euclidean distance between genome descriptors
* Cosine similarity between genome descriptors
* Feature-level comparison and interpretation
* Structured `GenomeComparison` results
* Sorted feature differences
* Comparison between real biological CDS sequences
* Curated fluorescent-protein example dataset
* Synthetic periodic control sequence
* `GenomeCollection` for multiple sequences
* Batch descriptor generation
* Euclidean distance matrices
* Cosine similarity matrices
* Structured `GenomeMatrix` results
* Matrix labels and metadata
* Matrix shape validation
* Matrix metric validation
* Label-based matrix value lookup
* Row-oriented matrix conversion
* Dictionary matrix conversion
* Defensive copies for converted matrix data
* Multi-genome comparison using real biological sequences

---

# Genome Representation

A genomic sequence is represented by the `Genome` class:

```python
from src.genome import Genome

genome = Genome(
    sequence="ACGTACGT",
    header=">Example sequence",
)
```

The class validates and normalizes the sequence when the object is created.

A genome can also be loaded directly from a FASTA file:

```python
genome = Genome.from_fasta(
    "data/fluorescent_proteins/aequorea_victoria_gfp_cds.fasta"
)
```

The current FASTA parser reads a single FASTA record, preserves its header and concatenates multiline sequence content.

---

# Mathematical Descriptors

The `Genome` class exposes individual descriptor methods that can be calculated and tested independently:

```python
genome.length()
genome.gc_content()
genome.at_content()
genome.nucleotide_frequencies()
genome.shannon_entropy()
genome.gc_skew()
genome.purine_content()
genome.pyrimidine_content()
genome.kmers(k=3)
genome.kmer_frequencies(k=3)
genome.kmer_diversity(k=3)
genome.kmer_entropy(k=3)
```

This separation keeps each mathematical calculation reusable and allows `descriptor()` to focus only on aggregating the results.

---

# GenomeDescriptor

The `GenomeDescriptor` object provides a structured mathematical representation of a genomic sequence.

A descriptor is created with:

```python
descriptor = genome.descriptor(k=3)
```

The current descriptor contains:

* `length`
* `gc_content`
* `at_content`
* `shannon_entropy`
* `gc_skew`
* `purine_content`
* `pyrimidine_content`
* `kmer_length`
* `kmer_diversity`
* `kmer_entropy`

The descriptor can be converted into a dictionary:

```python
descriptor.to_dict()
```

It can also be converted into a complete numerical vector:

```python
descriptor.to_vector()
```

For similarity and distance calculations, the descriptor provides a normalized vector:

```python
descriptor.to_normalized_vector()
```

The normalized vector excludes fields that are redundant, scale-dependent or configuration-related.

It currently contains:

* GC content
* Normalized Shannon entropy
* Normalized GC skew
* Purine content
* k-mer diversity
* Normalized k-mer entropy

This normalized six-dimensional vector is the mathematical foundation for genome comparison and future embedding construction.

---

# Descriptor Definitions

## Sequence Length

The total number of nucleotides in the sequence.

```text
length = number of nucleotides
```

Sequence length is included in the complete descriptor, but it is excluded from the normalized comparison vector because it has no universal maximum and could dominate distance calculations.

---

## GC Content

The proportion of guanine and cytosine nucleotides in the sequence.

```text
GC content = (G + C) / sequence length
```

---

## AT Content

The proportion of adenine and thymine nucleotides in the sequence.

```text
AT content = (A + T) / sequence length
```

For a valid DNA sequence:

```text
GC content + AT content = 1
```

AT content is retained in the complete descriptor but excluded from the normalized comparison vector because it is fully determined by GC content.

---

## Shannon Entropy

Shannon entropy measures the uncertainty of the nucleotide distribution.

```text
H = -Σ p(x) log₂ p(x)
```

For DNA, the maximum nucleotide entropy is `2 bits`, which occurs when:

```text
P(A) = P(C) = P(G) = P(T) = 0.25
```

A sequence containing only one nucleotide has an entropy of `0 bits`.

For normalized comparison:

```text
normalized Shannon entropy = Shannon entropy / 2
```

---

## GC Skew

GC skew measures the relative imbalance between guanine and cytosine.

```text
GC skew = (G - C) / (G + C)
```

GC skew ranges from `-1` to `1`.

If a sequence contains neither guanine nor cytosine, the GC skew is defined as `0.0`.

For normalized comparison:

```text
normalized GC skew = (GC skew + 1) / 2
```

This transforms the original range from `[-1, 1]` to `[0, 1]`.

---

## Purine Content

Purines are adenine and guanine.

```text
Purine content = (A + G) / sequence length
```

---

## Pyrimidine Content

Pyrimidines are cytosine and thymine.

```text
Pyrimidine content = (C + T) / sequence length
```

For a valid DNA sequence:

```text
Purine content + Pyrimidine content = 1
```

Pyrimidine content is retained in the complete descriptor but excluded from the normalized comparison vector because it is fully determined by purine content.

---

## k-mer Diversity

k-mer diversity measures how many distinct k-mers are observed relative to the maximum number that could be observed.

```text
k-mer diversity =
distinct observed k-mers
────────────────────────────
maximum observable k-mers
```

The maximum observable number of distinct k-mers is:

```text
min(total observed k-mers, 4ᵏ)
```

Therefore:

```text
k-mer diversity =
distinct observed k-mers
────────────────────────────────
min(total observed k-mers, 4ᵏ)
```

This normalization ensures that the value remains interpretable across sequences of different lengths.

A value close to `1.0` indicates that nearly all observable or theoretically possible k-mers are present.

A lower value indicates stronger repetition or reduced local sequence diversity.

---

## k-mer Entropy

k-mer entropy applies Shannon entropy to the frequency distribution of k-mers rather than individual nucleotides.

```text
Hₖ = -Σ p(k-mer) log₂ p(k-mer)
```

For DNA, the maximum possible k-mer entropy is:

```text
log₂(4ᵏ) = 2k
```

For normalized comparison:

```text
normalized k-mer entropy = k-mer entropy / (2k)
```

This descriptor captures local sequence organization that cannot be detected through nucleotide frequencies alone.

---

# Configurable k-mer Resolution

The current implementation uses `k = 3` as the default value in the demonstration program:

```python
DEFAULT_KMER_LENGTH = 3
```

This default provides a reproducible example, but the value is not fixed inside the core library.

Users can calculate descriptors and comparison matrices at different k-mer resolutions:

```python
descriptor = genome.descriptor(k=4)

distance_matrix = collection.euclidean_distance_matrix(
    labels=["Genome A", "Genome B"],
    k=4,
)

similarity_matrix = collection.cosine_similarity_matrix(
    labels=["Genome A", "Genome B"],
    k=4,
)
```

Different values of `k` capture different sequence properties:

* `k = 1` mainly describes nucleotide composition.
* `k = 2` captures dinucleotide-scale dependencies.
* `k = 3` captures triplet-scale organization.
* Higher values capture progressively more specific local patterns.

A triplet k-mer analysis does not automatically equal codon analysis because k-mers are currently extracted through a sliding window and do not necessarily follow a biological reading frame.

Increasing `k` also increases:

* The number of theoretically possible k-mers
* Feature sparsity
* Sensitivity to sequence length
* Computational cost
* The number of observations required for stable estimates

For these reasons, no single value of `k` should be assumed to be universally optimal.

---

# Genome Comparison

Genome descriptors can be compared through their normalized mathematical features.

The comparison methods operate on:

```python
descriptor.to_normalized_vector()
```

This avoids distortion caused by raw feature scales and excludes redundant dimensions.

---

## Euclidean Distance

Euclidean distance measures the geometric distance between two normalized descriptor vectors.

```text
d(A, B) = √Σ(Aᵢ - Bᵢ)²
```

It can be calculated with:

```python
distance = first_descriptor.euclidean_distance(second_descriptor)
```

A distance of `0.0` means that the two normalized descriptor vectors are identical.

Larger values indicate greater mathematical separation between the descriptor profiles.

Euclidean distance is symmetric:

```text
d(A, B) = d(B, A)
```

---

## Cosine Similarity

Cosine similarity measures the alignment between two normalized descriptor vectors.

```text
cosine similarity = (A · B) / (||A|| ||B||)
```

It can be calculated with:

```python
similarity = first_descriptor.cosine_similarity(second_descriptor)
```

A value of `1.0` indicates identical vector direction.

Lower values indicate increasingly different descriptor profiles.

Cosine similarity is symmetric:

```text
similarity(A, B) = similarity(B, A)
```

Cosine similarity must not be interpreted directly as a percentage of:

* Biological identity
* Sequence identity
* Protein identity
* Evolutionary relatedness
* Alignment similarity

It describes similarity only within the mathematical descriptor space currently implemented by the project.

---

## Feature Differences

A comparison can be explained feature by feature with:

```python
differences = first_descriptor.feature_differences(second_descriptor)
```

The method returns the absolute difference between corresponding normalized features:

```python
{
    "gc_content": 0.0590,
    "normalized_shannon_entropy": 0.0160,
    "normalized_gc_skew": 0.0047,
    "purine_content": 0.0052,
    "kmer_diversity": 0.0156,
    "normalized_kmer_entropy": 0.0164,
}
```

This makes each comparison interpretable by showing which mathematical properties distinguish the two sequences most strongly.

---

# GenomeComparison

The `GenomeComparison` object groups all pairwise comparison results into one structured representation.

It is created with:

```python
comparison = first_descriptor.compare(second_descriptor)
```

The object contains:

* `euclidean_distance`
* `cosine_similarity`
* `feature_differences`

Feature differences can be sorted from the largest to the smallest contribution:

```python
comparison.sorted_feature_differences()
```

Example:

```python
for feature_name, difference in comparison.sorted_feature_differences():
    print(f"{feature_name}: {difference:.4f}")
```

This architecture keeps comparison logic outside the presentation layer and allows comparison results to be reused by future ranking, export and visualization features.

---

# GenomeCollection

The `GenomeCollection` class manages multiple validated `Genome` objects.

A collection is created with:

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

The collection:

* Rejects non-list input
* Rejects empty collections
* Rejects values that are not `Genome` objects
* Preserves the original genome order
* Generates descriptors for all genomes
* Calculates Euclidean distance matrices
* Calculates cosine similarity matrices
* Returns structured `GenomeMatrix` objects

Descriptors are generated with:

```python
descriptors = collection.descriptors(k=3)
```

The returned descriptors preserve the same order as `collection.genomes`.

---

# GenomeMatrix

The `GenomeMatrix` object stores a multi-genome comparison matrix together with its metadata.

A matrix contains:

* `labels`
* `values`
* `metric`
* `kmer_length`

Example:

```python
matrix = collection.euclidean_distance_matrix(
    labels=[
        "Genome A",
        "Genome B",
        "Genome C",
    ],
    k=3,
)
```

The result is a `GenomeMatrix`:

```python
print(matrix.labels)
print(matrix.values)
print(matrix.metric)
print(matrix.kmer_length)
```

Example metadata:

```text
labels: ["Genome A", "Genome B", "Genome C"]
metric: "euclidean"
kmer_length: 3
```

The `values` field contains the numerical matrix:

```python
[
    [0.0, 0.12, 0.45],
    [0.12, 0.0, 0.39],
    [0.45, 0.39, 0.0],
]
```

`GenomeMatrix` validates that:

* Labels are not empty
* Matrix values are square
* The number of labels matches the matrix size
* The metric is supported

The currently supported metrics are:

* `euclidean`
* `cosine`

The object also provides:

* `get_value(row_label, column_label)`
* `to_rows()`
* `to_dict()`

These methods allow matrix data to be queried and converted without separating it from its labels and calculation metadata.

---

## Matrix Value Lookup

A single value can be retrieved by its row and column labels:

```python
distance = matrix.get_value(
    row_label="Genome A",
    column_label="Genome B",
)
```

The method resolves the requested labels and returns the corresponding matrix cell:

```python
print(distance)
```

```text
0.12
```

The lookup preserves row and column direction and does not assume that arbitrary matrices are symmetric.

If either label is unknown, the method raises:

```text
Unknown genome matrix label: Unknown.
```

---

## Row-Oriented Conversion

The matrix can be converted into a row-oriented representation:

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
    {
        "label": "Genome C",
        "values": [0.45, 0.39, 0.0],
    },
]
```

Each row contains:

* Its genome label
* A list of numerical values

The returned row values are copies.

Modifying the converted representation therefore does not modify the original `GenomeMatrix`.

---

## Dictionary Conversion

The complete matrix can be converted into a dictionary:

```python
matrix_dict = matrix.to_dict()
```

Example:

```python
{
    "labels": [
        "Genome A",
        "Genome B",
        "Genome C",
    ],
    "values": [
        [0.0, 0.12, 0.45],
        [0.12, 0.0, 0.39],
        [0.45, 0.39, 0.0],
    ],
    "metric": "euclidean",
    "kmer_length": 3,
}
```

The dictionary retains:

* Genome labels
* Numerical matrix values
* Metric name
* k-mer length

The returned labels and matrix rows are copies, preventing converted output from modifying the internal state of the original object.

This representation prepares the matrix for future JSON serialization and export functionality.

---

# Multi-Genome Comparison

A `GenomeCollection` can calculate pairwise metrics for every genome in the collection.

The labels must follow the same order as the genomes stored in the collection.

```python
labels = [
    "Genome A",
    "Genome B",
    "Genome C",
]
```

---

## Euclidean Distance Matrix

```python
distance_matrix = collection.euclidean_distance_matrix(
    labels=labels,
    k=3,
)
```

The resulting `GenomeMatrix` is:

* Square
* Symmetric
* Ordered consistently with `collection.genomes`
* Equal to `0.0` on the diagonal
* Associated with the `euclidean` metric
* Associated with the requested k-mer length

The numerical values are available through:

```python
distance_matrix.values
```

Example structure:

```text
                    Genome A   Genome B   Genome C
Genome A              0.0000     0.1200     0.4500
Genome B              0.1200     0.0000     0.3900
Genome C              0.4500     0.3900     0.0000
```

---

## Cosine Similarity Matrix

```python
similarity_matrix = collection.cosine_similarity_matrix(
    labels=labels,
    k=3,
)
```

The resulting `GenomeMatrix` is:

* Square
* Symmetric
* Ordered consistently with `collection.genomes`
* Equal to `1.0` on the diagonal
* Associated with the `cosine` metric
* Associated with the requested k-mer length

The numerical values are available through:

```python
similarity_matrix.values
```

Example structure:

```text
                    Genome A   Genome B   Genome C
Genome A              1.0000     0.9800     0.8100
Genome B              0.9800     1.0000     0.8300
Genome C              0.8100     0.8300     1.0000
```

---

# Real CDS-to-CDS Comparison

The primary pairwise example compares two real fluorescent-protein coding sequences:

* GFP CDS from *Aequorea victoria*
* GFP CDS from *Acropora millepora*

Using equivalent biological regions is important.

The project therefore compares:

```text
CDS versus CDS
```

rather than:

```text
complete mRNA versus CDS
```

This reduces differences caused solely by untranslated regions and makes the descriptor comparison more coherent.

Current example result for `k = 3`:

```text
Euclidean distance: 0.0656
Cosine similarity: 0.9996

Feature Differences:
gc_content: 0.0590
normalized_kmer_entropy: 0.0164
normalized_shannon_entropy: 0.0160
kmer_diversity: 0.0156
purine_content: 0.0052
normalized_gc_skew: 0.0047
```

These values indicate that the two coding sequences have highly similar profiles within the current six-dimensional descriptor space.

They do not establish sequence identity or evolutionary relatedness.

---

# Current Multi-Genome Result

The demonstration program compares:

* *Aequorea victoria* GFP CDS
* *Acropora millepora* GFP CDS
* *Discosoma* FP583 CDS
* A synthetic periodic control sequence

For `k = 3`, the Euclidean distance matrix is:

```text
                      Aequorea GFP   Acropora GFP   Discosoma FP583   Periodic control
Aequorea GFP                 0.0000         0.0656            0.1102             1.1187
Acropora GFP                 0.0656         0.0000            0.0795             1.1358
Discosoma FP583              0.1102         0.0795            0.0000             1.1393
Periodic control             1.1187         1.1358            1.1393             0.0000
```

The cosine similarity matrix is:

```text
                      Aequorea GFP   Acropora GFP   Discosoma FP583   Periodic control
Aequorea GFP                 1.0000         0.9996            0.9990             0.8071
Acropora GFP                 0.9996         1.0000            0.9993             0.8110
Discosoma FP583              0.9990         0.9993            1.0000             0.8189
Periodic control             0.8071         0.8110            0.8189             1.0000
```

Within the current descriptor space:

* The three real fluorescent-protein CDS sequences form a compact group.
* The synthetic periodic sequence is clearly separated.
* The Euclidean distance provides stronger numerical separation than cosine similarity.
* Cosine similarity remains very high between biological sequences because the current normalized descriptor dimensions are positive and relatively aggregated.

These observations are preliminary and should not be interpreted as biological validation.

The next validation step will require real non-fluorescent biological control sequences.

---

# Multiscale Representation

The project currently calculates descriptors using one configurable k-mer length at a time.

Future development will distinguish between two related concepts:

1. k-mer sensitivity analysis
2. Multiscale embeddings

---

## k-mer Sensitivity Analysis

Sensitivity analysis will calculate descriptors and comparison results independently across multiple values of `k`.

Conceptually:

```python
k_values = [1, 2, 3, 4, 5]
```

Each value will produce an independent descriptor space and independent comparison results.

This will make it possible to study:

* How descriptor values change with k-mer resolution
* How stable pairwise distances remain across scales
* Which values of `k` best separate biological groups
* How sequence length affects higher-order k-mer statistics
* Whether similarity rankings remain stable
* Whether clustering remains stable across parameter choices
* At which values sparsity begins to dominate
* How robust the representation is to methodological choices

Sensitivity analysis will not combine the scales.

Its purpose will be to measure how the method behaves when `k` changes.

---

## Multiscale Embeddings

A multiscale embedding will combine information derived from several k-mer lengths into a single mathematical representation.

Conceptually:

```text
global descriptors
        +
k = 1 descriptors
        +
k = 2 descriptors
        +
k = 3 descriptors
        +
higher-order k-mer descriptors
        =
multiscale genome embedding
```

A future API may resemble:

```python
embedding = genome.multiscale_embedding(
    k_values=[1, 2, 3, 4, 5],
)
```

Sensitivity analysis evaluates each scale separately.

Multiscale embedding combines complementary scales into one representation.

Potential advantages include:

* Preservation of both global and local sequence information
* Reduced dependence on a single parameter choice
* Improved discrimination between biologically similar sequences
* Better characterization of hierarchical sequence organization
* More robust representations across different sequence types

Future work will need to address:

* Feature normalization across scales
* Scale weighting
* Redundant information
* Dimensionality growth
* Increasing sparsity at larger values of `k`
* Dependence on sequence length
* Computational cost
* Biological interpretation
* Statistical validation

The goal is not to assume that one k-mer length is universally optimal, but to determine how information from different sequence resolutions can be combined transparently and reproducibly.

---

# Quick Example

```python
from pathlib import Path

from src.genome import Genome, GenomeCollection


PROJECT_ROOT = Path(__file__).resolve().parent

data_dir = PROJECT_ROOT / "data"

genomes = [
    Genome.from_fasta(
        data_dir
        / "fluorescent_proteins"
        / "aequorea_victoria_gfp_cds.fasta"
    ),
    Genome.from_fasta(
        data_dir
        / "fluorescent_proteins"
        / "acropora_millepora_gfp_cds.fasta"
    ),
    Genome.from_fasta(
        data_dir
        / "fluorescent_proteins"
        / "discosoma_fp583_cds.fasta"
    ),
    Genome.from_fasta(
        data_dir
        / "controls"
        / "periodic_sequence.fasta"
    ),
]

labels = [
    "Aequorea GFP",
    "Acropora GFP",
    "Discosoma FP583",
    "Periodic control",
]

collection = GenomeCollection(genomes)

descriptors = collection.descriptors(k=3)

pairwise_comparison = descriptors[0].compare(descriptors[1])

print(
    f"Euclidean distance: "
    f"{pairwise_comparison.euclidean_distance:.4f}"
)

print(
    f"Cosine similarity: "
    f"{pairwise_comparison.cosine_similarity:.4f}"
)

for feature_name, difference in (
    pairwise_comparison.sorted_feature_differences()
):
    print(f"{feature_name}: {difference:.4f}")

euclidean_matrix = collection.euclidean_distance_matrix(
    labels=labels,
    k=3,
)

cosine_matrix = collection.cosine_similarity_matrix(
    labels=labels,
    k=3,
)

selected_distance = euclidean_matrix.get_value(
    row_label="Aequorea GFP",
    column_label="Acropora GFP",
)

matrix_rows = euclidean_matrix.to_rows()
matrix_dict = euclidean_matrix.to_dict()

print(f"Selected distance: {selected_distance:.4f}")
print(f"First row: {matrix_rows[0]}")
print(f"Metric: {matrix_dict['metric']}")
print(f"k-mer length: {matrix_dict['kmer_length']}")
print(f"Labels: {matrix_dict['labels']}")

print(cosine_matrix.values)
```

---

# Running the Example

Run the demonstration program from the project root:

```bash
python main.py
```

The program:

1. Loads the GFP coding sequence from *Aequorea victoria*.
2. Loads the GFP coding sequence from *Acropora millepora*.
3. Loads the FP583 coding sequence from *Discosoma*.
4. Loads the synthetic periodic control sequence.
5. Creates a `GenomeCollection`.
6. Generates a mathematical descriptor for the reference CDS.
7. Prints basic sequence information.
8. Prints the complete descriptor.
9. Prints the raw descriptor vector.
10. Prints the normalized descriptor vector.
11. Displays the first observed k-mer frequencies.
12. Creates a detailed pairwise `GenomeComparison`.
13. Prints Euclidean distance and cosine similarity.
14. Prints feature differences in descending order.
15. Creates a structured Euclidean `GenomeMatrix`.
16. Creates a structured cosine `GenomeMatrix`.
17. Prints both matrices with their labels and k-mer resolution.
18. Retrieves a selected distance through `get_value()`.
19. Converts the Euclidean matrix through `to_rows()`.
20. Converts the Euclidean matrix through `to_dict()`.
21. Displays one converted row and the matrix metadata.

---

# Roadmap

## Genome representation

* [x] Genome class
* [x] FASTA parser
* [x] Sequence validation
* [x] Sequence normalization
* [x] Sequence length
* [x] Reverse complement

## Mathematical descriptors

* [x] Nucleotide frequencies
* [x] GC content
* [x] AT content
* [x] Shannon entropy
* [x] GC skew
* [x] Purine content
* [x] Pyrimidine content
* [x] k-mer extraction
* [x] k-mer frequencies
* [x] Normalized k-mer diversity
* [x] k-mer entropy
* [x] Configurable single k-mer length
* [x] GenomeDescriptor object
* [x] Descriptor dictionary conversion
* [x] Raw descriptor vector conversion
* [x] Normalized descriptor vector conversion

## Pairwise genome comparison

* [x] Descriptor normalization
* [x] Euclidean distance
* [x] Cosine similarity
* [x] Feature-level comparison
* [x] Explainable feature differences
* [x] GenomeComparison object
* [x] Sorted feature differences
* [x] Real CDS-to-CDS comparison

## Multi-genome comparison

* [x] GenomeCollection object
* [x] Collection validation
* [x] Batch descriptor generation
* [x] Descriptor-order preservation
* [x] Euclidean distance matrix
* [x] Cosine similarity matrix
* [x] Matrix symmetry validation
* [x] Matrix diagonal validation
* [x] GenomeMatrix object
* [x] Matrix labels
* [x] Matrix metric metadata
* [x] Matrix k-mer metadata
* [x] Matrix shape validation
* [x] Matrix label-count validation
* [x] Real multi-genome demonstration
* [x] Matrix value lookup
* [x] Matrix row conversion
* [x] Matrix dictionary conversion
* [x] Defensive matrix conversion copies
* [ ] Matrix export
* [ ] Biological negative controls
* [ ] Similarity ranking
* [ ] Clustering
* [ ] Heatmap visualization

## Multiscale representation

* [ ] k-mer sensitivity analysis
* [ ] Multi-k descriptor comparison
* [ ] Multi-k matrix comparison
* [ ] Cross-scale distance stability analysis
* [ ] Cross-scale similarity stability analysis
* [ ] Multiscale k-mer descriptors
* [ ] Multiscale genome embeddings
* [ ] Scale weighting
* [ ] Cross-scale normalization
* [ ] Cross-scale feature redundancy analysis
* [ ] Sequence-length sensitivity analysis
* [ ] Sparse high-k representation

## Future descriptors

* [ ] Conditional entropy
* [ ] Block entropy
* [ ] Entropy rate
* [ ] Mutual information
* [ ] Jensen-Shannon divergence
* [ ] Dinucleotide statistics
* [ ] Codon usage descriptors
* [ ] Reading-frame-aware descriptors
* [ ] Local window descriptors
* [ ] Graph-based descriptors
* [ ] Spectral descriptors
* [ ] Compression-based descriptors
* [ ] Fractal descriptors
* [ ] Number-theoretic descriptors

## Future platform development

* [ ] Genome embeddings
* [ ] Transcript embeddings
* [ ] Transcriptomic sample embeddings
* [ ] Configurable embedding parameters
* [ ] Embedding configuration object
* [ ] Parameter provenance
* [ ] Reproducible multiscale analysis
* [ ] Multiple FASTA record support
* [ ] Ambiguous nucleotide support
* [ ] RNA support
* [ ] Command-line interface
* [ ] Package distribution
* [ ] Descriptor export
* [ ] Comparison export
* [ ] Matrix export
* [ ] Embedding export
* [ ] Visualization tools
* [ ] Parallel processing
* [ ] Biological benchmark datasets
* [ ] Statistical validation
* [ ] External validation

---

# Genome Validation Rules

| Rule                                                      | Status     |
| --------------------------------------------------------- | ---------- |
| Sequence cannot be empty                                  | ✅          |
| Sequence must be a string                                 | ✅          |
| Sequence is automatically converted to uppercase          | ✅          |
| Only **A**, **C**, **G** and **T** are accepted           | ✅          |
| Invalid nucleotide positions are reported                 | ✅          |
| FASTA file cannot be empty                                | ✅          |
| FASTA header must begin with `>`                          | ✅          |
| k must be an integer                                      | ✅          |
| k must be greater than zero                               | ✅          |
| k cannot exceed the sequence length                       | ✅          |
| Descriptor comparisons require another `GenomeDescriptor` | ✅          |
| GenomeCollection input must be a list                     | ✅          |
| GenomeCollection cannot be empty                          | ✅          |
| GenomeCollection accepts only `Genome` objects            | ✅          |
| GenomeMatrix labels cannot be empty                       | ✅          |
| GenomeMatrix values must be square                        | ✅          |
| GenomeMatrix labels must match matrix size                | ✅          |
| GenomeMatrix metric must be supported                     | ✅          |
| Matrix order follows collection order                     | ✅          |
| Unknown lookup labels are rejected                        | ✅          |
| Converted matrix rows are copied                          | ✅          |
| Converted matrix dictionaries are copied                  | ✅          |
| Euclidean matrix diagonal equals `0.0`                    | ✅          |
| Cosine matrix diagonal equals `1.0`                       | ✅          |
| Ambiguous nucleotide support                              | 🚧 Planned |
| RNA support                                               | 🚧 Planned |

---

# Installation

Clone the repository:

```bash
git clone https://github.com/GuglielmoMarengo/genome-embeddings.git
```

Move into the project directory:

```bash
cd genome-embeddings
```

Install the dependencies:

```bash
pip install -r requirements.txt
```

---

# Running the Tests

Run the complete test suite from the project root:

```bash
python -m pytest
```

For more detailed output:

```bash
python -m pytest -v
```

The current test suite contains **71 tests** and covers:

* Sequence validation
* Sequence length
* GC content
* AT content
* Reverse complement
* FASTA parsing
* k-mer extraction
* k-mer frequencies
* Nucleotide frequencies
* Shannon entropy
* GC skew
* Purine content
* Pyrimidine content
* Descriptor properties
* k-mer diversity
* k-mer entropy
* GenomeDescriptor generation
* Raw descriptor vector conversion
* Normalized descriptor vector conversion
* Euclidean distance
* Euclidean distance symmetry
* Euclidean distance type validation
* Cosine similarity
* Cosine similarity symmetry
* Cosine similarity range
* Cosine similarity type validation
* Feature-level differences
* Feature-difference structure
* GenomeComparison generation
* GenomeComparison type validation
* Sorted feature differences
* GenomeCollection creation
* GenomeCollection input validation
* Empty collection validation
* Collection item validation
* Batch descriptor generation
* Descriptor-order preservation
* Invalid collection k-mer length handling
* GenomeMatrix creation
* GenomeMatrix metadata
* Empty matrix-label validation
* Square-matrix validation
* Matrix label-count validation
* Matrix metric validation
* Structured Euclidean matrix generation
* Euclidean matrix shape
* Euclidean matrix diagonal
* Euclidean matrix symmetry
* Euclidean matrix pairwise consistency
* Structured cosine matrix generation
* Cosine matrix shape
* Cosine matrix diagonal
* Cosine matrix symmetry
* Cosine matrix pairwise consistency
* Collection matrix label-count validation
* Label-based matrix value lookup
* Direction-preserving matrix lookup
* Unknown row-label validation
* Unknown column-label validation
* Row-oriented matrix conversion
* Row-conversion defensive copies
* Dictionary matrix conversion
* Dictionary-conversion defensive copies

---

# Project Structure

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

# Architecture

The project currently separates sequence analysis into six responsibilities.

## Genome

The `Genome` class:

* Stores the validated sequence
* Preserves an optional FASTA header
* Loads FASTA data
* Calculates individual mathematical descriptors
* Extracts and counts k-mers
* Produces a `GenomeDescriptor`

## GenomeDescriptor

The `GenomeDescriptor` class:

* Stores the calculated descriptor values
* Converts descriptors into dictionaries
* Converts descriptors into raw numerical vectors
* Produces normalized comparison vectors
* Calculates Euclidean distance
* Calculates cosine similarity
* Explains comparisons through feature differences
* Produces a structured `GenomeComparison`

## GenomeComparison

The `GenomeComparison` class:

* Stores Euclidean distance
* Stores cosine similarity
* Stores normalized feature differences
* Sorts feature differences from largest to smallest

## GenomeCollection

The `GenomeCollection` class:

* Stores multiple validated `Genome` objects
* Preserves genome order
* Generates descriptors for every genome
* Calculates Euclidean distance matrices
* Calculates cosine similarity matrices
* Produces structured `GenomeMatrix` objects

## GenomeMatrix

The `GenomeMatrix` class:

* Stores matrix labels
* Stores numerical matrix values
* Stores the comparison metric
* Stores the k-mer length
* Validates matrix dimensions
* Validates label count
* Validates the selected metric
* Retrieves values through row and column labels
* Converts the matrix into row-oriented data
* Converts the matrix into a dictionary
* Protects internal state through defensive copies
* Prepares matrix data for future export and visualization

## main.py

The demonstration program:

* Loads three real fluorescent-protein CDS records
* Loads one synthetic periodic control
* Creates a `GenomeCollection`
* Presents a detailed pairwise comparison
* Generates structured multi-genome matrices
* Prints matrix labels, values, metrics and k-mer resolution
* Demonstrates label-based value lookup
* Demonstrates row-oriented conversion
* Demonstrates dictionary conversion

The resulting architecture is:

```text
FASTA sequences
       │
       ▼
     Genome
       │
       ├── Individual descriptor methods
       │
       ▼
  descriptor(k)
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
                │
                └── sorted_feature_differences()

Multiple Genome objects
       │
       ▼
GenomeCollection
       │
       ├── descriptors(k)
       ├── euclidean_distance_matrix(labels, k)
       └── cosine_similarity_matrix(labels, k)
                          │
                          ▼
                    GenomeMatrix
                          │
                          ├── labels
                          ├── values
                          ├── metric
                          ├── kmer_length
                          ├── get_value()
                          ├── to_rows()
                          └── to_dict()
```

The next architectural steps may extend `GenomeMatrix` with:

* Matrix export
* Similarity ranking
* Nearest-neighbour retrieval
* Heatmap preparation
* Future support for multiple k-mer lengths

---

# Why This Project?

Most genomic embedding systems rely on complex neural networks trained on massive datasets.

Those models can be powerful, but their internal representations are often difficult to interpret.

Genome Embeddings explores a different approach: describing genomic and transcriptomic sequences through mathematical properties whose meaning remains explicit.

Rather than treating DNA or RNA as raw text, the project treats biological sequences as mathematical objects whose properties can be:

* Measured
* Interpreted
* Compared
* Explained
* Visualized
* Converted into vectors
* Organized into structured comparison matrices
* Queried through biological labels
* Converted into reusable data structures
* Analyzed across multiple scales
* Integrated with statistical models
* Integrated with machine-learning systems

The project emphasizes:

* Mathematical interpretability
* Biological interpretability
* Reproducibility
* Explainable descriptors
* Explainable comparisons
* Structured comparison results
* Multiscale sequence representation
* Modular architecture
* Automated testing
* Real biological datasets
* Comparisons between equivalent biological regions
* Extensibility
* Compatibility with future machine-learning applications

The central research question is:

> Can interpretable mathematical descriptors produce biologically meaningful, robust and reusable embeddings of genomic and transcriptomic data?

A complementary multiscale question is:

> Can information across multiple sequence resolutions be combined into more informative embeddings without sacrificing interpretability?

---

# Example Dataset

The repository contains a small curated dataset of fluorescent-protein nucleotide sequences.

## Fluorescent proteins

### `aequorea_victoria_gfp_cds.fasta`

* **Organism:** *Aequorea victoria*
* **Protein:** Green fluorescent protein
* **Gene:** GFP
* **Accession:** L29345.1
* **Protein ID:** AAA58246.1
* **Region:** CDS
* **Location:** 26..742
* **Sequence length:** 717 nt

This CDS is used as the primary reference sequence in `main.py`.

---

### `aequorea_victoria_gfp_mrna.fasta`

* **Organism:** *Aequorea victoria*
* **Protein:** Green fluorescent protein
* **Accession:** L29345.1
* **Region:** Complete mRNA record
* **Sequence length:** 922 nt

This file is retained as a reference record but is not used for the primary CDS-to-CDS comparison.

---

### `acropora_millepora_gfp_cds.fasta`

* **Organism:** *Acropora millepora*
* **Protein:** Green fluorescent protein
* **Accession:** AY646067.1
* **Protein ID:** AAU06846.1
* **Region:** CDS
* **Location:** 1..696
* **Sequence length:** 696 nt

This CDS is used as the primary pairwise comparison sequence in `main.py`.

---

### `discosoma_fp583_cds.fasta`

* **Organism:** *Discosoma* species
* **Protein:** Fluorescent protein FP583
* **Accession:** AF168419.2
* **Protein ID:** AAF03369.1
* **Region:** CDS
* **Location:** 54..731
* **Sequence length:** 678 nt

This CDS is included in the multi-genome distance and similarity matrices.

---

## Synthetic control

### `periodic_sequence.fasta`

This file contains a balanced artificial sequence constructed by repeating:

```text
ACGT
```

The current control contains 120 nucleotides.

It is retained as a methodological control for evaluating how descriptors respond to highly periodic local sequence organization.

The synthetic sequence is included in the multi-genome matrices but is not treated as a biological reference.

---

# Dataset Design

The biological comparisons use equivalent coding regions:

```text
Aequorea victoria GFP CDS
Acropora millepora GFP CDS
Discosoma FP583 CDS
```

This is more rigorous than mixing complete mRNA and CDS records because untranslated regions could affect:

* GC content
* Nucleotide entropy
* k-mer diversity
* k-mer entropy
* Descriptor distances

The retained complete mRNA record remains useful for future experiments exploring how genomic region selection affects mathematical descriptors.

The current dataset is intentionally small and exploratory.

Future validation will require:

* Non-fluorescent biological controls
* Unrelated CDS sequences
* Homologous and non-homologous sequence sets
* Larger gene families
* Taxonomically diverse datasets
* Perturbed and synthetic benchmark sequences
* Independent biological annotations
* Multi-k benchmark analyses
* Sequence-length sensitivity experiments

---

# Long-Term Vision

The project is divided conceptually into six phases.

## Phase 1: Core mathematical representation

The first phase builds reusable and interpretable measurements of genomic sequences.

Current examples include:

* Base composition
* Information entropy
* Sequence imbalance
* k-mer diversity
* k-mer entropy
* Normalized descriptor vectors
* Configurable single-scale k-mer analysis

## Phase 2: Explainable genome comparison

The second phase compares descriptor vectors through mathematical metrics.

Current capabilities include:

* Euclidean distance
* Cosine similarity
* Feature-level difference analysis
* Structured pairwise comparison objects
* Structured matrix comparison objects
* Label-based matrix lookup
* Reusable matrix conversions
* Real CDS-to-CDS comparison
* Multi-genome distance matrices
* Multi-genome similarity matrices

## Phase 3: Multiscale and advanced mathematical descriptors

This phase will investigate both independent scale analysis and combined multiscale representations.

Future capabilities may include:

* k-mer sensitivity analysis
* Multi-k descriptor comparison
* Multiscale genome embeddings
* Scale weighting
* Cross-scale normalization
* Conditional entropy
* Block entropy
* Entropy rate
* Mutual information
* Jensen-Shannon divergence
* Graph-based descriptors
* Spectral descriptors
* Compression-based descriptors
* Fractal descriptors
* Number-theoretic descriptors
* Local and multiscale descriptors

## Phase 4: Dataset-scale analysis

This phase will include:

* Matrix export
* Similarity ranking
* Nearest-neighbour retrieval
* Clustering
* Heatmaps
* Dimensionality reduction
* Batch processing
* Parallel execution
* Benchmark datasets
* Cross-scale stability analysis

## Phase 5: Transcriptomics

This phase may introduce:

* Transcript embeddings
* Isoform embeddings
* Expression-weighted embeddings
* Sample-level transcriptomic representations
* Co-expression graph descriptors
* Transcriptomic clustering
* Differential representation analysis
* Multiscale transcript representations

## Phase 6: Biological and translational validation

This phase will require:

* Biological benchmarks
* Statistical validation
* Ablation studies
* Robustness analysis
* Parameter sensitivity analysis
* External validation datasets
* Clinical research cohorts
* Uncertainty estimation
* Reproducible pipelines
* Auditability
* Domain-expert review

Any potential clinical or diagnostic application would require substantially more validation than the research-oriented functionality currently implemented.

---

# License

This project is released under the MIT License.

See the [LICENSE](LICENSE) file for details.