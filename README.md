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

The long-term goal is to build interpretable genome embeddings that can support downstream applications in:

* Bioinformatics
* Computational biology
* Comparative genomics
* Genome comparison
* Clustering
* Classification
* Anomaly detection
* Machine learning
* Artificial intelligence

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

The `GenomeComparison` object groups all comparison results into one structured representation.

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

This architecture keeps the comparison logic outside the presentation layer and allows comparison results to be reused by future matrix, ranking, export and visualization features.

---

# Real CDS-to-CDS Comparison

The primary example compares two real fluorescent-protein coding sequences:

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

# Quick Example

```python
from pathlib import Path

from src.genome import Genome


DATA_DIR = Path("data") / "fluorescent_proteins"

first_genome = Genome.from_fasta(
    DATA_DIR / "aequorea_victoria_gfp_cds.fasta"
)

second_genome = Genome.from_fasta(
    DATA_DIR / "acropora_millepora_gfp_cds.fasta"
)

first_descriptor = first_genome.descriptor(k=3)
second_descriptor = second_genome.descriptor(k=3)

comparison = first_descriptor.compare(second_descriptor)

print(f"First header: {first_genome.header}")
print(f"Second header: {second_genome.header}")

print(f"First sequence length: {first_descriptor.length} bp")
print(f"First GC content: {first_descriptor.gc_content * 100:.2f}%")
print(
    f"First Shannon entropy: "
    f"{first_descriptor.shannon_entropy:.4f} bits"
)
print(f"First k-mer diversity: {first_descriptor.kmer_diversity:.4f}")
print(f"First k-mer entropy: {first_descriptor.kmer_entropy:.4f} bits")

print(first_descriptor.to_dict())
print(first_descriptor.to_vector())
print(first_descriptor.to_normalized_vector())

print(f"Euclidean distance: {comparison.euclidean_distance:.4f}")
print(f"Cosine similarity: {comparison.cosine_similarity:.4f}")

for feature_name, difference in comparison.sorted_feature_differences():
    print(f"{feature_name}: {difference:.4f}")
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
3. Generates a mathematical descriptor for each CDS.
4. Prints basic information for the reference sequence.
5. Prints its complete mathematical descriptor.
6. Prints its raw descriptor vector.
7. Prints its normalized descriptor vector.
8. Displays the first observed k-mer frequencies.
9. Creates a structured `GenomeComparison`.
10. Calculates Euclidean distance.
11. Calculates cosine similarity.
12. Prints normalized feature differences in descending order.

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
* [x] GenomeDescriptor object
* [x] Descriptor dictionary conversion
* [x] Raw descriptor vector conversion
* [x] Normalized descriptor vector conversion

## Genome comparison

* [x] Descriptor normalization
* [x] Euclidean distance
* [x] Cosine similarity
* [x] Feature-level comparison
* [x] Explainable feature differences
* [x] GenomeComparison object
* [x] Sorted feature differences
* [x] Real CDS-to-CDS comparison
* [x] Curated fluorescent-protein dataset

## Future development

* [ ] Multiple-genome comparison
* [ ] GenomeCollection object
* [ ] Euclidean distance matrix
* [ ] Cosine similarity matrix
* [ ] Genome embeddings
* [ ] Similarity ranking
* [ ] Clustering
* [ ] Multiple FASTA record support
* [ ] Ambiguous nucleotide support
* [ ] RNA support
* [ ] Graph-based descriptors
* [ ] Spectral descriptors
* [ ] Compression-based descriptors
* [ ] Visualization tools
* [ ] Descriptor export
* [ ] Comparison export
* [ ] Embedding export

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

The current test suite contains 39 tests and covers:

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

The project currently separates sequence analysis into four responsibilities.

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

## main.py

The demonstration program:

* Loads two real fluorescent-protein CDS records
* Coordinates descriptor generation
* Presents raw and normalized vectors
* Creates a structured comparison
* Displays distance and similarity
* Displays feature-level explanations

The resulting architecture is:

```text
FASTA sequence
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
```

The next architectural step will introduce a collection abstraction for comparing more than two genomes.

---

# Why This Project?

Most genomic embedding systems rely on complex neural networks trained on massive datasets.

Those models can be powerful, but their internal representations are often difficult to interpret.

Genome Embeddings explores a different approach: describing genomic sequences through mathematical properties whose meaning remains explicit.

Rather than treating DNA as raw text, the project treats a genomic sequence as a mathematical object whose properties can be:

* Measured
* Interpreted
* Compared
* Explained
* Visualized
* Converted into vectors
* Integrated with machine-learning systems

The project emphasizes:

* Mathematical interpretability
* Reproducibility
* Explainable descriptors
* Explainable comparisons
* Modular architecture
* Automated testing
* Real biological datasets
* Comparisons between equivalent biological regions
* Extensibility
* Compatibility with future machine-learning applications

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

This CDS is used as the comparison sequence in `main.py`.

---

### `discosoma_fp583_cds.fasta`

* **Organism:** *Discosoma* species
* **Protein:** Fluorescent protein FP583
* **Accession:** AF168419.2
* **Protein ID:** AAF03369.1
* **Region:** CDS
* **Location:** 54..731
* **Sequence length:** 678 nt

This sequence is included for future multiple-genome comparison and distance-matrix analysis.

---

## Synthetic control

### `periodic_sequence.fasta`

This file contains a balanced artificial sequence constructed by repeating:

```text
ACGT
```

The current control contains 120 nucleotides.

It is retained as a methodological control for evaluating how descriptors respond to highly periodic local sequence organization.

The synthetic sequence is not used as the primary biological comparison.

---

# Dataset Design

The primary demonstration compares equivalent coding regions:

```text
Aequorea victoria GFP CDS
             versus
Acropora millepora GFP CDS
```

This is more rigorous than comparing a complete mRNA against a CDS because untranslated regions could affect:

* GC content
* Nucleotide entropy
* k-mer diversity
* k-mer entropy
* Descriptor distances

The retained complete mRNA record remains useful for future experiments exploring how genomic region selection affects mathematical descriptors.

---

# Long-Term Vision

The project is divided conceptually into three phases.

## Phase 1: Mathematical descriptors

The first phase builds reusable and interpretable measurements of genomic sequences.

Current examples include:

* Base composition
* Information entropy
* Sequence imbalance
* k-mer diversity
* k-mer entropy

Future examples may include:

* Graph-based descriptors
* Spectral descriptors
* Compression-based descriptors
* Fractal descriptors
* Additional information-theoretic descriptors

## Phase 2: Explainable genome comparison

The second phase compares descriptor vectors through mathematical metrics.

Current capabilities include:

* Descriptor normalization
* Euclidean distance
* Cosine similarity
* Feature-level difference analysis
* Structured comparison objects
* Real CDS-to-CDS comparison

Future comparison capabilities will include:

* Multiple-genome comparison
* Distance matrices
* Similarity matrices
* Similarity ranking
* Clustering
* Visualization
* Comparison export

## Phase 3: Genome embeddings

The third phase combines mathematical descriptors into reusable vector representations.

These embeddings will support:

* Genome similarity analysis
* Distance calculations
* Clustering
* Classification
* Anomaly detection
* Comparative genomics
* Machine-learning pipelines
* Artificial-intelligence applications

The goal is not to replace deep-learning approaches, but to complement them with transparent and explainable mathematical representations.

---

# License

This project is released under the MIT License.

See the [LICENSE](LICENSE) file for details.