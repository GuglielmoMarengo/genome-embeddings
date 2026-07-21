````markdown
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

Genome Embeddings explores a complementary direction: representing genomic sequences through transparent and reproducible mathematical descriptors.

The long-term goal is to build interpretable genome embeddings that can support downstream applications in:

- Bioinformatics
- Computational biology
- Genome comparison
- Clustering
- Classification
- Anomaly detection
- Machine learning
- Artificial intelligence

---

# Implemented Features

Current capabilities include:

- FASTA parsing
- DNA sequence validation
- Automatic uppercase normalization
- Sequence length calculation
- GC content
- AT content
- Reverse complement
- Nucleotide frequency analysis
- Shannon entropy
- GC skew
- Purine content
- Pyrimidine content
- k-mer extraction
- k-mer frequency analysis
- k-mer diversity
- k-mer entropy
- Genome descriptor generation
- Descriptor dictionary conversion
- Descriptor vector conversion

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
genome = Genome.from_fasta("data/gfp.fasta")
```

---

# Mathematical Descriptors

The `Genome` class exposes individual descriptor methods that can be calculated and tested independently.

```python
genome.length()
genome.gc_content()
genome.at_content()
genome.nucleotide_frequencies()
genome.shannon_entropy()
genome.gc_skew()
genome.purine_content()
genome.pyrimidine_content()
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

- `length`
- `gc_content`
- `at_content`
- `shannon_entropy`
- `gc_skew`
- `purine_content`
- `pyrimidine_content`
- `kmer_length`
- `kmer_diversity`
- `kmer_entropy`

The descriptor can be converted into a dictionary:

```python
descriptor.to_dict()
```

or into a numerical vector:

```python
descriptor.to_vector()
```

This vector representation is the foundation for future genome embeddings and similarity calculations.

---

# Descriptor Definitions

## Sequence Length

The total number of nucleotides in the sequence.

```text
length = number of nucleotides
```

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

---

## GC Skew

GC skew measures the relative imbalance between guanine and cytosine.

```text
GC skew = (G - C) / (G + C)
```

If a sequence contains neither guanine nor cytosine, the GC skew is defined as `0.0`.

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

---

## k-mer Diversity

k-mer diversity measures the proportion of distinct k-mers among all observed k-mers.

```text
k-mer diversity = distinct k-mers / total k-mers
```

A value close to `1.0` indicates high k-mer variety.

A lower value indicates that the sequence contains repeated k-mer patterns.

---

## k-mer Entropy

k-mer entropy applies Shannon entropy to the frequency distribution of k-mers rather than individual nucleotides.

```text
Hₖ = -Σ p(k-mer) log₂ p(k-mer)
```

This descriptor captures sequence organization that cannot be detected through nucleotide frequencies alone.

---

# Quick Example

```python
from src.genome import Genome

genome = Genome.from_fasta("data/gfp.fasta")
descriptor = genome.descriptor(k=3)

print(f"Header: {genome.header}")
print(f"Length: {descriptor.length} bp")
print(f"GC content: {descriptor.gc_content * 100:.2f}%")
print(f"AT content: {descriptor.at_content * 100:.2f}%")
print(f"Shannon entropy: {descriptor.shannon_entropy:.4f} bits")
print(f"GC skew: {descriptor.gc_skew:.4f}")
print(f"Purine content: {descriptor.purine_content * 100:.2f}%")
print(f"Pyrimidine content: {descriptor.pyrimidine_content * 100:.2f}%")
print(f"k-mer diversity: {descriptor.kmer_diversity:.4f}")
print(f"k-mer entropy: {descriptor.kmer_entropy:.4f} bits")

print(descriptor.to_dict())
print(descriptor.to_vector())
```

---

# Running the Example

Run the demonstration program from the project root:

```bash
python main.py
```

The program:

1. Loads the included GFP FASTA sequence.
2. Prints basic sequence information.
3. Generates its mathematical descriptor.
4. Prints the descriptor vector.
5. Displays the first observed k-mer frequencies.

---

# Roadmap

## Genome representation

- [x] Genome class
- [x] FASTA parser
- [x] Sequence validation
- [x] Sequence normalization
- [x] Sequence length
- [x] Reverse complement

## Mathematical descriptors

- [x] Nucleotide frequencies
- [x] GC content
- [x] AT content
- [x] Shannon entropy
- [x] GC skew
- [x] Purine content
- [x] Pyrimidine content
- [x] k-mer extraction
- [x] k-mer frequencies
- [x] k-mer diversity
- [x] k-mer entropy
- [x] GenomeDescriptor object
- [x] Descriptor dictionary conversion
- [x] Descriptor vector conversion

## Future development

- [ ] Descriptor normalization
- [ ] Genome embeddings
- [ ] Genome similarity metrics
- [ ] Euclidean distance
- [ ] Cosine similarity
- [ ] Multiple FASTA record support
- [ ] Ambiguous nucleotide support
- [ ] RNA support
- [ ] Graph-based descriptors
- [ ] Spectral descriptors
- [ ] Compression-based descriptors
- [ ] Visualization tools
- [ ] Descriptor export
- [ ] Embedding export

---

# Genome Validation Rules

| Rule | Status |
|---|---|
| Sequence cannot be empty | ✅ |
| Sequence must be a string | ✅ |
| Sequence is automatically converted to uppercase | ✅ |
| Only **A**, **C**, **G** and **T** are accepted | ✅ |
| Invalid nucleotide positions are reported | ✅ |
| FASTA file cannot be empty | ✅ |
| FASTA header must begin with `>` | ✅ |
| k must be an integer | ✅ |
| k must be greater than zero | ✅ |
| k cannot exceed the sequence length | ✅ |
| Ambiguous nucleotide support | 🚧 Planned |
| RNA support | 🚧 Planned |

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

The tests cover:

- Sequence validation
- Sequence length
- GC content
- AT content
- Reverse complement
- FASTA parsing
- k-mer extraction
- k-mer frequencies
- Nucleotide frequencies
- Shannon entropy
- GC skew
- Purine content
- Pyrimidine content
- Descriptor properties
- k-mer diversity
- k-mer entropy
- GenomeDescriptor generation

---

# Project Structure

```text
genome-embeddings/
├── data/
│   └── gfp.fasta
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

The project separates sequence analysis into three responsibilities.

## Genome

The `Genome` class:

- Stores the validated sequence
- Loads FASTA data
- Calculates individual mathematical descriptors
- Extracts and counts k-mers

## GenomeDescriptor

The `GenomeDescriptor` class:

- Stores the calculated descriptor values
- Converts descriptors into dictionaries
- Converts descriptors into numerical vectors

## main.py

The demonstration program:

- Loads the example dataset
- Coordinates the analysis
- Presents the calculated results

The resulting architecture is:

```text
Genomic sequence
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
       └── to_vector()
```

---

# Why This Project?

Most genomic embedding systems rely on complex neural networks trained on massive datasets.

Those models can be powerful, but their internal representations are often difficult to interpret.

Genome Embeddings explores a different approach: describing genomic sequences through mathematical properties whose meaning remains explicit.

Rather than treating DNA as raw text, the project treats a genomic sequence as a mathematical object whose properties can be:

- Measured
- Interpreted
- Compared
- Visualized
- Converted into vectors
- Integrated with machine-learning systems

The project emphasizes:

- Mathematical interpretability
- Reproducibility
- Explainable descriptors
- Modular architecture
- Automated testing
- Extensibility
- Compatibility with future machine-learning applications

---

# Example Dataset

The repository includes the following example sequence:

## `gfp.fasta`

- **Organism:** *Aequorea victoria*
- **Gene:** Green Fluorescent Protein
- **Abbreviation:** GFP
- **Sequence type:** mRNA, complete coding sequence
- **Sequence length:** 922 bp
- **Source:** NCBI GenBank
- **Accession:** L29345.1

This real biological sequence is used to demonstrate the implemented functionality.

---

# Long-Term Vision

The project is divided conceptually into two phases.

## Phase 1: Mathematical descriptors

The first phase builds reusable, interpretable measurements of genomic sequences.

Examples include:

- Base composition
- Information entropy
- Sequence imbalance
- k-mer diversity
- k-mer entropy
- Future graph, spectral and compression descriptors

## Phase 2: Genome embeddings

The second phase combines those measurements into numerical vector representations.

These embeddings will support:

- Genome similarity analysis
- Distance calculations
- Clustering
- Classification
- Anomaly detection
- Comparative genomics
- Machine-learning pipelines
- Artificial-intelligence applications

The goal is not to replace deep-learning approaches, but to complement them with transparent mathematical representations.

---

# License

This project is released under the MIT License.

See the [LICENSE](LICENSE) file for details.
````
