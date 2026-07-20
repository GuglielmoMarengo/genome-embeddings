# Genome Embeddings

[Python](https://www.python.org/)
[License](LICENSE)
[Last Commit](https://github.com/GuglielmoMarengo/genome-embeddings/commits/main)
[Status](https://github.com/GuglielmoMarengo/genome-embeddings)

**Turning genomes into mathematics.**

Genome Embeddings is an open-source Python project for representing genomic sequences through mathematically interpretable descriptors.

Rather than relying solely on machine learning, the project explores representations inspired by information theory, number theory, graph theory and statistics.

---



# Vision

Rather than relying exclusively on machine learning models, this project investigates mathematical descriptors derived from information theory, number theory, graph theory and statistics.

The long-term goal is to build interpretable genome embeddings that can support downstream bioinformatics, computational biology and AI applications.

---



# Implemented Features

Current capabilities include:

- FASTA parsing
- DNA sequence validation
- Sequence length
- GC content
- Reverse complement
- Nucleotide frequency analysis
- Shannon entropy
- k-mer extraction
- k-mer frequency analysis
- Genome descriptor objects

---



# Genome Descriptor

The current `GenomeDescriptor` object summarizes a sequence through mathematical features that are useful for downstream analysis and embedding construction.

Current descriptor fields include:

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

---



# Quick Example

```python
from src.genome import Genome

genome = Genome.from_fasta("data/gfp.fasta")
descriptor = genome.descriptor()

print(f"Length: {genome.length()} bp")
print(f"GC Content: {genome.gc_content() * 100:.2f}%")
print(f"Shannon Entropy: {genome.shannon_entropy():.4f} bits")

print(descriptor.to_dict())
print(descriptor.to_vector())
```

Example output:

```text
{
    "length": 922,
    "gc_content": 0.3623,
    "at_content": 0.6377,
    "shannon_entropy": 1.9403,
    "gc_skew": -0.0124,
    "purine_content": 0.4989,
    "pyrimidine_content": 0.5011,
    "kmer_length": 3,
    "kmer_diversity": 0.9857,
    "kmer_entropy": 4.1234
}
```

---



# Roadmap

- [x] Genome representation
- [x] FASTA parser
- [x] Sequence validation
- [x] GC content
- [x] Reverse complement
- [x] Nucleotide frequencies
- [x] Shannon entropy
- [x] k-mer extraction
- [x] k-mer frequencies
- [x] Genome descriptors
- [ ] Genome embeddings
- [ ] Genome similarity metrics
- [ ] Visualization tools
- [ ] Embedding export

---



# Genome Validation Rules


| Rule                                             | Status     |
| ------------------------------------------------ | ---------- |
| Sequence cannot be empty                         | ✅          |
| Sequence is automatically converted to uppercase | ✅          |
| Only **A**, **C**, **G** and **T** are accepted  | ✅          |
| FASTA header validation                          | ✅          |
| RNA support                                      | 🚧 Planned |


---



# Installation

```bash
git clone https://github.com/GuglielmoMarengo/genome-embeddings.git
cd genome-embeddings
pip install -r requirements.txt
```

---



# Running the Tests

```bash
python -m pytest
```

All current tests should pass successfully.

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



# Why this Project?

Most genomic embeddings rely on neural networks trained on massive datasets.

Genome Embeddings explores an alternative direction by representing genomes through mathematically interpretable descriptors.

Rather than treating a genome as raw text, this project aims to describe it as a mathematical object whose properties can be measured, compared and eventually embedded into a vector space.

The project emphasizes:

- Mathematical interpretability
- Reproducibility
- Explainable descriptors
- Extensibility
- Integration with Machine Learning

---



# Example Dataset

The repository currently includes:

**gfp.fasta**

- **Organism:** *Aequorea victoria*
- **Gene:** Green Fluorescent Protein (GFP)
- **Sequence length:** 922 bp
- **Source:** NCBI GenBank (L29345.1)

This real sequence is used to demonstrate every implemented feature of the library.

---



# Long-Term Vision

The project roadmap is intentionally divided into two phases.

The first phase focuses on building a collection of mathematically meaningful genome descriptors.

The second phase combines those descriptors into high-dimensional vector representations (genome embeddings) that can be used for:

- Genome similarity
- Clustering
- Classification
- Anomaly detection
- Machine Learning
- Artificial Intelligence applications

The goal is not to replace existing deep learning approaches, but to complement them with interpretable mathematical representations.

---



# License

This project is released under the MIT License.

See the [LICENSE](LICENSE) file for details.