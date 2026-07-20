# Genome Embeddings

[Python](https://www.python.org/)
[License](LICENSE)
[Last Commit](https://github.com/GuglielmoMarengo/genome-embeddings/commits/main)
[Status](https://github.com/GuglielmoMarengo/genome-embeddings)

**Turning genomes into mathematics.**

Genome Embeddings is an open-source Python project for representing genomic sequences through mathematical descriptors.

Rather than relying solely on machine learning, the project explores representations inspired by information theory, number theory, graph theory and statistics.

---

## Vision

Rather than relying exclusively on machine learning models, this project investigates mathematical descriptors derived from information theory, number theory, graph theory and statistics.

The long-term goal is to build interpretable genome embeddings that can support downstream bioinformatics, computational biology and AI applications.

---



## Implemented Features

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

---



## Quick Example

```python
from src.genome import Genome

genome = Genome.from_fasta("data/gfp.fasta")

print(f"Length: {genome.length()} bp")
print(f"GC Content: {genome.gc_content():.2f}%")
print(f"Shannon Entropy: {genome.shannon_entropy():.4f} bits")
```

---



## Roadmap

- [x] Genome representation
- [x] FASTA parser
- [x] Sequence validation
- [x] GC content
- [x] Reverse complement
- [x] Nucleotide frequencies
- [x] Shannon entropy
- [x] k-mer extraction
- [x] k-mer frequencies
- [ ] Genome descriptors
- [ ] Genome embeddings
- [ ] Genome similarity metrics
- [ ] Visualization tools
- [ ] Embedding export

---



## Genome Validation Rules


| Rule                                             | Status     |
| ------------------------------------------------ | ---------- |
| Sequence cannot be empty                         | ✅          |
| Sequence is automatically converted to uppercase | ✅          |
| Only **A**, **C**, **G** and **T** are accepted  | ✅          |
| RNA support                                      | 🚧 Planned |


---



## Installation

```bash
git clone https://github.com/GuglielmoMarengo/genome-embeddings.git
cd genome-embeddings
pip install -r requirements.txt
```

---



## Running the tests

```bash
pytest
```

---



## Project Structure

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



## Why this project?

Most genomic embeddings rely on neural networks trained on massive datasets.

Genome Embeddings explores an alternative direction: representing genomes through mathematically interpretable descriptors that can be analyzed, compared and eventually integrated with machine learning models.

The project emphasizes **interpretability**, **reproducibility**, and **mathematical insight** before predictive performance.

---



## Example Dataset

**gfp.fasta**

- **Organism:** *Aequorea victoria*
- **Gene:** Green Fluorescent Protein (GFP)
- **Source:** NCBI GenBank L29345.1

This real sequence is included as an example dataset for demonstrating the library's functionality.

---



## License

This project is released under the MIT License. See the [LICENSE](LICENSE) file for details.