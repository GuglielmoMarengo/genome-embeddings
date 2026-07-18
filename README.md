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

## Roadmap

- [x] Genome representation  
- [x] Sequence validation  
- [x] GC content  
- [x] Reverse complement  
- [ ] Shannon entropy  
- [ ] k-mer frequencies  
- [ ] Genome embeddings  
- [ ] Genome similarity metrics  
- [ ] FASTA parser  
- [ ] Visualization tools  
- [ ] Embedding export

---

## Genome Validation Rules


| Rule                                             | Status     |
| ------------------------------------------------ | ---------- |
| Sequence cannot be empty                         | ✅          |
| Sequence is automatically converted to uppercase | ✅          |
| Only **A**, **C**, **G** and **T** are accepted  | ✅          |
| Organism must be a non-empty string              | ✅          |
| Chromosome must be a non-empty string            | ✅          |
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

---

## Example dataset

**gfp.fasta**

- Organism: *Aequorea victoria*

- Gene: Green Fluorescent Protein (GFP)

- Source: NCBI GenBank L29345.1

---

## License

This project is released under the MIT License. See the [LICENSE](LICENSE) file for details.