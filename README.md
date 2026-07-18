# Genome Embeddings

Genome Embeddings is an open-source project that explores mathematical representations of genomic sequences.

## Vision

Instead of representing DNA using machine learning models, this project investigates mathematical descriptors derived from information theory, number theory, graph theory and statistics.

## Roadmap

- [ ] Sequence statistics

- [ ] GC Content

- [ ] Shannon Entropy

- [ ] k-mer frequencies

- [ ] Genome embeddings

- [ ] Genome comparison

## Genome Validation Rules

The `Genome` class validates all input data before creating a valid genome object.

Rules:

1. The sequence cannot be empty.
2. The sequence may be lowercase, uppercase or mixed case; it is automatically converted to uppercase.
3. The sequence must contain only valid DNA nucleotides: A, C, G and T.
4. The organism cannot be empty and must be a string.
5. The chromosome cannot be empty and must be a string.
6. Genome Embeddings v1 supports DNA only. RNA support is planned for a future release.

