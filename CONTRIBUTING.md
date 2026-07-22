# Contributing to Genome Embeddings

Thank you for considering a contribution to Genome Embeddings.

Contributions may include bug fixes, tests, documentation, mathematical descriptors, validation improvements, visualizations, scientific discussion and carefully scoped feature proposals.

Genome Embeddings is a research-oriented project. Contributions should preserve interpretability, reproducibility and clarity.

## Project Principles

The project follows these principles:

1. **Interpretability before complexity**
2. **Tests before or alongside implementation**
3. **Small, reviewable changes**
4. **Explicit mathematical definitions**
5. **Separation between analysis and visualization**
6. **No unsupported biological or clinical claims**
7. **Reproducible behavior across supported environments**
8. **One descriptor family or analytical concept at a time**

## Before Contributing

Before starting a substantial change:

1. search the existing issues;
2. review the roadmap in `README.md`;
3. open a feature request when the change affects the public API, mathematical model or project architecture;
4. describe the problem before proposing a large implementation.

Small corrections, tests and documentation improvements may be submitted directly through a pull request.

Security vulnerabilities must not be reported publicly. Follow `SECURITY.md` instead.

## Development Setup

Clone the repository:

```bash
git clone https://github.com/GuglielmoMarengo/genome-embeddings.git
cd genome-embeddings
```

Create and activate a virtual environment if desired:

```bash
python -m venv .venv
```

Windows PowerShell:

```powershell
.venv\Scripts\Activate.ps1
```

Linux or macOS:

```bash
source .venv/bin/activate
```

Install the dependencies:

```bash
python -m pip install -r requirements.txt
```

Run the test suite:

```bash
python -m pytest
```

Run the demonstration:

```bash
python main.py
```

## Development Workflow

The preferred workflow is:

1. write or update a test that describes the expected behavior;
2. confirm that the new test fails for the intended reason;
3. implement the smallest coherent change;
4. run the relevant tests;
5. run the complete test suite;
6. update documentation when public behavior changes;
7. review generated files and avoid committing local outputs.

Generated files under `outputs/` must not be committed.

## Code Organization

The current architecture separates responsibilities:

```text
src/genome.py
```

Contains sequence models, mathematical descriptors, comparisons, matrices and multiscale analysis.

```text
src/visualization.py
```

Contains plotting and figure-export functionality.

Core mathematical behavior should not depend on Matplotlib or graphical output.

New responsibilities should be placed in the most appropriate module rather than added automatically to `main.py`.

## Coding Style

Contributions should:

* use clear and descriptive names;
* include type annotations;
* prefer small functions with one responsibility;
* avoid unnecessary dependencies;
* preserve existing public behavior unless a breaking change is explicitly discussed;
* use explicit validation and informative error messages;
* avoid duplicated calculation logic;
* follow the formatting style already present in the project;
* keep lines and expressions readable;
* use UTF-8 encoded files.

Do not optimize prematurely. Correctness, clarity and testability take priority.

## Tests

Every behavioral change should include appropriate tests.

Tests should cover:

* expected behavior;
* invalid input;
* boundary conditions;
* ordering guarantees;
* returned data independence when relevant;
* numerical results using `pytest.approx`;
* regression cases for corrected bugs.

New visualization tests should use Matplotlib's non-interactive backend:

```python
matplotlib.use("Agg")
```

Run the complete suite before submitting a pull request:

```bash
python -m pytest
```

A pull request should not intentionally reduce test coverage without explanation.

## Mathematical and Scientific Contributions

New descriptors or metrics should include:

* a precise definition;
* the mathematical formula;
* input and output ranges;
* normalization rules;
* behavior for edge cases;
* interpretation and limitations;
* tests using simple sequences with predictable results;
* references when the method is derived from published work.

Avoid adding correlated or redundant variables without explaining why they are needed.

When introducing a new descriptor, distinguish between:

* raw sequence properties;
* normalized descriptors;
* full distribution vectors;
* distance or similarity metrics;
* biological interpretation.

## Biological Claims

Results produced by Genome Embeddings are exploratory.

Contributions must not claim that a descriptor, distance, similarity, clustering result or visualization proves:

* sequence identity;
* protein identity;
* evolutionary relatedness;
* functional equivalence;
* pathogenicity;
* clinical relevance;
* diagnostic performance.

Such claims require appropriate biological datasets, statistical analysis, external validation and, where applicable, regulatory review.

## Documentation

Update `README.md` when a contribution changes:

* installation;
* dependencies;
* public classes or methods;
* demonstration output;
* generated files;
* project structure;
* validation behavior;
* roadmap status;
* known limitations.

Documentation should explain both what a feature does and what it does not establish.

## Commit Messages

Use concise, imperative commit messages.

Examples:

```text
Add cross-scale trajectory analysis
Fix FASTA validation for empty records
Document matrix serialization
Refactor pair trajectory construction
```

Avoid vague messages such as:

```text
Update code
Changes
Fix stuff
```

## Pull Requests

A pull request should include:

* a clear title;
* the problem being addressed;
* a summary of the implementation;
* tests added or updated;
* documentation changes;
* known limitations;
* relevant issue references.

Keep pull requests focused. Unrelated changes should be submitted separately.

Before submission, verify:

```bash
python -m pytest
python main.py
```

Do not include generated files from `outputs/`, editor settings, virtual environments, caches or local credentials.

## Reporting Bugs

Use the bug report issue form.

Include:

* a minimal reproducible example;
* expected behavior;
* actual behavior;
* traceback or error output;
* Python version;
* operating system;
* relevant input characteristics.

Do not upload sensitive, proprietary, clinical or personally identifiable genomic data.

## Proposing Features

Use the feature request issue form.

Explain:

* the problem being solved;
* the proposed behavior;
* mathematical or scientific motivation;
* alternatives considered;
* expected tests;
* effect on interpretability;
* potential compatibility concerns.

## Code of Conduct

Participation in this project is governed by `CODE_OF_CONDUCT.md`.

By contributing, you agree to follow it.

## License

By submitting a contribution, you agree that your contribution may be distributed under the repository's MIT License.