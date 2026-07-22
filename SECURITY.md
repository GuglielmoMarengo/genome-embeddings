# Security Policy

Thank you for helping improve the security of Genome Embeddings.

Genome Embeddings is currently an early-stage research and software-development project. Security reports are welcome, particularly when they concern unsafe file handling, dependency vulnerabilities, unintended data exposure, path manipulation, denial-of-service risks or other behavior that could affect users of the library.

## Supported Versions

The project is under active development and has not yet reached a stable release.

Security fixes are currently applied only to the latest version available on the `main` branch.

| Version                        | Supported |
| ------------------------------ | --------- |
| Latest `main` branch           | Yes       |
| Older commits or forks         | No        |
| Unreleased local modifications | No        |

## Reporting a Vulnerability

Please do not report suspected security vulnerabilities through public GitHub issues, discussions or pull requests.

Use GitHub Private Vulnerability Reporting:

1. Open the repository's **Security** tab.
2. Select **Report a vulnerability**.
3. Provide the requested technical details.

Repository:

```text
https://github.com/GuglielmoMarengo/genome-embeddings
```

If private vulnerability reporting is temporarily unavailable, contact the project maintainer using the contact information published on the maintainer's GitHub profile.

## Information to Include

A useful vulnerability report should include, when applicable:

* a clear description of the vulnerability;
* the affected file, class, method or dependency;
* the conditions required to reproduce the issue;
* a minimal proof of concept;
* the expected security impact;
* the Python version and operating system used;
* any suggested mitigation or fix;
* whether the vulnerability has already been disclosed elsewhere.

Please avoid including sensitive personal, genomic, authentication or organizational data in the report.

Use synthetic or minimal test data whenever possible.

## What Qualifies as a Security Issue

Examples may include:

* unsafe processing of untrusted FASTA or other input files;
* path traversal or unintended file access;
* arbitrary code execution;
* command injection;
* insecure temporary-file handling;
* accidental exposure of secrets or private data;
* vulnerable dependencies;
* denial-of-service risks caused by crafted inputs;
* security-sensitive validation bypasses;
* unsafe behavior in future command-line, web or external-service integrations.

Ordinary calculation errors, unexpected scientific results, documentation issues and feature requests should normally be submitted through the public issue templates instead.

## Response Process

After receiving a report, the maintainer will aim to:

1. acknowledge the report;
2. assess whether the behavior represents a security vulnerability;
3. request additional information when necessary;
4. develop and test an appropriate correction;
5. coordinate disclosure with the reporter when appropriate;
6. publish a fix or security advisory when warranted.

Because the project is currently maintained on a best-effort basis, no guaranteed response or resolution time is promised.

## Coordinated Disclosure

Please allow reasonable time for investigation and remediation before publicly disclosing a vulnerability.

The maintainer may reject reports that:

* cannot be reproduced;
* do not affect the current project;
* concern unsupported forks or modified versions;
* describe expected behavior without a security impact;
* contain insufficient information for investigation.

## Research and Clinical Disclaimer

Genome Embeddings is not a clinical or diagnostic system.

Security reporting must not be interpreted as certification that the software is suitable for clinical, diagnostic, therapeutic or regulated use.