# regfans
Software for studying vector configurations and their triangulations. Originally developed for arXiv:2512.14817

## Installation
Computation of dual cones requires pplpy. Install via conda:
```
conda env create -f environment.yml
conda activate regfans
```
## Documentation

See [API_DOC.md](API_DOC.md) for full API reference.

To update documentation, just run `pydoc-markdown; py clean_API_DOC.py`
