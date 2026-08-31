![Fan flip graph](https://raw.githubusercontent.com/natemacfadden/regfans/main/images/fan_flip_graph.png)

# regfans
*[Nate MacFadden](https://github.com/natemacfadden), Liam McAllister Group, Cornell*

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.19406101.svg)](https://doi.org/10.5281/zenodo.19406101)

**regfans** is a Python library for studying lattice vector configurations.

**Used by** [CYTools](https://cy.tools/).

## Core Functionality

The library computes and modifies regular triangulations of vector configurations. These could also be known as "regular polyhedral fans" or "vex triangulations".

Key capabilities:
- Verify properties of vector configurations (solid, totally-cyclic) and fans (fine, regular, point-configuration-compatible),
- Construct regular triangulations via lifting,
- Compute all (regular) triangulations via flip graph traversal, and
- Perform efficient linear flips.

See [Triangulations: Structures for Algorithms and Applications](https://link.springer.com/book/10.1007/978-3-642-12971-1) by De Loera, Rambau, and Santos for a definitive reference on such topics. It is a wonderful book.

regfans has a lot of overlap with [TOPCOM](https://www.wm.uni-bayreuth.de/de/team/rambau_joerg/TOPCOM/); `regfans` will port over more computations to TOPCOM as beneficial for speed (via [triangulumancer](https://github.com/ariostas/triangulumancer)).

## Installation

Install via conda (most robust across platforms):

```
conda env create -f environment.yml
conda activate regfans
pip install -e .
```

Or via pip (see also [PyPI listing](https://pypi.org/project/regfans/)):

```
pip install regfans
```

**Note:** `regfans` now ships a C extension (the `grow4d` triangulation
backend), so installing from source needs a C compiler. Wheels carry it
prebuilt. Set `REGFANS_NATIVE=1` before `pip install -e .` to tune the build
to the local CPU, which is faster but not portable to other machines.

**Note:** Most methods require dual cone computation via [pplpy](https://pypi.org/project/pplpy/). On macOS and Linux x86_64, `pip install` pulls a prebuilt pplpy wheel automatically (no system libraries needed). On other platforms (Windows, or ARM/musl Linux on older Pythons) there is no wheel and pip falls back to a source build that needs the ppl/gmp C libraries; the conda environment supplies those reliably.

## Tests

Install the test extra and run the suite:

```
pip install -e ".[test]"
pytest
```

## Primary Interface

The main class is `VectorConfiguration`:

```python
from regfans import VectorConfiguration

pts = [[1, -2, -1, -1], [1, 1, -1, 2], [-2, 0, 0, -1],
       [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]]
vc = VectorConfiguration(pts)

# construct a regular triangulation via lifting
fan = vc.subdivide()
print(fan.is_fine(), fan.is_regular())

# compute all triangulations and the flip graph
all_fans = vc.all_triangulations()
G, fans, labels = vc.flip_graph(compute_node_labels=True)
```

See [documentation/api.md](documentation/api.md) for the full API reference and the [tutorials directory](tutorials/) for annotated examples.

## Citation

If you are using `regfans` for toric-geometric or Calabi-Yau applications, please cite [Calabi-Yau Threefolds from Vex Triangulations](https://arxiv.org/abs/2512.14817):

```bibtex
@article{MacFadden:2512.14817,
  author  = {MacFadden, Nate and Sheridan, Elijah},
  title   = {Calabi-{Y}au Threefolds from Vex Triangulations},
  doi     = {10.48550/arXiv.2512.14817},
  url     = {https://arxiv.org/abs/2512.14817},
}
```

Otherwise, please cite this repository:

```bibtex
@software{regfans,
  author  = {MacFadden, Nate},
  title   = {regfans},
  doi     = {10.5281/zenodo.19406101},
  url     = {https://github.com/natemacfadden/regfans},
  orcid   = {0000-0002-8481-3724},
}
```

This package was developed for constructing toric varieties in [Calabi-Yau Threefolds from Vex Triangulations](https://arxiv.org/abs/2512.14817), supported in part by NSF grant PHY-2309456. Toric-geometric computations are provided by [CYTools](https://github.com/LiamMcAllisterGroup/cytools), which extends regfans via a `vector_config` module.

## Organization

```
regfans/
├── src/regfans/
│   ├── vectorconfig.py   # VectorConfiguration class: triangulations, flip graph, VC properties
│   ├── fan.py            # Fan class: triangulation of a vector configuration
│   ├── circuits.py       # Circuit and Circuits classes
│   ├── grow4d.pyx        # Cython wrapper for the grow4d enumeration backend
│   ├── grow4d.h          # header-only C kernel: all fans by growth + backtracking
│   └── util.py           # shared utilities for vector configuration computations
├── tests/
│   ├── test_vectorconfig.py   # tests for VectorConfiguration
│   ├── test_fan.py            # tests for Fan
│   ├── test_circuits.py       # tests for circuits
│   └── test_util.py           # tests for utilities
├── tutorials/
│   ├── simple_lifting_ex.py   # constructing a regular triangulation via lifting
│   └── all_fans_ex.py         # computing all triangulations and the flip graph
├── documentation/
│   ├── api.md                 # full API reference
│   └── clean_api.py           # script to generate api.md
├── setup.py                   # builds the grow4d C extension
└── pyproject.toml
```

## License

[GPLv3](LICENSE). Copyright (c) 2026 Nate MacFadden.
