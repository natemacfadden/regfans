![Fan flip graph](images/fan_flip_graph.png)

# regfans
Software for studying vector configurations defined over the lattice vectors. This includes
- constructing regular triangulations (i.e., polyhedral fans) of such vector configurations via lifting,
- constructing all (regular) triangulations via computation of flip graphs,
- verification of various properties of the vector configuration/fan, and
- efficient linear flipping.

See [Triangulations: Structures for Algorithms and Applications](https://link.springer.com/book/10.1007/978-3-642-12971-1) by De Loera, Rambau, and Santos for a definitive resource on such topics.

This package, `regfans`, was originally developed for constructing toric varieties in the work [Calabi-Yau Threefolds from Vex Triangulations](https://arxiv.org/abs/2512.14817). Said work was supported in part by NSF grant PHY-2309456. All toric-geometric computations are isolated to [CYTools](https://github.com/LiamMcAllisterGroup/cytools), which has an extension `vector_config` building off of `regfans`.

## Installation
`regfans` can be installed using either conda or pip. To install `regfans` using conda, please see/use the provided `environment.yml` file:
```
conda env create -f environment.yml
conda activate regfans
```
To install `regfans` using pip, either run (to install the most recent release; also see [PyPI listing](https://pypi.org/project/regfans/))
```
pip install regfans
```
or (to install a local version)
```
pip install .
```
N.B.: many methods in `regfans` require computation of dual cones (i.e., the generators of a cone defined via hyperplanes or vice-versa). Currently, this requires [pplpy](https://pypi.org/project/pplpy/) which cannot be automatically installed via pip.

## API

See [api.md](documentation/api.md) for full API reference.

(To update documentation, just run `pydoc-markdown; py documentation/clean_api.py`)

## Tutorials

See the [tutorials directory](tutorials/) for some commented example scripts showing how to construct a vector configuration, check properties of it, construct fans from it, and check the properties of said fans.
