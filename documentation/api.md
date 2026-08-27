<a id="vectorconfig"></a>

---


# vectorconfig

<a id="vectorconfig.VectorConfiguration"></a>

---


## VectorConfiguration Objects

```python
class VectorConfiguration()
```

This class handles definition/operations on vector configurations. It is
analogous to CYTools' Polytope class. This object can be triangulated,
making a simplicial fan.

Constructs a `VectorConfiguration` object describing a lattice vector
configuration.

Parameters
----------
vectors : ArrayLike
    The vectors defining the VC in row format. I.e., vectors[i], the
    ith row, is treated as the ith vector in the config.
labels : Iterable[int] | None
    A list of labels for the vectors. Only integral labels are allowed.
eps : float
    Threshold for checking for non-integral vectors.
gale_basis : Iterable[int] | None
    An optional basis for the gale transform. If provided, the gale
    transform is put in a basis such that the submatrix given by these
    labels equals the identity. These labels must index a unimodular
    submatrix (determinant +/-1).

<a id="vectorconfig.VectorConfiguration.__init__"></a>

---


#### \_\_init\_\_

```python
def __init__(vectors: ArrayLike,
             labels: Iterable[int] | None = None,
             eps: float = 1e-4,
             gale_basis: Iterable[int] | None = None) -> None
```

Initializes a `VectorConfiguration` object.

Parameters
----------
vectors : ArrayLike
    The vectors defining the VC in row format. I.e., vectors[i],
    the ith row, is treated as the ith vector in the config.
labels : Iterable[int] | None, optional
    A list of integer labels for the vectors. Only integral labels
    are allowed. Defaults to None.
eps : float, optional
    Threshold for checking for non-integral vectors. Defaults to 0.0001.
gale_basis : Iterable[int] | None, optional
    An optional basis for the gale transform. If provided, the gale
    transform is put in a basis such that the submatrix given by
    these labels equals the identity. These labels must index a
    unimodular submatrix (determinant +/-1). Defaults to None.

<a id="vectorconfig.VectorConfiguration.__repr__"></a>

---


#### \_\_repr\_\_

```python
def __repr__() -> str
```

String representation of the VectorConfiguration.
(more detailed than __str__)

Returns
-------
out : str
    String representation of the object.

<a id="vectorconfig.VectorConfiguration.__str__"></a>

---


#### \_\_str\_\_

```python
def __str__() -> str
```

String description of the VectorConfiguration.
(less detailed than __repr__ but more readable)

Returns
-------
out : str
    String description of the object.

<a id="vectorconfig.VectorConfiguration.__hash__"></a>

---


#### \_\_hash\_\_

```python
def __hash__() -> int
```

Hash for the vector configuration. Defined by hashing a dictionary from
labels to vectors.

Returns
-------
out : int
    The hash.

<a id="vectorconfig.VectorConfiguration.__eq__"></a>

---


#### \_\_eq\_\_

```python
def __eq__(o: VectorConfiguration) -> bool
```

Equality checking between two VectorConfiguration objects.

Parameters
----------
o : VectorConfiguration
    The other VectorConfiguration to compare against.

Returns
-------
out : bool
    True if self==o. False if self!=o.

<a id="vectorconfig.VectorConfiguration.__ne__"></a>

---


#### \_\_ne\_\_

```python
def __ne__(o: VectorConfiguration) -> bool
```

Inequality checking between two VectorConfiguration objects.

Parameters
----------
o : VectorConfiguration
    The other VectorConfiguration to compare against.

Returns
-------
out : bool
    True if self!=o. False if self==o.

<a id="vectorconfig.VectorConfiguration.copy"></a>

---


#### copy

```python
def copy() -> VectorConfiguration
```

Copy method.

Returns
-------
out : VectorConfiguration
    A copy of the vector configuration.

<a id="vectorconfig.VectorConfiguration.labels"></a>

---


#### labels

```python
@property
def labels() -> tuple[int]
```

Returns the labels of the vectors in the VC.

Returns
-------
out : tuple[int]
    The labels of the vectors in the VC.

<a id="vectorconfig.VectorConfiguration.labels_to_inds_dict"></a>

---


#### labels\_to\_inds\_dict

```python
@property
def labels_to_inds_dict() -> dict[int, int]
```

Returns a dictionary mapping vector labels to their indices in the
vector configuration.

Returns
-------
out : dict[int, int]
    The mapping from labels to indices.

<a id="vectorconfig.VectorConfiguration.size"></a>

---


#### size

```python
@property
def size() -> int
```

Returns the number of the vectors in the VC.

Returns
-------
out : int
    The number of the vectors in the VC.

<a id="vectorconfig.VectorConfiguration.ambient_dim"></a>

---


#### ambient\_dim

```python
@property
def ambient_dim() -> int
```

Returns the ambient dimension of the VC.

Returns
-------
out : int
    The ambient dimension of the VC.

<a id="vectorconfig.VectorConfiguration.dim"></a>

---


#### dim

```python
@property
def dim() -> int
```

Returns the dimension of the VC. I.e., the dimension of the subspace
spanned by the vectors.

Returns
-------
out : int
    The dimension of the VC.

<a id="vectorconfig.VectorConfiguration.vectors"></a>

---


#### vectors

```python
def vectors(which: int | Iterable[int] | None = None) -> ArrayLike
```

Returns the vectors, optionally only those with given labels.

Parameters
----------
which : int | Iterable[int] | None, optional
    Either a single label, for which the single corresponding vector
    will be returned, or a list of labels. If not provided, then all
    vectors are returned. Defaults to None.

Returns
-------
out : ArrayLike
    The corresponding vector(s), in order specified by which.

<a id="vectorconfig.VectorConfiguration.vectors_to_labels"></a>

---


#### vectors\_to\_labels

```python
def vectors_to_labels(vectors: ArrayLike) -> int | list[int]
```

Maps the vectors to their corresponding labels

Parameters
----------
vectors : ArrayLike
    Either a single vector, for which the single corresponding label
    will be returned, or a list of vectors.

Returns
-------
out : int | list[int]
    The corresponding label(s).

<a id="vectorconfig.VectorConfiguration.labels_to_inds"></a>

---


#### labels\_to\_inds

```python
def labels_to_inds(labels: Iterable[int],
                   ambient_labels: Iterable[int] | None = None,
                   offset: int = 0) -> int | Iterable[int]
```

Maps the labels to their indices in ambient_labels, optionally with a
fixed offset.

Parameters
----------
labels : Iterable[int]
    The labels of interest.
ambient_labels : Iterable[int] | None, optional
    The ambient labels to get the indices in. If None, use all labels
    of the VectorConfiguration. Defaults to None.
offset : int, optional
    Return i+offset for i the index of a label in ambient_labels.
    Defaults to 0.

Returns
-------
out : int | Iterable[int]
    The indices of the labels.

<a id="vectorconfig.VectorConfiguration.is_solid"></a>

---


#### is\_solid

```python
def is_solid() -> bool
```

Return whether or not the VC is full-dimensional.

Returns
-------
out : bool
    True if the VC is full-dimensional. False otherwise.

<a id="vectorconfig.VectorConfiguration.is_totally_cyclic"></a>

---


#### is\_totally\_cyclic

```python
def is_totally_cyclic() -> bool
```

Return whether or not the VC is totally cyclic. That is, whether
self.conv() equals the subspace containing it (the supporting
hyperplane).

Only implemented for solid (full-dimensional) VCs, the intended use
case; raises NotImplementedError otherwise.

Returns
-------
out : bool
    True if the VC is totally cyclic. False otherwise.

<a id="vectorconfig.VectorConfiguration.is_acyclic"></a>

---


#### is\_acyclic

```python
def is_acyclic() -> bool
```

Return whether or not the VC is acyclic. That is, whether there exists
some direction psi such that
    psi.vi > 0 for all vi.

This is equivalent to defining the cone {x: vi.x >= 0} and checking if
it is full-dimensional.

Returns
-------
out : bool
    True if the VC is acyclic. False otherwise.

<a id="vectorconfig.VectorConfiguration.support"></a>

---


#### support

```python
def support() -> ArrayLike
```

Get the support of the vector configuration as a hyperplane
representation.

Returns
-------
out : ArrayLike
    The hyperplanes defining the support.

<a id="vectorconfig.VectorConfiguration.cone_contains"></a>

---


#### cone\_contains

```python
def cone_contains(cone_labels: Iterable[int],
                  vec_label: Iterable[int],
                  strict: bool = False) -> bool
```

Check if a cone, specified by cone_labels, contains the ray specified
by vec_label.

I.e., if
    H = self.cone(cone_labels).hyperplanes()
    v = self.vectors(vec_label)
    H@v >= int(strict)

Parameters
----------
cone_labels : Iterable[int]
    The labels of vectors defining the cone.
vec_label : Iterable[int]
    The label of the vector to check.
strict : bool, optional
    Whether to check if the vector is in the strict interior. Defaults
    to False.

Returns
-------
out : bool
    Whether the associated cone contains the vector.

<a id="vectorconfig.VectorConfiguration.gale"></a>

---


#### gale

```python
def gale(set_basis: bool = False) -> ArrayLike
```

Compute the gale transform of the config.

I.e., a basis of the null-space of the vectors.

Parameters
----------
set_basis : bool, optional
    Whether to set a particular basis of the Gale transform. Defaults
    to False.

Returns
-------
out : ArrayLike
    The gale transform.

<a id="vectorconfig.VectorConfiguration.project"></a>

---


#### project

```python
def project(vec: ArrayLike) -> ArrayLike
```

Project down a vector from height-space to chamber-space.

Parameters
----------
vec : ArrayLike
    The height-space vector.

Returns
-------
out : ArrayLike
    The chamber-space vector.

<a id="vectorconfig.VectorConfiguration.jorp"></a>

---


#### jorp

```python
def jorp(vec: ArrayLike) -> ArrayLike
```

Undo a projection from height-space to chamber-space.

I.e., map from chamber-space to height-space

Parameters
----------
vec : ArrayLike
    The chamber-space vector.

Returns
-------
out : ArrayLike
    The height-space vector.

<a id="vectorconfig.VectorConfiguration.triangulate"></a>

---


#### triangulate

```python
def triangulate(heights: ArrayLike | None = None,
                cells: ArrayLike | None = None,
                tol: float = 1e-14,
                backend: str | None = None,
                make_fine: bool = None,
                check_heights: bool = True,
                cure_heights: bool = True,
                verbosity: int = 0) -> Fan
```

Subdivide the vector configuration either by specified cells/simplices
or by heights.

Parameters
----------
heights : ArrayLike | None, optional
    The heights to lift the vectors by. Defaults to None.
cells : ArrayLike | None, optional
    The cells to use in the triangulation. Defaults to None.
tol : float, optional
    Numerical tolerance used for curing negative heights. Defaults to
    1e-14.
backend : str | None, optional
    The lifting backend. Currently allowed to be "cgal" or "ppl".
    Defaults to None.
make_fine : bool, optional
    Convert the triangulation to a fine triangulation, if not already
    fine. Defaults to None.
check_heights : bool, optional
    Whether to check that the heights land in the secondary cone of the
    output triangulation. Defaults to True.
cure_heights : bool, optional
    If the heights do not land in the secondary cone, try to cure them
    by linear flipping towards the heights. Defaults to True.
verbosity : int, optional
    The verbosity level. Higher is more verbose. Defaults to 0.

Returns
-------
out : Fan
    The resultant subdivision.

<a id="vectorconfig.VectorConfiguration.all_triangulations"></a>

---


#### all\_triangulations

```python
def all_triangulations(only_fine: bool = False,
                       only_regular: bool = True,
                       verbosity: int = 0) -> list[Fan]
```

Generate all triangulations of this vector configuration via taking
flips from some regular triangulation.

NOTE: In theory, this might miss an irregular triangulation that is
disconnected from the regular triangulations.

Such irregular triangulations exist (see "A Point Set Whose Space of
Triangulations is Disconnected" by Santos) but are likely exceedingly
rare. E.g., it is unknown whether such cases can occur in 4D.

Could instead compute this via computing incidence vectors but that'd
be *much* slower. Roughly, this would be to
    1) compute all possible simplices
    2) if there are N possible simplices, construct an N-dim space
    3) define all 0/1-vectors. For each 0/1-vector, check if it defines
       a valid triangulation. If so, save it
The incidence vector strategy is analogous to rejection sampling and
will be much slower than the flip-based method, but it would see *all*
triangulations.

Parameters
----------
only_fine : bool, optional
    Whether to restrict to fine triangulations. Defaults to False.
only_regular : bool, optional
    Whether to restrict to regular triangulations. Defaults to True.
verbosity : int, optional
    The verbosity level. Higher is more verbose. Defaults to 0.

Returns
-------
out : list[Fan]
    A list of Fan objects, one for each triangulation of the VC.

<a id="vectorconfig.VectorConfiguration.random_triangulations_fast"></a>

---


#### random\_triangulations\_fast

```python
def random_triangulations_fast(
        method: str = "delaunay",
        h0: ArrayLike | None = None,
        sigma: float = 0.1,
        N: int | None = None,
        as_list: bool = False,
        attempts_per_triang: int = 1000,
        backend: str | None = None,
        seed: int = 0,
        verbosity: int = 0) -> Generator[Fan] | list[Fan]
```

Generate random regular triangulations by picking random heights.

Parameters
----------
method : str, optional
    Either "delaunay" or "isotropic". The former picks heights around
    some input height (e.g., the Delaunay heights). The latter picks
    heights isotropically. Defaults to "delaunay".
h0 : ArrayLike | None, optional
    The reference heights, for Delaunay method. Defaults to None.
sigma : float, optional
    How big of a distribution to study around h0. Defaults to 0.1.
N : int | None, optional
    The number of triangulations to generate. If as_list, then code
    will keep track of all triangulations, retrying at most
    attempts_per_triang tries to get a new triangulation until the list
    has N triangs. O/w, then the first N height vectors are used
    (regardless of duplicates). Defaults to None.
as_list : bool, optional
    Whether to return the triangulations as a list, or as a generator.
    Defaults to False.
attempts_per_triang : int, optional
    Quit if we can't generate a new triangulation after this many
    tries. Defaults to 1000.
backend : str | None, optional
    The lifting backend. See `VectorConfiguration.triangulate`. Defaults
    to None.
seed : int, optional
    A random number seed. Defaults to 0.
verbosity : int, optional
    The verbosity level. Higher is more verbose. Defaults to 0.

Returns
-------
out : Generator[Fan] | list[Fan]
    The random triangulations.

<a id="vectorconfig.VectorConfiguration.circuit"></a>

---


#### circuit

```python
def circuit(labels: Iterable[int],
            lmbda: Iterable | None = None,
            set_non_dependencies: bool = True,
            save_circuits: bool = True,
            enforce_positive: int | None = None) -> Circuit
```

Format/compute the circuit corresponding to the specified labels.

Parameters
----------
labels : Iterable[int]
    Labels indicating the vectors in the circuit.
lmbda : Iterable | None, optional
    Vector demonstrating the dependence. Defaults to None.
set_non_dependencies : bool, optional
    Whether to update our list of non-circuits. Defaults to True.
save_circuits : bool, optional
    Whether to save circuits... best to keep True for most
    circumstances. Defaults to True.
enforce_positive : int | None, optional
    A label to enforce is in Zpos, overriding the default orientation
    by support size. Returns None if that label has a vanishing
    coefficient. What is cached stays in the default orientation, so
    this may be varied between calls. Defaults to None.

Returns
-------
out : Circuit
    Circuit object containing
        - the support of the circuit as property 'Z',
        - the signed circuit as property 'Zpos' and 'Zneg',
        - the dependency as property 'lmbda', and
        - the signature as property 'signature'.

<a id="vectorconfig.VectorConfiguration.circuits"></a>

---


#### circuits

```python
def circuits(verbosity: int = 0) -> list[Circuit]
```

Compute all possible circuits of this vector configuration.

Parameters
----------
verbosity : int, optional
    The verbosity level. Higher is more verbose. Defaults to 0.

Returns
-------
out : list[Circuit]
    A list of Circuit objects.

<a id="vectorconfig.VectorConfiguration.flip_graph"></a>

---


#### flip\_graph

```python
def flip_graph(max_flips: int | None = None,
               only_fine: bool = False,
               only_regular: bool = True,
               only_pc_triang: bool = False,
               compute_node_labels: bool = False,
               verbosity: int = 0) -> (nx.Graph, list[Fan], list[dict])
```

Compute the flip graph. Wrapper for flip_subgraph.

Parameters
----------
max_flips : int | None, optional
    The maximum number of flips to take from the seed. If none is
    provided, then the entire flip graph is calculated. Defaults to
    None.
only_fine : bool, optional
    Whether to only compute fine triangulations. Defaults to False.
only_regular : bool, optional
    Whether to only compute regular triangulations. Note, we never will
    see irregular triangulations that are not connected to regular
    ones. Defaults to True.
only_pc_triang : bool, optional
    Whether to only compute triangulations that also correspond to star
    triangulations of the underlying point config. Defaults to False.
compute_node_labels : bool, optional
    Whether to check whether each node is fine, regular, and a PC
    triangulation. Defaults to False.
verbosity : int, optional
    The verbosity level. Higher is more verbose. Defaults to 0.

Returns
-------
out : (nx.Graph, list[Fan], list[dict])
    - The flip graph as a networkx.Graph object.
    - A list of the triangulations
    - A list of the labels for each triangulation (labels are a
      dictionary from the property to a bool)

<a id="vectorconfig.VectorConfiguration.secondary_fan"></a>

---


#### secondary\_fan

```python
def secondary_fan(only_fine: bool = False,
                  formal_fan: bool = False,
                  verbosity: int = 0)
```

Compute the secondary fan of the vector configuration.

Parameters
----------
only_fine : bool, optional
    Restrict to fine triangulations. Defaults to False.
formal_fan : bool, optional
    Save as a formal Fan object. Defaults to False.
verbosity : int, optional
    The verbosity level. Higher is more verbose. Defaults to 0.

Returns
-------
out : tuple[list, list[Fan]]
    The secondary fan triangulations (secondary cones and Fan list).

<a id="vectorconfig.VectorConfiguration.central_fan"></a>

---


#### central\_fan

```python
def central_fan() -> Fan
```

Generate the central fan of the vector configuration. Can be defined
as lifting each vector by a height of 1.

Returns
-------
out : Fan
    The central fan.

<a id="util"></a>

---


# util

<a id="util.gcd"></a>

---


#### gcd

```python
def gcd(vals: list[float], max_denom: float = 10**6) -> float
```

Computes the 'GCD' of a collection of floating point numbers.
This is the smallest number, g, such that g*values is integral.

This is computed by
    1) converting `values` to be rational [n0/d0, n1/d1, ...],
    2) computing the LCM, lcm, of [d0, d1, ...],
    3) computing the GCD, g', of [lcm*n0/d0, lcm*n1/d1, ...], and then
    4) returning g=g'/lcm.

Parameters
----------
vals : list[float]
    The numbers to compute the GCD of.
max_denom : float, optional
    Assert |di| <= max_denom. Defaults to 10 ** 6.

Returns
-------
gcd : float
    The minimum number g' such that g'*vals is integral.

<a id="util.primitive"></a>

---


#### primitive

```python
def primitive(vec: list[float], max_denom=10**6)
```

Computes the primitive vector associated to the input ray {c*vec: c>=0}.
Very similar to the gcd function.

This is equivalent to
    vec/gcd(vec)
but just uses a rational representation.

Parameters
----------
vec : list[float]
    A vector defining the ray {c*vec: c>=0}
max_denom : float, optional
    Assert |di| <= max_denom. Defaults to 10 ** 6.

Returns
-------
vec : list[int]
    The primitive vector along the ray.

<a id="util.lerp"></a>

---


#### lerp

```python
def lerp(p0: ArrayLike, p1: ArrayLike, t: float) -> ArrayLike
```

Computes the point specified by t along the line passing through p0 and p1.

Particular values:
    -) t=0   -> p0
    -) t=0.5 -> (p0+p1)/2
    -) t=1   -> p1

Parameters
----------
p0 : ArrayLike
    One point.
p1 : ArrayLike
    The other point.
t : float
    Parameter specifying where along the line Conv({p0, p1}) to return.

Returns
-------
pt : ArrayLike
    The point p0 + t*(p1-p0).

<a id="util.first_hit"></a>

---


#### first\_hit

```python
def first_hit(p0: ArrayLike,
              p1: ArrayLike,
              H: ArrayLike,
              verbosity: int = 0) -> (int, float)
```

Given a point p0 in a convex cone {x: Hx>=0}, find the first hyperplane hit
along the direction (p1-p0). I.e., the first intersection of the ray
{p0+t*(p1-p0): t>=0} with the cones bounding hyperplanes.

Allow violated hyperplanes (i.e., n such that n.p0 < 0) but ignore them.

Parameters
----------
p0 : ArrayLike
    One point.
p1 : ArrayLike
    The other point.
H : ArrayLike
    An array of hyperplane normals (as rows).
verbosity : int, optional
    The verbosity level. Higher is more verbose. Defaults to 0.

Returns
-------
out : (int, float)
    The index, i, of the first-hit hyperplane.
    The distance, t, such that dot(H[i], lerp(p0,p1,t)) = 0.

<a id="util.dual_cone"></a>

---


#### dual\_cone

```python
def dual_cone(data: ArrayLike) -> ArrayLike
```

Compute the data of the cone dual to the input 'primal' cone.

This can be thought of in a couple of equivalent ways, summarized in the
following table. E.g., if rays of the primal are input, then the
hyperplanes of the primal are output (or, equivalently, the rays of the
dual).

INPUT       | PRIMAL OUTPUT | DUAL OUTPUT
-----------------------------------------
rays        | hyperplanes   | rays
hyperplanes | rays          | hyperplanes


For simplicity in the following discussion, take the convention that one
maps hyperplanes of the primal to rays of the primal.

Parameters
----------
data : ArrayLike
    An array whose rows represent rays of the primal cone. (see table)

Returns
-------
out : ArrayLike
    An array whose rows represent hyperplanes of the primal cone.
    (see table)

<a id="util.cone_dim"></a>

---


#### cone\_dim

```python
def cone_dim(*, R: ArrayLike = None, H: ArrayLike = None) -> int
```

Return the dimension of the cone.

The cone is either specified via rays,
    {R.T @ lambda: lambda>=0},
or via hyperplanes,
    {x: H @ x>=0}.

Parameters
----------
R : ArrayLike, optional
    The rays of the cone as rows. Defaults to None.
H : ArrayLike, optional
    The hyperplanes defining the cone. Defaults to None.

Returns
-------
dim : int
    The dimension of the cone

<a id="util.is_full_rank"></a>

---


#### is\_full\_rank

```python
def is_full_rank(R: ArrayLike) -> bool
```

Return whether the rows of R are linearly independent (exact).

Assumes integer entries (consistent with the rest of the library).

Parameters
----------
R : ArrayLike
    The rays as rows.

Returns
-------
out : bool
    True iff the rows of R are linearly independent.

<a id="util.is_solid"></a>

---


#### is\_solid

```python
def is_solid(*, R: ArrayLike = None, H: ArrayLike = None) -> bool
```

Return whether the cone is full-dimensional.

The cone is either specified via rays,
    {R.T @ lambda: lambda>=0},
or via hyperplanes,
    {x: H @ x>=0}.

Parameters
----------
R : ArrayLike, optional
    The rays of the cone as rows. Defaults to None.
H : ArrayLike, optional
    The hyperplanes defining the cone. Defaults to None.

Returns
-------
out : bool
    Whether the cone is full-dimensional.

<a id="util.contains"></a>

---


#### contains

```python
def contains(*,
             p: ArrayLike,
             R: ArrayLike = None,
             H: ArrayLike = None) -> bool
```

Return if the point p is contained in the cone.

The cone is either specified via rays,
    {R.T @ lambda: lambda>=0},
or via hyperplanes,
    {x: H @ x>=0}.

Parameters
----------
p : ArrayLike
    The query point.
R : ArrayLike, optional
    The rays of the cone as rows. Defaults to None.
H : ArrayLike, optional
    The hyperplanes defining the cone. Defaults to None.

Returns
-------
out : bool
    Whether p is contained in the cone.

<a id="util.find_interior_point"></a>

---


#### find\_interior\_point

```python
def find_interior_point(*,
                        R: ArrayLike = None,
                        H: ArrayLike = None,
                        stretching: float = 1,
                        nonneg: bool = False,
                        verbosity: int = 0) -> ArrayLike | None
```

Returns a point p in the strict interior of a cone. The cone can be
specified either via its rays or its generators.

If no point p exists, return `None`.

Modified from CYTools' `Cone.find_interior_point`.

Parameters
----------
R : ArrayLike, optional
    Generators defining the cone. Defaults to None.
H : ArrayLike, optional
    Hyperplanes defining the cone. Defaults to None.
stretching : float, optional
    How far p must be from any hyperplane. Defaults to 1.
nonneg : bool, optional
    Whether to restrict to non-negative vectors. Defaults to False.
verbosity : int, optional
    The verbosity level. Defaults to 0.

Returns
-------
p : ArrayLike | None
    A point p in the strict interior.

<a id="fan"></a>

---


# fan

<a id="fan.Fan"></a>

---


## Fan Objects

```python
class Fan()
```

This class handles definition/operations on fans. It is analogous to
CYTools' Triangulation class.

Constructs a `Fan` object describing a triangulation of a lattice vector
configuration.

This class is *not* intended to be called directly. Instead, it is meant to
be called through VectorConfiguration.triangulate.

Parameters
----------
vc : VectorConfiguration
    The ambient vector configuration that this fan is over.
cones : list[list[int]]
    The cones defining the fan. Each cone is a collection of integer labels.
heights : list[float] | None
    The heights defining the fan, if it is regular. Can be computed later.

<a id="fan.Fan.__init__"></a>

---


#### \_\_init\_\_

```python
def __init__(vc: VectorConfiguration,
             cones: list[list[int]],
             heights: list[float] | None = None) -> None
```

Initializes a `Fan` object.

Parameters
----------
vc : VectorConfiguration
    The ambient vector configuration that this fan is over.
cones : list[list[int]]
    The cones defining the fan. Each cone is a collection of
    integer labels.
heights : list[float] | None, optional
    The heights defining the fan, if it is regular. Can be
    computed later. Defaults to None.

<a id="fan.Fan.__repr__"></a>

---


#### \_\_repr\_\_

```python
def __repr__() -> str
```

String representation of the Fan.
(more detailed than __str__)

Returns
-------
out : str
    String representation of the object.

<a id="fan.Fan.__str__"></a>

---


#### \_\_str\_\_

```python
def __str__() -> str
```

String description of the Fan.
(less detailed than __repr__ but more readable)

Returns
-------
out : str
    String description of the object.

<a id="fan.Fan.__hash__"></a>

---


#### \_\_hash\_\_

```python
def __hash__() -> int
```

Hash for the fan. Defined by hashing vector configuration and the cones.

Returns
-------
out : int
    The hash.

<a id="fan.Fan.__eq__"></a>

---


#### \_\_eq\_\_

```python
def __eq__(o: Fan) -> bool
```

Equality checking between two Fan objects.

Parameters
----------
o : Fan
    The other Fan to compare against.

Returns
-------
out : bool
    True if self==o. False if self!=o.

<a id="fan.Fan.__ne__"></a>

---


#### \_\_ne\_\_

```python
def __ne__(o: Fan) -> bool
```

Inequality checking between two Fan objects.

Parameters
----------
o : Fan
    The other Fan to compare against.

Returns
-------
out : bool
    True if self!=o. False if self==o.

<a id="fan.Fan.vector_config"></a>

---


#### vector\_config

```python
@property
def vector_config() -> VectorConfiguration
```

Returns the associated vector configuration.

Returns
-------
vc : VectorConfiguration
    The associated vector configuration.

<a id="fan.Fan.labels"></a>

---


#### labels

```python
@property
def labels() -> tuple[int]
```

Returns the labels of the vectors in the VC.

Returns
-------
out : tuple[int]
    The labels of the vectors in the VC.

<a id="fan.Fan.used_labels"></a>

---


#### used\_labels

```python
@property
def used_labels() -> tuple[int]
```

Returns the labels of the vectors in the VC used by cones in the Fan.

Returns
-------
out : tuple[int]
    The labels of the vectors in the VC used by cones in the Fan.

<a id="fan.Fan.labels_to_cones"></a>

---


#### labels\_to\_cones

```python
@property
def labels_to_cones() -> dict[int, set[tuple[int]]]
```

Returns a map from vector labels to the cones the vector appears in.

Returns
-------
out : dict[int, set[tuple[int]]]
    A map from vector label to a set of cones (as tuples of labels) that
    the vector appears in.

<a id="fan.Fan.ambient_dim"></a>

---


#### ambient\_dim

```python
@property
def ambient_dim() -> int
```

Returns the ambient dimension of the VC.

Returns
-------
out : int
    The ambient dimension of the VC.

<a id="fan.Fan.dim"></a>

---


#### dim

```python
@property
def dim() -> int
```

Returns the dimension of the VC. I.e., the dimension of the subspace
spanned by the vectors.

Returns
-------
out : int
    The dimension of the VC.

<a id="fan.Fan.vectors"></a>

---


#### vectors

```python
def vectors(which: int | Iterable[int] | None = None,
            lifted: bool = False) -> ArrayLike
```

Returns the vectors, optionally only those with given labels. Also,
optionally, give the vectors lifted by the heights (if the Fan is
regular).

Parameters
----------
which : int | Iterable[int] | None, optional
    Either a single label, for which the single corresponding
    vector will be returned, or a list of labels. If not
    provided, then all vectors are returned. Defaults to None.
lifted : bool, optional
    Whether to give the lifted vectors. Defaults to False.

Returns
-------
out : ArrayLike
    The corresponding vector(s), in order specified by which.

<a id="fan.Fan.cones"></a>

---


#### cones

```python
def cones(dim: int = None,
          as_rays: bool = False,
          as_hyps: bool = False,
          as_inds: bool = False,
          ind_offset: int = 0) -> tuple[tuple[int]] | list[ArrayLike]
```

Returns the cones in the fan in a variety of formats. They are:
    - (default) as a tuple of labels
    - (as_rays=True) as an array whose rows are the generators
    - (as_hyps=True) as an array whose rows are hyperplane normals
    - (as_inds=True) as a tuple of indices
Optionally, allow an offset to the indices.

By default the maximal cones are returned. If `dim` is set, then return
the `dim`-dimensional cones (faces of maximal ones). Only implemented
for simplicial fans currently

Parameters
----------
dim : int, optional
    If set, return the `dim`-dimensional sub-cones Only implemented
    for simplicial fans. Defaults to None.
as_rays : bool, optional
    Whether to return the cones as their generators. Defaults to False.
as_hyps : bool, optional
    Whether to return the cones as their hyperplanes. Defaults to False.
as_inds : bool, optional
    Whether to return the cones as indices (not labels).
    Defaults to False.
ind_offset : int, optional
    An additive offset for the indices. Defaults to 0.

Returns
-------
out : tuple[tuple[int]] | list[ArrayLike]
    The cones, specified according to the input flags.

<a id="fan.Fan.facets"></a>

---


#### facets

```python
def facets() -> dict[tuple[int], list[tuple[int]]]
```

Returns the facets of the cones. Save them as a dictionary from facet
labels to a list of containing cones, stored by their labels.

Only implemented for simplicial fans (triangulations), the intended
use case; raises NotImplementedError otherwise.

Returns
-------
out : dict[tuple[int], list[tuple[int]]]
    A dictionary from facet labels to a list of containing cones.

<a id="fan.Fan.is_valid"></a>

---


#### is\_valid

```python
def is_valid(verbosity: int = 0) -> bool
```

Return whether or not the cones define a valid polyhedral fan.

Follows cor. 4.5.13 of "Triangulations" by De Loera, Rambau, Santos.
Implements the MaxMP, MaxAdjHP, and (single-point) IPP conditions; the
MaxAdjLP condition is a no-op for triangulations and is unimplemented
for subdivisions.

This can only check full-dimensional, integral triangulations
(simplicial fans).

Parameters
----------
verbosity : int, optional
    The verbosity level. Higher is more verbose. Defaults to 0.

Returns
-------
out : bool
    True if the cones define a valid fan. False otherwise.

<a id="fan.Fan.respects_ptconfig"></a>

---


#### respects\_ptconfig

```python
def respects_ptconfig(via_circuits=False) -> bool
```

Return whether or not the fan also defines a (star) subdivision of the
original underlying point configuration.

Only implemented for regular fans, the intended use case; raises
NotImplementedError otherwise.

Parameters
----------
via_circuits : bool, optional
    This method uses secondary cones. Allow construction of the
    secondary cone via circuits. That's unsafe if the fan is
    potentially irregular. Defaults to False.

Returns
-------
out : bool
    True if the fan defines a subdivision of the point
    configuration. False otherwise.

<a id="fan.Fan.is_triangulation"></a>

---


#### is\_triangulation

```python
def is_triangulation() -> bool
```

Return whether or not the fan is a triangulation (not a subdivision).

Only checks that every cone is a simplex; it does not separately
re-verify fan validity (see is_valid).

Returns
-------
out : bool
    True if the fan is a triangulation. False otherwise.

<a id="fan.Fan.is_fine"></a>

---


#### is\_fine

```python
def is_fine() -> bool
```

Return whether or not the fan is fine.

Returns
-------
out : bool
    True if the fan is fine. False otherwise.

<a id="fan.Fan.is_regular"></a>

---


#### is\_regular

```python
def is_regular() -> bool
```

Return whether or not the fan is regular.

Returns
-------
out : bool
    True if the fan is regular. False otherwise.

<a id="fan.Fan.heights"></a>

---


#### heights

```python
def heights() -> list[float] | None
```

Return some heights defining the cone, if it is regular. Else, return
None.

Returns
-------
out : list[float] | None
    The heights defining the fan, if it is regular.

<a id="fan.Fan.contains"></a>

---


#### contains

```python
def contains(c: Iterable[int] | Iterable[Iterable[int]]) -> bool
```

Check if any cone (specified by its labels) is contained in the fan.
The cone need not be solid. Can also be called for a collection of
cones, in which case the check is if all cones are contained in the fan.

Parameters
----------
c : Iterable[int] | Iterable[Iterable[int]]
    The cone(s). Either a single collection of cone, specified
    by an iterable of labels, or a collection of cones, each
    specified by an iterable of labels.

Returns
-------
out : bool
    Whether (all) cone(s) are contained in the fan.

<a id="fan.Fan.circuit"></a>

---


#### circuit

```python
def circuit(labels: Iterable[int] | None = None,
            enforce_positive: int | None = None,
            lmbda: Iterable[float] | None = None,
            check_containment: bool = True,
            save_circuits_in_vc: bool = False,
            verbosity: int = 0) -> Circuit
```

Format/compute the circuit corresponding to the specified labels.

Parameters
----------
labels : Iterable[int] | None, optional
    Labels indicating the vectors in the circuit. Defaults to None.
enforce_positive : int | None, optional
    A label to enforce is in Zpos. Defaults to None.
lmbda : Iterable[float] | None, optional
    A dependency demonstrating the circuit. Defaults to None.
check_containment : bool, optional
    Whether to check that this fan contains every cone in the
    positive triangulation, Tpos. Defaults to True.
save_circuits_in_vc : bool, optional
    Whether to save circuits... best to keep True for most
    circumstances. Defaults to False.
verbosity : int, optional
    The verbosity level. Higher is more verbose. Defaults to 0.

Returns
-------
out : Circuit
    Circuit object containing
        - the support of the circuit as property 'Z',
        - the signed circuit as property 'Zpos' and 'Zneg',
        - the dependency as property 'lmbda', and
        - the signature as property 'signature'.

<a id="fan.Fan.circuits"></a>

---


#### circuits

```python
def circuits(facets: dict[Iterable[int], Iterable[Iterable[int]]]
             | None = None,
             verbosity: int = 0) -> list[Circuit]
```

Compute all circuits associated to this fan (i.e., those 'embedded' in
this fan). All will be oriented such that the positive triangulation
(i.e., Tpos/T_+) is embedded in the fan. This enables us to directly
interpret lambda as the normal in the secondary cone.

Parameters
----------
facets : dict[Iterable[int], Iterable[Iterable[int]]] | None, optional
    The facets of the fan (not just the VC...). I.e., codim-1
    cones. Defaults to None.
verbosity : int, optional
    The verbosity level. Higher is more verbose. Defaults to 0.

Returns
-------
out : list[Circuit]
    A list of Circuit objects for all circuits embedded in the fan.

<a id="fan.Fan.star"></a>

---


#### star

```python
def star(cell: Iterable[int], old_way: bool = False) -> Iterable[tuple[int]]
```

Compute the star of some cell. This is the subcomplex of all cones
containing the cell (and their faces)

Parameters
----------
cell : Iterable[int]
    The cell of interest.
old_way : bool, optional
    Whether to do the computation in an old/slow manner.
    Defaults to False.

Returns
-------
out : Iterable[tuple[int]]
    A list of all solid cones (as tuples of ints) containing the cell.

<a id="fan.Fan.link"></a>

---


#### link

```python
def link(cell: Iterable[int]) -> list[tuple[int]]
```

Compute the link of some cell. This is the subcomplex of all cones in
the star that don't intersect the cell.

Parameters
----------
cell : Iterable[int]
    The cell of interest.

Returns
-------
out : list[tuple[int]]
    The link.

<a id="fan.Fan.embed"></a>

---


#### embed

```python
def embed(cell: Iterable[int],
          link: Iterable[Iterable[int]]) -> list[tuple[int]]
```

Embed some cell into the Fan by combining it with each cell in the link.

Parameters
----------
cell : Iterable[int]
    The cell of interest.
link : Iterable[Iterable[int]]
    The link of said cell.

Returns
-------
out : list[tuple[int]]
    A list of solid cones representing the embedding of the cell
    into the Fan via the link.

<a id="fan.Fan.flip"></a>

---


#### flip

```python
def flip(circ: Circuit,
         formal: bool = True,
         verbosity: int = 0) -> Fan | tuple[tuple[int]]
```

Make a flip across a circuit.

Parameters
----------
circ : Circuit
    The circuit to flip through.
formal : bool, optional
    Whether to return a formal Fan (otherwise, just a tuple of
    cones). Defaults to True.
verbosity : int, optional
    The verbosity level. Higher is more verbose. Defaults to 0.

Returns
-------
out : Fan | tuple[tuple[int]]
    The flipped Fan.

<a id="fan.Fan.flip_linear"></a>

---


#### flip\_linear

```python
def flip_linear(
    h_target: Iterable[float] | None = None,
    direction: Iterable[float] | None = None,
    h_init: Iterable[float] | None = None,
    max_N_flips: int | None = None,
    stop_at_deletion: bool = True,
    stop_at_pct: bool = False,
    check_regularity: bool = True,
    record_fans: bool = False,
    record_circs: bool = False,
    hook_init: Callable | None = None,
    hook_halt: Callable | None = None,
    hook_flip: Callable | None = None,
    eps: float = 1e-8,
    verbosity: int = 0
) -> list[int | Exception, ArrayLike, Fan, ArrayLike, int]
```

Compute all flips along the linear height homotopy
    t*h_target + (1-t)*h_init
for t=0 increasing to t=1.

Allow early stops of this homotopy at a certain number `max_N_flips` of
flips. Also allow early stopping upon the following conditions
    - (default True) reaching a deletion flip or
    - (default False) hitting a fan that respects the point config.

Parameters
----------
h_target : Iterable[float] | None, optional
    The target heights. Defaults to None.
direction : Iterable[float] | None, optional
    The direction to travel. Defaults to None.
h_init : Iterable[float] | None, optional
    The initial heights (regular triangulations don't have
    unique heights, even up to scaling... any h in the secondary
    cone is valid. If this is left unset, then arbitrary valid
    heights are chosen) (early stopping). Defaults to None.
max_N_flips : int | None, optional
    The maximum number of flips allowed. Defaults to None.
stop_at_deletion : bool, optional
    Whether to early-terminate the homotopy at any deletion flip
    seen. Defaults to True.
stop_at_pct : bool, optional
    Whether to early-terminate the homotopy at any fan that
    respects the point configuration. (sanity checks). Defaults
    to False.
check_regularity : bool, optional
    This method is inherently regular (it uses heights...). We
    can check the regularty of the initial fan. (record
    keeping). Defaults to True.
record_fans : bool, optional
    Whether to record the fans seen along the homotopy. Defaults
    to False.
record_circs : bool, optional
    Whether to record the circuits flipped along the homotopy.
    (numerical parameters). Defaults to False.
eps : float, optional
    A small number for an allowed violation of heights landing
    outside the secondary fan (in case the heights 'truly'
    landed on a wall of the secondary fan). Such violations are
    naturally resolved by pulling heights back into the secondary
    fan. (diagnostics). Defaults to 1e-08.
verbosity : int, optional
    The verbosity level. Higher is more verbose. Defaults to 0.

Returns
-------
out : list[int | Exception, ArrayLike, Fan, ArrayLike, int]
    - The status of the homotopy. Either 1 (if successful) or
    an Exception.
    - The current heights at the end of the homotopy. Not
    always h_target.
    - The associated fan at the end of the homotopy.
    - The hyperplanes of the secondary cone at the end of the homotopy.
    - The number of flips taken.

<a id="fan.Fan.neighbors"></a>

---


#### neighbors

```python
def neighbors(
        only_fine: bool = False,
        formal: bool = True,
        verbosity: int = 0
) -> tuple[list[Fan | tuple[tuple[int]]], list[Circuit]]
```

Compute the neighboring fans (those reachable by a single flip).

Allow restrictions to only fine fans.

Parameters
----------
only_fine : bool, optional
    Whether to only compute/return fine neighbors. Defaults to False.
formal : bool, optional
    Whether to return the neighbors as formal fans (if False,
    just return cones). Defaults to True.
verbosity : int, optional
    The verbosity level. Higher is more verbose. Defaults to 0.

Returns
-------
out : tuple[list[Fan | tuple[tuple[int]]], list[Circuit]]
    - The neighbors, either as formal Fan objects or as collections of
      cones (each cone a collection of labels)
    - The circuits flipped to get the corresponding neighbors.

<a id="fan.Fan.secondary_cone_hyperplanes"></a>

---


#### secondary\_cone\_hyperplanes

```python
def secondary_cone_hyperplanes(via_circuits: bool = False,
                               verbosity: int = 0) -> ArrayLike
```

Compute the hyperplanes of the secondary cone associated to this fan.
This cone has the interpretation:
    for a regular fan, a height h generates the fan iff it is in the
    relative interior of the secondary cone.

Irregular fans do not have heights generating them and thus do not have
secondary cones. One way to check regularity of a simplicial fan (i.e.,
a triangulation) is to attempt to construct the secondary cone. This
should be solid (i.e., full-dimensional). If the output cone is
non-solid, then the fan is irregular.

IRREGULARITY CHECKING ONLY WORKS IF `via_circuits=False`. WHEN
ATTEMPTING TO COMPUTE THE SECONDARY CONE OF AN IRREGULAR FAN USING
CIRCUITS, ONE CAN GET A FULL-DIMENSIONAL CONE!!!

Parameters
----------
via_circuits : bool, optional
    Whether to use circuits to compute the secondary cone.
    Should always be correct if the fan is regular but
    dangerous/not correct for checking irregularity...
    Alternative is local folding. Defaults to False.
verbosity : int, optional
    The verbosity level. Higher is more verbose. Defaults to 0.

Returns
-------
H : ArrayLike
    An array of hyperplanes, H, defining the cone as {x: Hx>=0}

<a id="fan.make_fine"></a>

---


#### make\_fine

```python
def make_fine(fan: Fan) -> Fan
```

Convert a fan to a fine fan by linear flipping

Parameters
----------
fan : Fan
    The initial fan.

Returns
-------
out : Fan
    A fine fan, from linearly flipping the original fan.

<a id="fan.flip_subgraph"></a>

---


#### flip\_subgraph

```python
def flip_subgraph(
        seed,
        max_flips: int | None = None,
        only_fine: bool = False,
        only_regular: bool = True,
        only_pc_triang: bool = False,
        compute_node_labels: bool = False,
        verbosity: int = 0) -> tuple[nx.Graph, list[Fan], list[dict]]
```

Compute the flip graph centered at some input 'seed' triangulation.

Optionally, allow restrictions including only allowing triangulations
    - that are fewer than `max_flips` from the seed,
    - that are fine (use all vectors),
    - that are regular, and
    - that consist of triangulations which 'respect the point configuration'
      (i.e., also correspond to a fine, star triangulation of the
      associated point configuration).
If any such restrictions are applied but the seed doesn't obey them, then an
empty graph will be output.

Parameters
----------
seed : Fan
    The seed triangulation (center of flip graph).
max_flips : int | None, optional
    Max number of flips to consider from seed. Defaults to None.
only_fine : bool, optional
    Whether to restrict to fine triangulations. Defaults to False.
only_regular : bool, optional
    Whether to restrict to regular triangulations. Defaults to True.
only_pc_triang : bool, optional
    Whether to restrict to triangulations which 'respect the point
    configuration'. Defaults to False.
compute_node_labels : bool, optional
    Whether to compute 'labels' for the nodes indicating whether the
    triangulation is fine, regular, and/or respects the point
    configuration. Defaults to False.
verbosity : int, optional
    The verbosity level. Higher is more verbose. Defaults to 0.

Returns
-------
out : tuple[nx.Graph, list[Fan], list[dict]]
    - The flip graph as a networkx.Graph object.
    - A list of the triangulations
    - A list of the labels for each triangulation (labels are a
      dictionary from the property to a bool)

<a id="circuits"></a>

---


# circuits

<a id="circuits.Circuit"></a>

---


## Circuit Objects

```python
class Circuit()
```

This class is a helper data structure to contain a single circuit of some
vector configuration.

Constructs a `Circuit` object describing a circuit of a vector
configuration.

Parameters
----------
vc : VectorConfiguration
    The ambient vector configuration.
Z : Iterable[int]
    The support of the circuit.
Zpos : Iterable[int]
    The 'positive' side of the circuit.
Zneg : Iterable[int]
    The 'negative' side of the circuit.
lmbda : ArrayLike
    A dependency vector demonstrating the circuit.
signature : tuple[int, int]
    The signature (|Zpos|, |Zneg|) of the circuit.

<a id="circuits.Circuit.__init__"></a>

---


#### \_\_init\_\_

```python
def __init__(vc, Z, Zpos, Zneg, lmbda, signature) -> None
```

Initializes a `Circuit` object.

Parameters
----------
vc : VectorConfiguration
    The ambient vector configuration.
Z : Iterable[int]
    The support of the circuit.
Zpos : Iterable[int]
    The 'positive' side of the circuit.
Zneg : Iterable[int]
    The 'negative' side of the circuit.
lmbda : ArrayLike
    A dependency vector demonstrating the circuit.
signature : tuple[int, int]
    The signature (|Zpos|, |Zneg|) of the circuit.

<a id="circuits.Circuits"></a>

---


## Circuits Objects

```python
class Circuits()
```

This class is a helper data structure to contain the circuits of some
vector configuration.

Constructs a `Circuits` object describing all circuits of some VC.

<a id="circuits.Circuits.__init__"></a>

---


#### \_\_init\_\_

```python
def __init__() -> None
```

Initializes a `Circuits` object.

<a id="circuits.Circuits.__getitem__"></a>

---


#### \_\_getitem\_\_

```python
def __getitem__(label_inds: Iterable[int]) -> Circuit | int
```

Get the circuit corresponding to the indicated indices.

Parameters
----------
label_inds : Iterable[int]
    The iterable of vector/label indices.

Returns
-------
out : Circuit | int
    Cases
        - if indices correspond to known circuit -> the `Circuit`
        - if indices correspond to non-circuit   -> -1
        - if indices aren't known                -> 0

<a id="circuits.Circuits.set_circuit"></a>

---


#### set\_circuit

```python
def set_circuit(circuit: Circuit, verbosity: int = 0) -> None
```

Set the circuit properties corresponding to the indicated indices.

Parameters
----------
circuit : Circuit
    A Circuit object.
verbosity : int, optional
    The verbosity level. Defaults to 0.

<a id="circuits.Circuits.set_non_dependency"></a>

---


#### set\_non\_dependency

```python
def set_non_dependency(label_inds: Iterable[int], verbosity: int = 0) -> None
```

Record a set of points that is not dependent

Parameters
----------
label_inds : Iterable[int]
    The iterable of vector/label indices.
verbosity : int, optional
    The verbosity level. Defaults to 0.

<a id="circuits.Circuits.values"></a>

---


#### values

```python
def values() -> Iterable[Circuit]
```

Get the values (the actual circuits)

Returns
-------
out : Iterable[Circuit]
    The circuits.

<a id="circuits.Circuits.copy"></a>

---


#### copy

```python
def copy() -> Circuits
```

Copy the circuits object

Returns
-------
out : Circuits
    A copy of the circuits.

<a id="circuits.Circuits.pop"></a>

---


#### pop

```python
def pop(*args, **kwargs)
```

Pop an element from the circuits dict

<a id="circuits.Circuits.encode"></a>

---


#### encode

```python
def encode(label_inds: Iterable[int]) -> int
```

Convert an iterable of integers to a binary vector, b, such that
    b_i = 1 <=> i in label_inds

Parameters
----------
label_inds : Iterable[int]
    The iterable of integers.

Returns
-------
out : int
    The encoding

<a id="circuits.Circuits.decode"></a>

---


#### decode

```python
def decode(encoding) -> list[int]
```

Convert a binary vector b to a list of integers such that
    b_i = 1 <=> i in label_inds

Parameters
----------
encoding : int
    The encoding to map to label indices

Returns
-------
out : list[int]
    The label indices

<a id="circuits.Circuits.is_superset"></a>

---


#### is\_superset

```python
def is_superset(setA, setB) -> bool
```

Check if the set encoded by setA is a superset of setB.

Parameters
----------
setA : int
    The candidate-superset encoding.
setB : int
    The candidate-subset encoding.

Returns
-------
out : bool
    Whether setA is a superset of setB.

<a id="circuits.Circuits.is_subset"></a>

---


#### is\_subset

```python
def is_subset(setA: int, setB: int) -> bool
```

Check if the set encoded by setA is a subset of setB.

Parameters
----------
setA : int
    The candidate-subset encoding.
setB : int
    The candidate-superset encoding.

Returns
-------
out : bool
    Whether setA is a subset of setB.

<a id="__init__"></a>

---


# \_\_init\_\_

