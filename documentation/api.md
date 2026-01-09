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

**Description:**
Constructs a `Circuit` object describing a circuit of a vector
configuration. This is handled by the hidden [`__init__`](#__init__)
function.

**Arguments:**
- `vc`:        The ambient vector configuration.
- `Z`:         The support of the circuit.
- `Zpos`:      The 'positive' side of the circuit.
- `Zpos`:      The 'negative' side of the circuit.
- `lmbda`:     A dependency vector demonstrating the circuit.
- `signature`: The signature (|Zpos|, |Zneg|) of the circuit.

**Returns:**
Nothing.

<a id="circuits.Circuit.__init__"></a>

---


#### \_\_init\_\_

```python
def __init__(vc, Z, Zpos, Zneg, lmbda, signature)
```

**Description:**
Initializes a `Circuit` object.

**Arguments:**
- `vc`:        The ambient vector configuration.
- `Z`:         The support of the circuit.
- `Zpos`:      The 'positive' side of the circuit.
- `Zpos`:      The 'negative' side of the circuit.
- `lmbda`:     A dependency vector demonstrating the circuit.
- `signature`: The signature (|Zpos|, |Zneg|) of the circuit.

**Returns:**
Nothing.

<a id="circuits.Circuits"></a>

---


## Circuits Objects

```python
class Circuits()
```

This class is a helper data structure to contain the circuits of some
vector configuration.

**Description:**
Constructs a `Circuits` object describing all circuits of some VC. This
is handled by the hidden [`__init__`](#__init__) function.

**Arguments:**
None.

**Returns:**
Nothing.

<a id="circuits.Circuits.__init__"></a>

---


#### \_\_init\_\_

```python
def __init__()
```

**Description:**
Initializes a `Circuits` object.

**Arguments:**
None.

**Returns:**
Nothing.

<a id="circuits.Circuits.__getitem__"></a>

---


#### \_\_getitem\_\_

```python
def __getitem__(label_inds: Iterable[int]) -> Union["Circuit", int]
```

**Description:**
Get the circuit corresponding to the indicated indices.

**Arguments:**
- `label_inds`: The iterable of vector/label indices.

**Returns:**
Cases
- if indices correspond to known circuit -> the `Circuit`
- if indices correspond to non-circuit   -> -1
- if indices aren't known                -> 0

<a id="circuits.Circuits.set_circuit"></a>

---


#### set\_circuit

```python
def set_circuit(circuit: "Circuit", verbosity: int = 0) -> None
```

**Description:**
Set the circuit properties corresponding to the indicated indices.

**Arguments:**
- `circuit`:   Dict describing the circuit.
- `verbosity`: The verbosity level.

**Returns:**
Nothing.

<a id="circuits.Circuits.set_non_dependency"></a>

---


#### set\_non\_dependency

```python
def set_non_dependency(label_inds: Iterable[int], verbosity: int = 0) -> None
```

**Description:**
Record a set of points that is not dependent

**Arguments:**
- `label_inds`: The iterable of vector/label indices.
- `verbosity`:  The verbosity level.

**Returns:**
Nothing.

<a id="circuits.Circuits.values"></a>

---


#### values

```python
def values() -> Iterable["Circuit"]
```

**Description:**
Get the values (the actual circuits)

**Arguments:**
None

**Returns:**
The circuits.

<a id="circuits.Circuits.copy"></a>

---


#### copy

```python
def copy() -> "Circuits"
```

**Description:**
Copy the circuits object

**Arguments:**
None

**Returns:**
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

**Description:**
Convert an iterable of integers to a binary vector, b, such that
b_i = 1 <=> i in label_inds

**Arguments:**
- `label_inds`: The iterable of integers.

**Returns:**
The encoding

<a id="circuits.Circuits.decode"></a>

---


#### decode

```python
def decode(encoding) -> list[int]
```

**Description:**
Convert a binary vector b to a list of of integers such that
b_i = 1 <=> i in label_inds

**Arguments:**
- `encoding`: The encoding to map to label indices

**Returns:**
The label indices

<a id="circuits.Circuits.is_superset"></a>

---


#### is\_superset

```python
def is_superset(setA, setB) -> bool
```

**Description:**
Check if the set encoded by setA is a superset of setB.

**Arguments:**
- `setA`: The candidate-superset encoding.
- `setB`: The candidate-subset encoding.

**Returns:**
Whether setA is a superset of setB.

<a id="circuits.Circuits.is_subset"></a>

---


#### is\_subset

```python
def is_subset(setA: int, setB: int) -> bool
```

**Description:**
Check if the set encoded by setA is a subset of setB.

**Arguments:**
- `setA`: The candidate-superset encoding.
- `setB`: The candidate-subset encoding.

**Returns:**
Whether setA is a subset of setB.

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

**Description:**
Constructs a `Fan` object describing a triangulation of a lattice vector
configuration. This is handled by the hidden [`__init__`](#__init__)
function.

This class is *not* intended to be called directly. Instead, it is meant to
be called through VectorConfiguration.triangulate.

**Arguments:**
- `vc`:      The ambient vector configuration that this fan is over.
- `cones`:   The cones defining the fan. Each cone is a collection of
integer labels.
- `heights`: The heights defining the fan, if it is regular. Can be
computed later.

**Returns:**
Nothing.

<a id="fan.Fan.__init__"></a>

---


#### \_\_init\_\_

```python
def __init__(vc: "VectorConfiguration",
             cones: list[list[int]],
             heights: list[float] = None)
```

**Description:**
Initializes a `Fan` object.

**Arguments:**
- `vc`:      The ambient vector configuration that this fan is over.
- `cones`:   The cones defining the fan. Each cone is a collection of
integer labels.
- `heights`: The heights defining the fan, if it is regular. Can be
computed later.

**Returns:**
Nothing.

<a id="fan.Fan.__repr__"></a>

---


#### \_\_repr\_\_

```python
def __repr__() -> str
```

**Description:**
String representation of the Fan.
(more detailed than __str__)

**Arguments:**
None.

**Returns:**
String representation of the object.

<a id="fan.Fan.__str__"></a>

---


#### \_\_str\_\_

```python
def __str__() -> str
```

**Description:**
String description of the Fan.
(less detailed than __repr__ but more readable)

**Arguments:**
None.

**Returns:**
String description of the object.

<a id="fan.Fan.__hash__"></a>

---


#### \_\_hash\_\_

```python
def __hash__() -> int
```

**Description:**
Hash for the fan. Defined by hashing vector configuration and the cones.

**Arguments:**
None.

**Returns:**
The hash.

<a id="fan.Fan.__eq__"></a>

---


#### \_\_eq\_\_

```python
def __eq__(o: "Fan") -> bool
```

**Description:**
Equality checking between two Fan objects.

**Arguments:**
- `o`: The other Fan to compare against.

**Returns:**
True if self==o. False if self!=o.

<a id="fan.Fan.__ne__"></a>

---


#### \_\_ne\_\_

```python
def __ne__(o: "Fan") -> bool
```

**Description:**
Inequality checking between two Fan objects.

**Arguments:**
- `o`: The other Fan to compare against.

**Returns:**
True if self!=o. False if self==o.

<a id="fan.Fan.vector_config"></a>

---


#### vector\_config

```python
@property
def vector_config() -> "VectorConfiguration"
```

**Description:**
Returns the associated vector configuration.

**Arguments:**
None.

**Returns:**
The associated vector configuration.

<a id="fan.Fan.labels"></a>

---


#### labels

```python
@property
def labels() -> tuple[int]
```

**Description:**
Returns the labels of the vectors in the VC.

**Arguments:**
None.

**Returns:**
The labels of the vectors in the VC.

<a id="fan.Fan.used_labels"></a>

---


#### used\_labels

```python
@property
def used_labels() -> tuple[int]
```

**Description:**
Returns the labels of the vectors in the VC used by cones in the Fan.

**Arguments:**
None.

**Returns:**
The labels of the vectors in the VC used by cones in the Fan.

<a id="fan.Fan.labels_to_cones"></a>

---


#### labels\_to\_cones

```python
@property
def labels_to_cones() -> dict[int, set[tuple[int]]]
```

**Description:**
Returns a map from vector labels to the cones the vector appears in.

**Arguments:**
None.

**Returns:**
A map from vector label to a set of cones (as tuples of indices) that
the vector appears in.

<a id="fan.Fan.ambient_dim"></a>

---


#### ambient\_dim

```python
@property
def ambient_dim() -> int
```

**Description:**
Returns the ambient dimension of the VC.

**Arguments:**
None.

**Returns:**
The ambient dimension of the VC.

<a id="fan.Fan.dim"></a>

---


#### dim

```python
@property
def dim() -> int
```

**Description:**
Returns the dimension of the VC. I.e., the dimension of the subspace
spanned by the vectors.

**Arguments:**
None.

**Returns:**
The dimension of the VC.

<a id="fan.Fan.vectors"></a>

---


#### vectors

```python
def vectors(which: int | Iterable[int] = None,
            lifted: bool = False) -> "ArrayLike"
```

**Description:**
Returns the vectors, optionally only those with given labels. Also,
optionally, give the vectors lifted by the heights (if the Fan is
regular).

**Arguments:**
- `which`:  Either a single label, for which the single corresponding
vector will be returned, or a list of labels. If not
provided, then all vectors are returned.
- `lifted`: Whether to give the lifted vectors.

**Returns:**
The corresponding vector(s), in order specified by which.

<a id="fan.Fan.cones"></a>

---


#### cones

```python
def cones(as_rays: bool = False,
          as_hyps: bool = False,
          as_inds: bool = False,
          ind_offset: int = 0) -> Union[tuple[tuple[int]], list["ArrayLike"]]
```

**Description:**
Returns the cones in the fan in a variety of formats. They are:
- (default) as a tuple of labels
- (as_rays=True) as an array whose rows are the generators
- (as_hyps=True) as an array whose rows are hyperplane normals
- (as_inds=True) as a tuple of indices
Optionally, allow an offset to the indices.

**Arguments:**
- `as_rays`:    Whether to return the cones as their generators.
- `as_hyps`:    Whether to return the cones as their hyperplanes.
- `as_inds`:    Whether to return the cones as indices (not labels).
- `ind_offset`: An additive offset for the indices

**Returns:**
The corresponding vector(s), in order specified by which.

<a id="fan.Fan.facets"></a>

---


#### facets

```python
def facets() -> dict[tuple[int], list[tuple[int]]]
```

**Description:**
Returns the facets of the cones. Save them as a dictionary from facet
labels to a list of containing cones, stored by their labels.

**Arguments:**
None.

**Returns:**
A dictionary from facet labels to a list of containing cones.

<a id="fan.Fan.is_valid"></a>

---


#### is\_valid

```python
def is_valid(verbosity: int = 0) -> bool
```

**Description:**
Return whether or not the cones define a valid polyhedral fan.

Follows cor. 4.5.13 of "Triangulations" by De Loera, Rambau, Santos.

**Arguments:**
- `verbosity`: The verbosity level. Higher is more verbose.

**Returns:**
True if the cones define a valid fan. False otherwise.

<a id="fan.Fan.respects_ptconfig"></a>

---


#### respects\_ptconfig

```python
def respects_ptconfig() -> bool
```

**Description:**
Return whether or not the fan also defines a (star) subdivision of the
original underlying point configuration.

**Arguments:**
None.

**Returns:**
True if the fan defines a subdivision of the point configuration. False
otherwise.

<a id="fan.Fan.is_triangulation"></a>

---


#### is\_triangulation

```python
def is_triangulation() -> bool
```

**Description:**
Return whether or not the fan is a triangulation (not a subdivision).

**Arguments:**
None.

**Returns:**
True if the fan is a triangulation. False otherwise.

<a id="fan.Fan.is_fine"></a>

---


#### is\_fine

```python
def is_fine() -> bool
```

**Description:**
Return whether or not the fan is fine.

**Arguments:**
None.

**Returns:**
True if the fan is fine. False otherwise.

<a id="fan.Fan.is_regular"></a>

---


#### is\_regular

```python
def is_regular() -> bool
```

**Description:**
Return whether or not the fan is regular.

**Arguments:**
None.

**Returns:**
True if the fan is regular. False otherwise.

<a id="fan.Fan.heights"></a>

---


#### heights

```python
def heights() -> list[float] | None
```

**Description:**
Return some heights defining the cone, if it is regular. Else, return
None.

**Arguments:**
None.

**Returns:**
True if the fan is regular. False otherwise.

<a id="fan.Fan.contains"></a>

---


#### contains

```python
def contains(c: Iterable[int] | Iterable[Iterable[int]]) -> bool
```

**Description:**
Check if any cone (specified by its labels) is contained in the fan.
The cone need not be solid. Can also be called for a collection of
cones, in which case the check is if all cones are contained in the fan.

**Arguments:**
- `c`: The cone(s). Either a single collection of cone, specified by
an iterable of labels, or a collection of cones, each specified
by an iterable of labels.

**Returns:**
Whether (all) cone(s) are contained in the fan.

<a id="fan.Fan.circuit"></a>

---


#### circuit

```python
def circuit(labels: Iterable[int] = None,
            enforce_positive: int = None,
            lmbda: Iterable[float] = None,
            check_containment: bool = True,
            save_circuits_in_vc: bool = False,
            verbosity: int = 0) -> "Circuit"
```

**Description:**
Format/compute the circuit corresponding to the specified labels.

**Arguments:**
- `labels`:              Labels indicating the vectors in the circuit.
- `enforce_positive`:    A label to enforce is in Zpos.
- `lmbda`:               A dependency demonstrating the circuit.
- `check_containment`:   Whether to check that this fan contains every
cone in the positive triangulation, Tpos.
- `save_circuits_in_vc`: Whether to save circuits... best to keep True
for most circumstances.
- `verbosity`:           The verbosity level. Higher is more verbose.

**Returns:**
Circuit object containing
- the support of the circuit as property 'Z',
- the signed circuit as property 'Zpos' and 'Zned',
- the dependency as property 'lmbda', and
- the signature as property 'signature'.

<a id="fan.Fan.circuits"></a>

---


#### circuits

```python
def circuits(facets: dict[Iterable[int], Iterable[Iterable[int]]] = None,
             verbosity: int = 0) -> list["Circuit"]
```

**Description:**
Compute all circuits associated to this fan (i.e., those 'embedded' in
this fan). All will be oriented such that the positive triangulation
(i.e., Tpos/T_+) is embedded in the fan. This enables us to directly
interpret lambda as the normal in the secondary cone.

**Arguments:**
- `facets`:    The facets of the fan (not just the VC...). I.e., codim-1
cones.
- `verbosity`: The verbosity level. Higher is more verbose.

**Returns:**
A list of Circuit objects for all circuits embedded in the fan.

<a id="fan.Fan.star"></a>

---


#### star

```python
def star(cell: Iterable[int], old_way: bool = False) -> Iterable[tuple[int]]
```

**Description:**
Compute the star of some cell. This is the subcomplex of all cones
containing the cell (and their faces)

**Arguments:**
- `cell`:    The cell of interest.
- `old_way`: Whether to do the computation in an old/slow manner.

**Returns:**
A list of all solid cones (as tuples of ints) containing the cell.

<a id="fan.Fan.link"></a>

---


#### link

```python
def link(cell: Iterable[int]) -> list[tuple[int]]
```

**Description:**
Compute the link of some cell. This is the subcomplex of all cones in
the star that don't intersect the cell.

**Arguments:**
- `cell`: The cell of interest.

**Returns:**
A list of all solid cones (as tuples of ints) containing the cell.

<a id="fan.Fan.embed"></a>

---


#### embed

```python
def embed(cell: Iterable[int],
          link: Iterable[Iterable[int]]) -> list[tuple[int]]
```

**Description:**
Embed some cell into the Fan bu combining it with each cell in the link.

**Arguments:**
- `cell`: The cell of interest.
- `link`: The link of said cell.

**Returns:**
A list of solid cones representing the embedding of the cell into the
Fan via the link.

<a id="fan.Fan.flip"></a>

---


#### flip

```python
def flip(circ: "Circuit",
         formal: bool = True,
         verbosity: int = 0) -> Union["Fan", tuple[tuple[int]]]
```

**Description:**
Make a flip across a circuit.

**Arguments:**
- `circ`:      The circuit to flip through.
- `formal`:    Whether to return a formal Fan (otherwise, just a tuple
of cones).
- `verbosity`: The verbosity level. Higher is more verbose.

**Returns:**
The flipped Fan.

<a id="fan.Fan.flip_linear"></a>

---


#### flip\_linear

```python
def flip_linear(
    h_target: Iterable[float] = None,
    direction: Iterable[float] = None,
    h_init: Iterable[float] = None,
    max_N_flips: int = None,
    stop_at_deletion: bool = True,
    stop_at_pct: bool = False,
    check_regularity: bool = True,
    record_fans: bool = False,
    record_circs: bool = False,
    hook_init: Callable = None,
    hook_flip: Callable = None,
    eps: float = 1e-8,
    verbosity: int = 0
) -> list[int | Exception, "ArrayLike", "Fan", "ArrayLike", int]
```

**Description:**
Compute all flips along the linear height homotopy
t*h_target + (1-t)*h_init
for t=0 increasing to t=1.

Allow early stops of this homotopy at a certain number `max_N_flips` of
flips. Also allow early stopping upon the following conditions
- (default True) reaching a deletion flip or
- (default False) hitting a fan that respects the point config.

**Arguments:**
(defining the homotopy)
- `h_target`:         The target heights.
- `direction`:        The direction to travel.
- `h_init`:           The initial heights (regular triangulations don't
have unique heights, even up to scaling... any h
in the secondary cone is valid. If this is left
unset, then arbitrary valid heights are chosen)
(early stopping)
- `max_N_flips`:      The maximum number of flips allowed.
- `stop_at_deletion`: Whether to early-terminate the homotopy at any
deletion flip seen.
- `stop_at_pct`:      Whether to early-terminate the homotopy at any
fan that respects the point configuration.
(sanity checks)
- `check_regularity`: This method is inherently regular (it uses
heights...). We can check the regularty of the
initial fan.
(record keeping)
- `record_fans`:      Whether to record the fans seen along the
homotopy.
- `record_circs`:     Whether to record the circuits flipped along the
homotopy.
(numerical parameters)
- `eps`:              A small number for an allowed violation of heights
landing outside the secondary fan (in case the
heights 'truly' landed on a wall of the secondary
fan). Such violations are naturally resolved by
pulling heights back into the secondary fan.
(diagnostics)
- `verbosity`:        The verbosity level. Higher is more verbose.

**Returns:**
- The status of the homotopy. Either 1 (if successful) or an Exception.
- The current heights at the end of the homotopy. Not always h_target.
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
) -> tuple[list[Union["Fan", tuple[tuple[int]]]], list["Circuit"]]
```

**Description:**
Compute the neighboring fans (those reachable by a single flip).

Allow restrictions to only fine fans.

**Arguments:**
- `only_fine`: Whether to only compute/return fine neighbors
- `formal`:    Whether to return the neighbors as formal fans (if
False, just return cones).
- `verbosity`: The verbosity level. Higher is more verbose.

**Returns:**
- The neighbors, either as formal Fan objects or as collections of
cones (each cone a collection of labels)
- The circuits flipped to get the corresponding neighbors.

<a id="fan.Fan.secondary_cone_hyperplanes"></a>

---


#### secondary\_cone\_hyperplanes

```python
def secondary_cone_hyperplanes(via_circuits: bool = False,
                               verbosity: int = 0) -> "ArrayLike"
```

**Description:**
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

**Arguments:**
- `via_circuits`:      Whether to use circuits to compute the secondary
cone. Should always be correct if the fan is
regular but dangerous/not correct for checking
irregularity... Alternative is local folding.
- `verbosity`:         The verbosity level. Higher is more verbose.

**Returns:**
An array of hyperplanes, H, defining the cone as {x: Hx>=0}

<a id="fan.flip_subgraph"></a>

---


#### flip\_subgraph

```python
def flip_subgraph(
        seed,
        max_flips: int = None,
        only_fine: bool = False,
        only_regular: bool = True,
        only_pc_triang: bool = False,
        compute_node_labels: bool = False,
        verbosity: int = 0
) -> tuple["networkx.Graph", list["Fan"], list[dict]]
```

**Description:**
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

**Arguments:**
- `seed`:                The seed triangulation (center of flip graph).
- `max_flips`:           Max number of flips to consider from seed.
- `only_fine`:           Whether to restrict to fine triangulations.
- `only_regular`:        Whether to restrict to regular triangulations.
- `only_pc_triang`:      Whether to restrict to triangulations which
'respect the point configuration'.
- `compute_node_labels`: Whether to compute 'labels' for the nodes
indicating whether the triangulation is fine,
regular, and/or respects the point configuration.
- `verbosity`:           The verbosity level. Higher is more verbose.

**Returns:**
- The flip graph as a networkx.Graph object.
- A list of the triangulations
- A list of the labels for each triangulation (labels are a dictionary from
the property to a bool)

<a id="util"></a>

---


# util

<a id="util.gcd"></a>

---


#### gcd

```python
def gcd(vals: list[float], max_denom: float = 10**6) -> float
```

**Description:**
Computes the 'GCD' of a collection of floating point numbers.
This is the smallest number, g, such that g*values is integral.

This is computed by
1) converting `values` to be rational [n0/d0, n1/d1, ...],
2) computing the LCM, l, of [d0, d1, ...],
3) computing the GCD, g', of [l*n0/d0, l*n1/d1, ...], and then
4) returning g=g'/l.

**Arguments:**
- `vals`:      The numbers to compute the GCD of.
- `max_denom`: Assert |di| <= max_denom

**Returns:**
The minimum number g' such that g'*vals is integral.

<a id="util.primitive"></a>

---


#### primitive

```python
def primitive(vec: list[float], max_denom=10**10)
```

**Description:**
Computes the primitive vector associated to the input ray {c*vec: c>=0}.
Very similar to the gcd function.

This is equivalent to
vec/gcd(vec)
but just uses a rational representation.

**Arguments:**
- `vec`:       A vector defining the ray {c*vec: c>=0}
- `max_denom`: Assert |di| <= max_denom

**Returns:**
The primitive vector along the ray.

<a id="util.lerp"></a>

---


#### lerp

```python
def lerp(p0: "ArrayLike", p1: "ArrayLike", t: float) -> "ArrayLike"
```

**Description:**
Computes the point specified by t along the line passing through p0 and p1.

Particular values:
-) t=0   -> p0
-) t=0.5 -> (p0+p1)/2
-) t=1   -> p1

**Arguments:**
- `p0`: One point.
- `p1`: The other point.
- `t`:  Parameter specifing where along the line Conv({p0, p1}) to return.

**Returns:**
The point p0 + t*(p1-p0).

<a id="util.hyperplane_inter"></a>

---


#### hyperplane\_inter

```python
def hyperplane_inter(p0: "ArrayLike", p1: "ArrayLike",
                     n: "ArrayLike") -> float
```

**Description:**
Computes the distance, t, along the line p0 + t*(p1-p0) that intersects the
hyperplane {x: n.x=0}.

This is simply computed as
0 = n.[p0 + t*(p1-p0)]
0 = n.p0 + t*n.(p1-p0)
t = -n.p0 / n.(p1-p0)

**Arguments:**
- `p0`: One point.
- `p1`: The other point.
- `n`:  The hyperplane normal.

**Returns:**
The distance t such that p0 + t*(p1-p0) lands upon the hyperplane.

<a id="util.first_hit"></a>

---


#### first\_hit

```python
def first_hit(p0: "ArrayLike",
              p1: "ArrayLike",
              H: "ArrayLike",
              max_dist: float = 1,
              verbosity: int = 0) -> (int, float)
```

**Description:**
Given a point p0 in a convex cone {x: Hx>=0}, find the first hyperplane hit
along the direction (p1-p0). I.e, the first intersection of the ray
{p0+t*(p1-p0): t>=0} with the cones bounding hyperplanes.

Allow violated hyperplanes (i.e., n such that n.p0 < 0) but ignore them.

**Arguments:**
- `p0`:         One point.
- `p1`:         The other point.
- `H`:          An array of hyperplane normals (as rows).
- `max_dist`:   Only consider intersections along the line segment
[p0, p0+max_dist*(p1-p0)]
- `verbosity`:  The verbosity level. Higher is more verbose.

**Returns:**
The index, i, of the first-hit hyperplane.
The distance, t, such that H[i].lerp(p0,p1,t) = 0.

<a id="util.dual_cone"></a>

---


#### dual\_cone

```python
def dual_cone(data: "ArrayLike") -> "ArrayLike"
```

**Description:**
Compute the data of the cone dual to the input 'primal' cone.

This can be thought of in a couple of equivalent ways, summarized in the
following table. E.g., if rays of the primal are input, then the
hyperplanes of the primal are output (or, equivalently, the rays of the
dual).

INPUT       | PRIMAL OUTPUT | DUAL OUTPUT
-----------------------------------------
rays        | hyperplanes   | rays
hyperplanes | rays          | hyperplanes


For simplicitly in the following discussion, take the convention that one
maps hyperplanes of the primal to rays of the primal.

**Arguments:**
- `data`: An array whose rows represent rays of the primal cone. (see table)

**Returns:**
An array whose rows represent hyperplanes of the primal cone. (see table)

<a id="util.cone_dim"></a>

---


#### cone\_dim

```python
def cone_dim(*, R: "ArrayLike" = None, H: "ArrayLike" = None) -> int
```

**Description:**
Return the dimension of the cone.

The cone is either specified via rays,
{R.T @ lambda: lambda>=0},
or via hyperplanes,
{x: H @ x>=0}.

**Arguments:**
- `R`: The rays of the cone as rows.
- `H`: The hyperplanes defining the cone.

**Returns:**
The dimension of the cone

<a id="util.is_solid"></a>

---


#### is\_solid

```python
def is_solid(*, R: "ArrayLike" = None, H: "ArrayLike" = None) -> int
```

**Description:**
Return whether the cone is full-dimensional.

The cone is either specified via rays,
{R.T @ lambda: lambda>=0},
or via hyperplanes,
{x: H @ x>=0}.

**Arguments:**
- `R`: The rays of the cone as rows.
- `H`: The hyperplanes defining the cone.

**Returns:**
The dimension of the cone

<a id="util.contains"></a>

---


#### contains

```python
def contains(*,
             p: "ArrayLike",
             R: "ArrayLike" = None,
             H: "ArrayLike" = None) -> bool
```

**Description:**
Return if the point p is contained in the cone.

The cone is either specified via rays,
{R.T @ lambda: lambda>=0},
or via hyperplanes,
{x: H @ x>=0}.

**Arguments:**
- `R`: The rays of the cone as rows.
- `H`: The hyperplanes defining the cone.

**Returns:**
Whether p is contained in the cone.

<a id="util.find_interior_point"></a>

---


#### find\_interior\_point

```python
def find_interior_point(*,
                        R: "ArrayLike" = None,
                        H: "ArrayLike" = None,
                        stretching: float = 1,
                        nonneg: bool = False,
                        verbosity: int = 0) -> Union["ArrayLike", None]
```

**Description:**
Returns a point p in the relative interior of a cone. The cone can be
specified either via its rays or its generators.

If no point p exists, return `None`.

Modified from CYTools' `Cone.find_interior_point`.

**Arguments:**
- `R`:          Generators defining the cone.
- `H`:          Hyperplanes defining the cone.
- `stretching`: How far p must be from any hyperplane.
- `nonneg`:     Whether to restrict to non-negative vectors.
- `verbosity`:  The verbosity level.

**Returns:**
A point p in the strict interior.

<a id="__init__"></a>

---


# \_\_init\_\_

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

**Description:**
Constructs a `VectorConfiguration` object describing a lattice vector
configuration. This is handled by the hidden [`__init__`](#__init__) function.

**Arguments:**
- `vectors`:    The vectors defining the VC.
- `labels`:     A list of labels for the vectors. Only integral labels are
allowed.
- `eps`:        Threshold for checking for non-integral vectors.
- `gale_basis`: An optional basis for the gale transform. If provided, then
the gale transform will be put a basis such that the
submatrix given by these labels equals the identity.

**Returns:**
Nothing.

<a id="vectorconfig.VectorConfiguration.__init__"></a>

---


#### \_\_init\_\_

```python
def __init__(vectors: "ArrayLike",
             labels: Iterable[int] = None,
             eps: float = 1e-4,
             gale_basis: Iterable[int] = None) -> None
```

**Description:**
Initializes a `VectorConfiguration` object.

**Arguments:**
- `vectors`:    The vectors defining the VC.
- `labels`:     A list of integer labels for the vectors. Only integral
labels are allowed.
- `eps`:        Threshold for checking for non-integral vectors.
- `gale_basis`: An optional basis for the gale transform. If provided,
then the gale transform will be put a basis such that
the submatrix given by these labels equals the identity.

**Returns:**
Nothing.

<a id="vectorconfig.VectorConfiguration.__repr__"></a>

---


#### \_\_repr\_\_

```python
def __repr__() -> str
```

**Description:**
String representation of the VectorConfiguration.
(more detailed than __str__)

**Arguments:**
None.

**Returns:**
String representation of the object.

<a id="vectorconfig.VectorConfiguration.__str__"></a>

---


#### \_\_str\_\_

```python
def __str__() -> str
```

**Description:**
String description of the VectorConfiguration.
(less detailed than __repr__ but more readable)

**Arguments:**
None.

**Returns:**
String description of the object.

<a id="vectorconfig.VectorConfiguration.__hash__"></a>

---


#### \_\_hash\_\_

```python
def __hash__() -> int
```

**Description:**
Hash for the vector configuration. Defined by hashing a dictionary from
labels to vectors.

**Arguments:**
None.

**Returns:**
The hash.

<a id="vectorconfig.VectorConfiguration.__ne__"></a>

---


#### \_\_ne\_\_

```python
def __ne__(o: "VectorConfiguration") -> bool
```

**Description:**
Inequality checking between two VectorConfiguration objects.

**Arguments:**
- `o`: The other VectorConfiguration to compare against.

**Returns:**
True if self!=o. False if self==o.

<a id="vectorconfig.VectorConfiguration.__eq__"></a>

---


#### \_\_eq\_\_

```python
def __eq__(o: "VectorConfiguration") -> bool
```

**Description:**
Equality checking between two VectorConfiguration objects.

**Arguments:**
- `o`: The other VectorConfiguration to compare against.

**Returns:**
True if self==o. False if self!=o.

<a id="vectorconfig.VectorConfiguration.copy"></a>

---


#### copy

```python
def copy() -> "VectorConfiguration"
```

**Description:**
Copy method.

**Arguments:**
None.

**Returns:**
A copy of the vector configuration.

<a id="vectorconfig.VectorConfiguration.labels"></a>

---


#### labels

```python
@property
def labels() -> tuple[int]
```

**Description:**
Returns the labels of the vectors in the VC.

**Arguments:**
None.

**Returns:**
The labels of the vectors in the VC.

<a id="vectorconfig.VectorConfiguration.labels_to_inds_dict"></a>

---


#### labels\_to\_inds\_dict

```python
@property
def labels_to_inds_dict() -> dict[int, int]
```

**Description:**
Returns the a dictionary mapping vector labels to their indices in the
vector configuration.

**Arguments:**
None.

**Returns:**
The mapping from labels to indices.

<a id="vectorconfig.VectorConfiguration.size"></a>

---


#### size

```python
@property
def size() -> int
```

**Description:**
Returns the number of the vectors in the VC.

**Arguments:**
None.

**Returns:**
The number of the vectors in the VC.

<a id="vectorconfig.VectorConfiguration.ambient_dim"></a>

---


#### ambient\_dim

```python
@property
def ambient_dim() -> int
```

**Description:**
Returns the ambient dimension of the VC.

**Arguments:**
None.

**Returns:**
The ambient dimension of the VC.

<a id="vectorconfig.VectorConfiguration.dim"></a>

---


#### dim

```python
@property
def dim() -> int
```

**Description:**
Returns the dimension of the VC. I.e., the dimension of the subspace
spanned by the vectors.

**Arguments:**
None.

**Returns:**
The dimension of the VC.

<a id="vectorconfig.VectorConfiguration.vectors"></a>

---


#### vectors

```python
def vectors(which: int | Iterable[int] = None) -> "ArrayLike"
```

**Description:**
Returns the vectors, optionally only those with given labels.

**Arguments:**
- `which`: Either a single label, for which the single corresponding
vector will be returned, or a list of labels. If not
provided, then all vectors are returned.

**Returns:**
The corresponding vector(s), in order specified by which.

<a id="vectorconfig.VectorConfiguration.vectors_to_labels"></a>

---


#### vectors\_to\_labels

```python
def vectors_to_labels(vectors: "ArrayLike") -> int | list[int]
```

**Description:**
Maps the vectors to their corresponding labels

**Arguments:**
- `vectors`: Either a single vector, for which the single corresponding
label will be returned, or a list of vectors.

**Returns:**
The corresponding label(s).

<a id="vectorconfig.VectorConfiguration.labels_to_inds"></a>

---


#### labels\_to\_inds

```python
def labels_to_inds(labels: Iterable[int],
                   ambient_labels: Iterable[int] = None,
                   offset: int = 0) -> int | Iterable[int]
```

**Description:**
Maps the labels to their indices in ambient_labels, optionally with a
fixed offset.

**Arguments:**
- `labels`:         The labels of interest.
- `ambient_labels`: The ambient labels to get the indices in. If None,
use all labels of the VectorConfiguration.
- `offset`:         Return i+offset for i the index of a label in
ambient_labels.

**Returns:**
The indices of the labels.

<a id="vectorconfig.VectorConfiguration.is_solid"></a>

---


#### is\_solid

```python
def is_solid() -> bool
```

**Description:**
Return whether or not the VC is full-dimensional.

**Arguments:**
None.

**Returns:**
True if the VC is full-dimensional. False otherwise.

<a id="vectorconfig.VectorConfiguration.is_totally_cyclic"></a>

---


#### is\_totally\_cyclic

```python
def is_totally_cyclic() -> bool
```

**Description:**
Return whether or not the VC is totally cyclic. That is, whether
self.conv() equals the subspace containing it (the supporting
hyperplane).

**Arguments:**
None.

**Returns:**
True if the VC is totally cyclic. False otherwise.

<a id="vectorconfig.VectorConfiguration.is_acyclic"></a>

---


#### is\_acyclic

```python
def is_acyclic() -> bool
```

**Description:**
Return whether or not the VC is acyclic. That is, whether there exists
some direction psi such that
psi.vi > 0 for all vi.

This is equivalent to defining the cone {x: vi.x >= 0} and checking if
it is full-dimensional.

**Arguments:**
None.

**Returns:**
True if the VC is acyclic. False otherwise.

<a id="vectorconfig.VectorConfiguration.support"></a>

---


#### support

```python
def support() -> "ArrayLike"
```

**Description:**
Get the support of the vector configuration as a hyperplane
representation.

**Arguments:**
None.

**Returns:**
The hyperplanes defining the support.

<a id="vectorconfig.VectorConfiguration.cone_contains"></a>

---


#### cone\_contains

```python
def cone_contains(cone_labels: Iterable[int],
                  vec_label: Iterable[int],
                  strict: bool = False) -> bool
```

**Description:**
Check if a cone, specified by cone_labels, contains a the ray specified
by vec_label.

I.e., if
H = self.cone(cone_labels).hyperplanes()
v = self.vectors(vec_label)
H@v >= int(strict)

**Arguments:**
- `cone_labels`: The labels of vectors defining the cone.
- `vec_label`:   The label of the vector to check.
- `strict`:      Whether to check if the vector is in the strict
interior.

**Returns:**
Whether the associated cone contains the vector.

<a id="vectorconfig.VectorConfiguration.gale"></a>

---


#### gale

```python
def gale(set_basis: bool = False) -> "ArrayLike"
```

**Description:**
Compute the gale transform of the config.

I.e., a basis of the null-space of the vectors.

**Arguments:**
- `set_basis`: Whether to set a particular basis of the Gale transform.

**Returns:**
The gale transform.

<a id="vectorconfig.VectorConfiguration.project"></a>

---


#### project

```python
def project(vec: "ArrayLike") -> "ArrayLike"
```

**Description:**
Project down a vector from height-space to chamber-space.

**Arguments:**
- `vec`: The height-space vector.

**Returns:**
The chamber-space vector.

<a id="vectorconfig.VectorConfiguration.jorp"></a>

---


#### jorp

```python
def jorp(vec: "ArrayLike") -> "ArrayLike"
```

**Description:**
Undo a projection from height-space to chamber-space.

I.e., map from chamber-space to height-space

**Arguments:**
- `vec`: The chamber-space vector.

**Returns:**
The chamber-space vector.

<a id="vectorconfig.VectorConfiguration.triangulate"></a>

---


#### triangulate

```python
def triangulate(heights: "ArrayLike" = None,
                cells: "ArrayLike" = None,
                tol: float = 1e-14,
                backend: str = None,
                check_heights: bool = True,
                verbosity: int = 0) -> "Fan"
```

**Description:**
Subdivide the vector configuration either by specified cells/simplices
or by heights.

**Arguments:**
- `heights`:       The heights to lift the vectors by.
- `cells`:         The cells to use in the triangulation.
- `tol`:           Numerical tolerance used for curing negative heights
- `backend`:       The lifting backend. Currently allowed to be "cgal"
or "ppl".
- `check_heights`: Whether to check that the heights land in the
secondary cone of the output triangulation.
- `verbosity`:     The verbosity level. Higher is more verbose

**Returns:**
The resultant subdivision.

<a id="vectorconfig.VectorConfiguration.all_triangulations"></a>

---


#### all\_triangulations

```python
def all_triangulations(only_fine: bool = False,
                       only_regular: bool = True,
                       verbosity: int = 0) -> list["Fan"]
```

**Description:**
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


**Arguments:**
- `only_fine`:    Whether to restrict to fine triangulations
- `only_regular`: Whether to restrict to regular triangulations
- `verbosity`:    The verbosity level. Higher is more verbose.

**Returns:**
A list of Fan objects, one for each triangulation of the VC.

<a id="vectorconfig.VectorConfiguration.random_triangulations_fast"></a>

---


#### random\_triangulations\_fast

```python
def random_triangulations_fast(
        method: str = "delaunay",
        h0: "ArrayLike" = None,
        sigma: float = 0.1,
        N: int = None,
        as_list: bool = False,
        attempts_per_triang: int = 1000,
        backend: str = None,
        seed: int = 0,
        verbosity: int = 0) -> Generator["Fan"] | list["Fan"]
```

**Description:**
Generate random regular triangulations by picking random heights.

**Arguments:**
- `method`:             Either "delaunay" or "isotropic". The former
picks heights around some input height (e.g.,
the Deulaunay heights). The latter picks
heights isotropically
- `h0`:                 The reference heights, for Delaunay method.
- `sigma`:              How big of a distribution to study around h0.
- `N`:                  The number of triangulations to generate. If
as_list, then code will keep track of all
triangulations, retrying at most
attempts_per_triang tries to get a new
triangulation until the list has N triangs.
O/w, then the first N height vectors are used
(regardless of duplicates).
- `as_list`:            Whether to return the triangulations as a list,
or as a generator.
- `attempts_per_triang`:Quit if we can't generate a new triangulation
after this many tries.
- `backend`:            The lifting backend. See
`VectorConfiguration.triangulate`.
- `seed`:               A random number seed.
- `verbosity`:          The verbosity level. Higher is more verbose.

**Returns:**
The random triangulations.

<a id="vectorconfig.VectorConfiguration.circuit"></a>

---


#### circuit

```python
def circuit(labels: Iterable[int],
            lmbda: Iterable = None,
            set_non_dependencies: bool = True,
            save_circuits: bool = True) -> "Circuit"
```

**Description:**
Format/compute the circuit corresponding to the specified labels.

**Arguments:**
- `labels`:               Labels indicating the vectors in the circuit.
- `lmbda`:                Vector demonstrating the dependence.
- `set_non_dependencies`: Whether to update our list of non-circuits.
- `save_circuits`:        Whether to save circuits... best to keep True
for most circumstances.

**Returns:**
Circuit object containing
- the support of the circuit as property 'Z',
- the signed circuit as property 'Zpos' and 'Zned',
- the dependency as property 'lmbda', and
- the signature as property 'signature'.

<a id="vectorconfig.VectorConfiguration.circuits"></a>

---


#### circuits

```python
def circuits(verbosity: int = 0) -> list["Circuit"]
```

**Description:**
Compute all possible circuits of this vector configuration.

**Arguments:**
- `verbosity`: The verbosity level. Higher is more verbose.

**Returns:**
A list of Circuit objects.

<a id="vectorconfig.VectorConfiguration.flip_graph"></a>

---


#### flip\_graph

```python
def flip_graph(max_flips: int = None,
               only_fine: bool = False,
               only_regular: bool = True,
               only_pc_triang: bool = False,
               compute_node_labels: bool = False,
               verbosity: int = 0) -> (nx.Graph, list["Fan"], list[str])
```

**Description:**
Compute the flip graph. Wrapper for flip_subgraph.

**Arguments:**
- `max_flips`:           The maximum number of flips to take from the
seed. If none is provided, then the entire flip
graph is calculated.
- `only_fine`:           Whether to only compute fine triangulations.
- `only_regular`:        Whether to only compute regular triangulations.
Note, we never will see irregular
triangulations that are not connected to
regular ones.
- `only_pc_triang`:      Whether to only compute triangulations that also
correspond to star triangulations of the
underlying point config.
- `compute_node_labels`: Whether to check whether each node is fine,
regular, and a PC triangulation.
- `verbosity`:           The verbosity level. Higher is more verbose.

**Returns:**
The Graph object, whose nodes correspond to (and have values equal to)
Fan objects.  The edge between nodes correspond to flips, and have
labels equal to the corresponding circuit.

<a id="vectorconfig.VectorConfiguration.secondary_fan"></a>

---


#### secondary\_fan

```python
def secondary_fan(only_fine: bool = False,
                  formal_fan: bool = False,
                  verbosity: int = 0)
```

**Description:**
Compute the secondary fan of the vector configuration.

**Arguments:**
- `only_fine`:  Restrict to fine triangulations.
- `formal_fan`: Save as a formal Fan object.
- `verbosity`:  The verbosity level. Higher is more verbose

**Returns:**
The secondary fan triangulations.

<a id="vectorconfig.VectorConfiguration.central_fan"></a>

---


#### central\_fan

```python
def central_fan() -> "Fan"
```

**Description:**
Generate the central fan of the vector configuration. Can be defined
as lifting each vector by a height of 1.

**Arguments:**
None.

**Returns:**
The central fan.

