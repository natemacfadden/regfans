# =============================================================================
#    Copyright (C) 2025  Nate MacFadden for the Liam McAllister Group
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
# =============================================================================
#
# -----------------------------------------------------------------------------
# Description:  This module contains a class designed to perform computations
#               on vector configurations.
# -----------------------------------------------------------------------------

# external imports
from collections.abc import Generator, Iterable
import copy
import flint
import itertools
import networkx as nx
import numpy as np
import scipy as sp
import triangulumancer

# local imports
from . import util, circuits, fan


class VectorConfiguration:
    """
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
    """
    def __init__(
        self,
        vectors: "ArrayLike",
        labels: Iterable[int] = None,
        eps: float = 1e-4,
        gale_basis: Iterable[int] = None,
    ) -> None:
        """
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
        """
        # sanitize vectors
        # ----------------
        self._vectors = np.array(vectors)

        # check if vectors are integral
        if np.issubdtype(self._vectors.dtype, np.integer):
            # vectors are of an integral type... automatically OK
            pass
        else:
            # vectors are not obviously integral... check them
            rounded_vecs = np.rint(self._vectors)
            if np.any(np.abs(self._vectors - rounded_vecs) > eps):
                raise ValueError("Only integral vectors are allowed")
            else:
                self._vectors = rounded_vecs.astype(int)

        # delete origin if it's included
        norms = np.linalg.norm(self._vectors, ord=1, axis=1)
        small_norm = np.where(norms < 0.5)[0]

        if len(small_norm):
            print(
                f"The vectors {[self._vectors[i] for i in small_norm]} "
                "all had too-small norms... deleting them..."
            )

            good_norm = [
                i for i in range(len(self._vectors)) if i not in small_norm
            ]

            self._vectors = self._vectors[good_norm, :]
            if labels is not None:
                labels = [labels[i] for i in good_norm]

        # get the labels
        # --------------
        if labels is None:
            # start labelling at 1
            # (to support construction of point configurations (in, e.g.,
            #  CYTools) and their associated triangulations from this VC, it's
            #  nice to reserve label 0 for the origin)
            labels = [i + 1 for i in range(len(self._vectors))]
        
        self._labels = tuple([label for label in labels])
        if not all([isinstance(l,int) for l in self._labels]):
            raise ValueError("Labels must be integral")

        self._standard_labels = (self._labels == tuple(range(1, self.size + 1)))

        # construct useful maps
        # ---------------------
        self._labels_to_vectors = {
            label: vec for label, vec in zip(self._labels, self._vectors)
        }
        self._vectors_to_labels = {
            tuple(vec): label for label, vec in zip(self._labels, self._vectors)
        }
        self._labels_to_inds = None

        # initialize other variables
        # --------------------------
        self._dim = None

        self._circuits = circuits.Circuits()
        self._computed_all_circuits = False
        self._refinements = dict()

        self._poly = dict()

        self._flip_graphs = dict()

        # allow setting of a particular basis of the Gale transform
        self._gale_basis    = None

        # allow caching of the Gale transform
        self._gale_in_basis = None
        self._gale          = None

    # defaults
    # ========
    def __repr__(self) -> str:
        """
        **Description:**
        String representation of the VectorConfiguration.
        (more detailed than __str__)

        **Arguments:**
        None.

        **Returns:**
        String representation of the object.
        """
        vecs = self.vectors().tolist()

        return (
            f"A {self.dim}-dimensional vector configuration consisting of the "
            f"following #{self.size} vectors: {vecs} "
            f"with labels: {self.labels}"
        )

    def __str__(self) -> str:
        """
        **Description:**
        String description of the VectorConfiguration.
        (less detailed than __repr__ but more readable)

        **Arguments:**
        None.

        **Returns:**
        String description of the object.
        """
        return (
            f"A {self.dim}-dimensional vector configuration consisting of "
            f"#{self.size} vectors"
        )

    def __hash__(self) -> int:
        """
        **Description:**
        Hash for the vector configuration. Defined by hashing a dictionary from
        labels to vectors.

        **Arguments:**
        None.

        **Returns:**
        The hash.
        """
        # immutable dictionary-like object mapping labels to vectors
        l2v_immut = [
            (label, tuple(self.vector(label))) for label in sorted(self._labels)
        ]
        l2v_immut = tuple(l2v_immut)

        return hash(l2v_immut)

    def __ne__(self, o: "VectorConfiguration") -> bool:
        """
        **Description:**
        Inequality checking between two VectorConfiguration objects.

        **Arguments:**
        - `o`: The other VectorConfiguration to compare against.

        **Returns:**
        True if self!=o. False if self==o.
        """
        # check type
        if (self.__class__.__name__   != o.__class__.__name__):
            return False
        if (self.__class__.__module__ != o.__class__.__module__):
            return False

        # check that labels and vectors identically match
        if self.labels != o.labels:
            return True
        elif (self.vectors() != o.vectors()).any():
            return True

        # all checks passed
        return False

    def __eq__(self, o: "VectorConfiguration") -> bool:
        """
        **Description:**
        Equality checking between two VectorConfiguration objects.

        **Arguments:**
        - `o`: The other VectorConfiguration to compare against.

        **Returns:**
        True if self==o. False if self!=o.
        """
        return not self.__ne__(o)

    def copy(self) -> "VectorConfiguration":
        """
        **Description:**
        Copy method.

        **Arguments:**
        None.

        **Returns:**
        A copy of the vector configuration.
        """
        return copy.deepcopy(self)

    # getters
    # =======
    @property
    def labels(self) -> tuple[int]:
        """
        **Description:**
        Returns the labels of the vectors in the VC.

        **Arguments:**
        None.

        **Returns:**
        The labels of the vectors in the VC.
        """
        return self._labels

    @property
    def labels_to_inds_dict(self) -> dict[int, int]:
        """
        **Description:**
        Returns the a dictionary mapping vector labels to their indices in the
        vector configuration.

        **Arguments:**
        None.

        **Returns:**
        The mapping from labels to indices.
        """
        if self._labels_to_inds is None:
            self._labels_to_inds = {
                label:ind for ind,label in enumerate(self.labels)
            }

        return self._labels_to_inds

    @property
    def size(self) -> int:
        """
        **Description:**
        Returns the number of the vectors in the VC.

        **Arguments:**
        None.

        **Returns:**
        The number of the vectors in the VC.
        """
        return self._vectors.shape[0]

    @property
    def ambient_dim(self) -> int:
        """
        **Description:**
        Returns the ambient dimension of the VC.

        **Arguments:**
        None.

        **Returns:**
        The ambient dimension of the VC.
        """
        return self._vectors.shape[1]

    @property
    def dim(self) -> int:
        """
        **Description:**
        Returns the dimension of the VC. I.e., the dimension of the subspace
        spanned by the vectors.

        **Arguments:**
        None.

        **Returns:**
        The dimension of the VC.
        """
        if self._dim is None:
            self._dim = util.cone_dim(R=self.vectors())

        return self._dim

    # vectors
    # -------
    def vectors(self, which: int | Iterable[int] = None) -> "ArrayLike":
        """
        **Description:**
        Returns the vectors, optionally only those with given labels.

        **Arguments:**
        - `which`: Either a single label, for which the single corresponding
                   vector will be returned, or a list of labels. If not
                   provided, then all vectors are returned.

        **Returns:**
        The corresponding vector(s), in order specified by which.
        """
        # if no labels are provided, return all vectors
        if which is None:
            which = self.labels

        # cast to iterable
        single_vec = not isinstance(which, Iterable)
        if single_vec:
            which = (which,)

        # return
        out = np.array([self._labels_to_vectors[label] for label in which])
        if single_vec:
            out = out[0]
        return out

    # aliases
    vector = vectors

    def vectors_to_labels(self, vectors: "ArrayLike") -> int | list[int]:
        """
        **Description:**
        Maps the vectors to their corresponding labels

        **Arguments:**
        - `vectors`: Either a single vector, for which the single corresponding
                     label will be returned, or a list of vectors.

        **Returns:**
        The corresponding label(s).
        """
        # input sanitization
        vectors = np.array(vectors)

        # ensure that vectors is a 2D array
        if len(vectors.shape) == 1:
            vectors = vectors.reshape(1, -1)
            return_list = False
        else:
            return_list = True

        # map vectors to labels
        out = [self._vectors_to_labels.get(tuple(v), None) for v in vectors]

        # return
        if not return_list:
            out = out[0]
        return out

    def labels_to_inds(self,
                      labels: Iterable[int],
                      ambient_labels: Iterable[int] = None,
                      offset: int = 0) -> int | Iterable[int]:
        """
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
        """
        # optimization for standard labels 1, ..., N
        if (ambient_labels is None) and self._standard_labels:
            if not isinstance(labels, Iterable):
                return labels - 1 + offset
            else:
                return tuple([i - 1 + offset for i in labels])

        # get default labels
        # construct dict mapping label to index
        if ambient_labels is None:
            _labels_to_inds = self.labels_to_inds
        else:
            _labels_to_inds = {label:ind for ind,label in enumerate(ambient_labels)}

        # either return a single index, or a tuple of indices
        if not isinstance(labels, Iterable):
            return _labels_to_inds[labels] + offset
        else:
            return tuple([_labels_to_inds[i] + offset for i in labels])

    # aliases
    label_to_ind = labels_to_inds

    # basic properties
    # ================
    def is_solid(self) -> bool:
        """
        **Description:**
        Return whether or not the VC is full-dimensional.

        **Arguments:**
        None.

        **Returns:**
        True if the VC is full-dimensional. False otherwise.
        """
        return self.ambient_dim == self.dim

    # aliases
    is_full_dim = is_solid

    def is_totally_cyclic(self) -> bool:
        """
        **Description:**
        Return whether or not the VC is totally cyclic. That is, whether
        self.conv() equals the subspace containing it (the supporting
        hyperplane).

        **Arguments:**
        None.

        **Returns:**
        True if the VC is totally cyclic. False otherwise.
        """
        if not self.is_solid():
            # could definitely be generalized to non-solid
            # likely just check if
            # len(dual_cone(self.vectors())) == 2*(ambient-dim)
            raise NotImplementedError

        return len(util.dual_cone(self.vectors())) == 0

    def is_acyclic(self) -> bool:
        """
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
        """
        return util.is_solid(H=self.vectors())

    def support(self) -> "ArrayLike":
        """
        **Description:**
        Get the support of the vector configuration as a hyperplane
        representation.

        **Arguments:**
        None.

        **Returns:**
        The hyperplanes defining the support.
        """
        return util.dual_cone(self.vectors())

    # cones
    # -----
    # ray containment
    def cone_contains(self,
                      cone_labels: Iterable[int],
                      vec_label: Iterable[int],
                      strict: bool = False) -> bool:
        """
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
        """
        # combine all of the labels
        labs = list(cone_labels)+[vec_label]

        # get the circuit in the VC
        # this is a bit of a misnomer... here, we are really just computing a
        # dependency (not necessary a circuit)
        #
        # for use as in lemma 4.1.11 of DRS
        circ = self.circuit(labs)

        # not even a dependency - cone definitely doesn't contain the vec
        if circ is None:
            return False

        # check if this circuit is insertion/deletion
        circ_type = circ.signature
        if circ_type[0] == 1:
            Zsmall = circ.Zpos
            Zlarge = circ.Zneg
        elif circ_type[1] == 1:
            Zsmall = circ.Zneg
            Zlarge = circ.Zpos
        else:
            # not insertion/deletion
            return False

        # check if the vector is what's being inserted/deleted
        if vec_label != Zsmall[0]:
            return False

        # check if the containment is strict
        if strict:
            return set(Zlarge) == set(cone_labels)

        # all checks passed... return True
        return True

    # regularity
    # ==========
    def gale(self, set_basis: bool = False) -> "ArrayLike":
        """
        **Description:**
        Compute the gale transform of the config.

        I.e., a basis of the null-space of the vectors.

        **Arguments:**
        - `set_basis`: Whether to set a particular basis of the Gale transform.

        **Returns:**
        The gale transform.
        """
        # compute it
        if set_basis:
            assert self._gale_basis is not None

        if set_basis and (self._gale_in_basis is not None):
            return self._gale_in_basis
        elif (not set_basis) and (self._gale is not None):
            return self._gale
        
        # compute null space
        A = self.vectors().T.tolist()
        B, nullity = flint.fmpz_mat(A).nullspace()

        # map to a numpy array
        B = np.array(B.tolist()).astype(int)
        B = B.T[:nullity]
        B = B//np.gcd.reduce(B,axis=1).reshape(-1,1)

        if set_basis:
            # change basis
            P = np.zeros(shape=B.shape, dtype=int)

            gale_basis_inds = self.labels_to_inds(self._gale_basis)
            for i,j in enumerate(gale_basis_inds):
                P[i,j] = 1
        
            C = np.linalg.inv(P @ B.T)
            B = (B.T @ C).T
        
            # map back to integral
            B = np.rint(B).astype(int)

        # save/return
        if set_basis:
            self._gale_in_basis = B.T
            return self._gale_in_basis
        else:
            self._gale = B.T
            return self._gale

    def project(self, vec: "ArrayLike") -> "ArrayLike":
        """
        **Description:**
        Project down a vector from height-space to chamber-space.

        **Arguments:**
        - `vec`: The height-space vector.

        **Returns:**
        The chamber-space vector.
        """
        return self.gale().T@vec

    # aliases
    proj = project

    def jorp(self, vec: "ArrayLike") -> "ArrayLike":
        """
        **Description:**
        Undo a projection from height-space to chamber-space.

        I.e., map from chamber-space to height-space

        **Arguments:**
        - `vec`: The chamber-space vector.

        **Returns:**
        The chamber-space vector.
        """
        return np.linalg.lstsq(self.gale().T, vec, rcond=None)[0]

    # generating triangulations
    # =========================
    def subdivide(
        self,
        heights: "ArrayLike" = None,
        cells: "ArrayLike" = None,
        tol: float = 1e-14,
        seed: int = 0,
        verbosity: int = 0,
    ) -> "Fan":
        """
        **Description:**
        Subdivide the vector configuration either by specified cells/simplices
        or by heights.

        **Arguments:**
        - `heights`:   The heights to lift the vectors by.
        - `cells`:     The cells to use in the triangulation.
        - `backend`:   The lifting backend. Use 'qhull'.
        - `tol`:       Numerical tolerance used.
        - `verbosity`: The verbosity level. Higher is more verbose

        **Returns:**
        The resultant subdivision.
        """
        # triangulate via cells
        # =====================
        if cells is not None:
            return fan.Fan(self,
                           cones=cells,
                           heights=heights)

        # triangulate via heights
        # =======================
        # (if no heights are provided, compute Delaunay triangulation)
        # (need to add noise to ensure it is a *triangulation* and not a
        #  subdivision)
        # (allow retrying in case the noise brings the heights outside the
        #  secondary fan... exceedingly unlikely though)
        if heights is None:
            heights = np.sum(self.vectors()*self.vectors(), axis=1)

        # ensure the heights are non-negative
        already_nonneg = all([h_i >= 0 for h_i in heights])
        if not already_nonneg:
            if verbosity >= 1:
                print("Heights must be cured from negative components...")

            # more detailed check... see DRS 4.1.39
            B = self.gale(set_basis=False)
            Bh = B.T@heights
            heights_new, res = sp.optimize.nnls(B.T, Bh)
            if res>tol:
                print(f"Residuals {res} > tol {tol}...")
                raise ValueError

            # do the check
            if False: # too slow... instead just check if we found coeffs...
                H = sigma.hyperplanes()

                if np.any(H@Bh < 0):
                    msg =   "Heights outside support of secondary fan! "
                    msg += f"{H@B.T@heights}"
                    raise ValueError(msg)

            # get non-negative heights
            if verbosity >= 1:
                msg =   "Check coeffs: B.T@heights_new == Bh?"
                msg += f"{np.allclose(B.T@heights_new, Bh)}"
                print(msg)
                
            if verbosity >= 3:
                # check that heights differ by linear evaluation of vectors
                c, res, *_ = np.linalg.lstsq(self.vectors(),heights-heights_new)
                print('differ by linear eval of A?', max(res)<tol)

            if heights_new is None:
                raise ValueError(f"Heights outside support of secondary fan!")
            heights = heights_new

        # lift & compute simplices
        # ------------------------
        # check for heights=0
        if np.max(heights)==0:
            return self.subdivide(cells=[self.labels])

        # nonzero heights -> lift via a point configuration
        pts = np.vstack( [np.zeros((1,self.dim), dtype=int), self.vectors()] )
        pc  = triangulumancer.PointConfiguration(pts)

        # adjust heights for PC such that the triangulation is star...
        # just ensure that 0 is in all simplices
        #
        # this should always be true by construction. Maybe perturbations of
        # the heights in the various backend for odd heights like those which
        # give subdivisions cause the origin to be skipped...
        height_norm = np.linalg.norm(heights)
        height_orig = -height_norm
        while True:
            heights_pc  = np.concatenate(([height_orig], heights))
            simp_pcinds = pc.triangulate_with_heights(heights_pc).simplices()

            # lower the height of the origin if not star
            if not all([0 in simp for simp in simp_pcinds]):
                height_orig -= height_norm
                continue

            # star :)
            break

        # read the simplices as labels
        simp_vcinds = [[pti-1 for pti in s if pti!=0] for s in simp_pcinds]
        simp_labels = [[self.labels[vci] for vci in s] for s in simp_vcinds]

        return self.subdivide(cells=simp_labels)

    # aliases
    triangulate = subdivide

    def all_triangulations(
        self,
        only_fine: bool = False,
        only_regular: bool = True,
        verbosity: int = 0
    ) -> list["Fan"]:
        """
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
        """
        G, triangs, labs = self.flip_graph(
            only_fine=only_fine, only_regular=only_regular, verbosity=verbosity
        )

        return triangs

    def random_triangulations_fast(
        self,
        method: str="delaunay",
        h0: "ArrayLike" = None,
        sigma: float = 0.1,  # for delaunay
        N: int = None,
        as_list: bool = False,
        attempts_per_triang: int = 1000,
        backend: str = "qhull",
        verbosity: int = 0,
    ) -> Generator["Fan"] | list["Fan"]:
        """
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
        - `backend`:            The lifting backend.
        - `verbosity`:          The verbosity level. Higher is more verbose.

        **Returns:**
        The random triangulations.
        """
        # set default height
        if method == "delaunay":
            if h0 is None:
                h0 = np.sum(self.vectors()*self.vectors(), axis=1)
        elif method == "isotropic":
            if not hasattr(self, "_vector_norms"):
                self._vector_norms = np.linalg.norm(self.vectors(), axis=1)
        else:
            raise ValueError(f"Unrecognized method = '{method}'")

        if as_list:
            # get the generator
            gen = self.random_triangulations_fast(  # high=high,
                h0=h0,
                sigma=sigma,
                N=None,
                as_list=False,
                backend=backend,
                verbosity=verbosity,
            )

            # main object of interest
            triangs = set()

            # fill until done
            num_Ts = 0
            while num_Ts < N:
                for _ in range(attempts_per_triang):
                    if verbosity >= 1:
                        print(
                            f"Constructing triangulation #{num_Ts} "
                            f"(out of {N})... "
                            f"(attempt #{_} for this triangulation)",
                            end="\r",
                        )
                    # try generating a new triangulation...
                    triangs.add(next(gen))

                    if len(triangs) > num_Ts:
                        # actually new!
                        num_Ts += 1
                        break
                else:
                    # hit limit on attempts/triang... quitting!
                    return list(triangs)

            return list(triangs)

        def gen():
            # define iterator that can handle infinite looping (if N is None)
            if N is None:
                iterator = iter(int, 1)
            else:
                iterator = range(N)

            # generate the triangulations
            for _ in iterator:
                if method == "delaunay":
                    # generate triangulations near Delaunay
                    while True:
                        h = h0 + np.random.normal(scale=sigma, size=len(h0))
                        if all(h >= 0):
                            # valid heights
                            break
                elif method == "isotropic":
                    # pick random heights with non-negative components
                    h = np.random.normal(size=self.size)
                    h = np.multiply(h, np.sign(h))

                    # multiply by vector norms
                    # (think: vector norms are meaningless for VC... these
                    #  heights make the most sense when all vectors are unit
                    #  norm... just scale accordingly)
                    h = np.multiply(h, self._vector_norms)

                try:
                    t = self.subdivide(heights=h, backend=backend)
                except sp.spatial.qhull.QhullError:
                    # QHull error :(
                    if verbosity >= 0:
                        print(f"QHull error for heights = {h}... :( skipping!")
                    continue

                if t.is_triangulation():
                    if verbosity >= 1:
                        print(f"Yielding triangulation via heights = {h}!")
                    yield t

        return gen()

    # flips
    # -----
    def circuit(self,
                labels: Iterable[int],
                lmbda: Iterable = None,
                set_non_dependencies: bool=True,
                save_circuits: bool=True) -> "Circuit":
        """
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
        """
        labels = tuple(sorted(labels))

        # return the answer if known
        circ = self._circuits[labels]
        if circ not in (0, -1):
            # this is the circuit!
            return circ

        # if no dependency is given, check that labels define a circuit
        if lmbda is None:
            # check that this is actually a circuit
            dim = np.linalg.matrix_rank(self.vectors(labels))
            if dim != (len(labels) - 1):
                if set_non_dependencies:
                    self._circuits.set_non_dependency(labels)
                return None

            # compute the dependence
            A = self.vectors(labels).T.tolist()
            X, nullity = flint.fmpz_mat(A).nullspace()
            assert nullity == 1
            lmbda = np.array([int(X[i, 0]) for i in range(X.nrows())])
            lmbda = lmbda//np.gcd.reduce(lmbda)

        # else check the data type
        elif lmbda.dtype != int:
            # dependencies must be integral
            raise ValueError()

        lmbda = tuple(lmbda.tolist())

        # split labels by sign (discard 0s...)
        rel_labels     = []
        rel_dependence = []

        Zpos, Zneg = [], []
        for label, coeff in zip(labels, lmbda):
            if coeff > 0:
                Zpos.append(label)
            elif coeff < 0:
                Zneg.append(label)
            else:
                # l==0... skip!
                continue

            # save the relevant label, dependence
            rel_labels.append(label)
            rel_dependence.append(coeff)

        # reorient if |Zpos| < |Zneg|
        if len(Zpos) < len(Zneg):
            rel_dependence = tuple([-coeff for coeff in rel_dependence])
            Zpos, Zneg = Zneg, Zpos

        # get the type
        Ztype = [
            sum(coeff > 0 for coeff in rel_dependence),
            sum(coeff < 0 for coeff in rel_dependence),
        ]

        # save, return the circuit
        circ = circuits.Circuit(self,
                                Z=tuple(rel_labels),
                                Zpos=tuple(Zpos),
                                Zneg=tuple(Zneg),
                                lmbda=tuple(rel_dependence),
                                signature=tuple(Ztype))
        if save_circuits:
            self._circuits.set_circuit(circ)

        return circ

    def circuits(self, verbosity: int = 0) -> list["Circuit"]:
        """
        **Description:**
        Compute all possible circuits of this vector configuration.

        **Arguments:**
        - `verbosity`: The verbosity level. Higher is more verbose.

        **Returns:**
        A list of Circuit objects.
        """
        # return answer if known
        if self._computed_all_circuits:
            # maybe we should copy? Is mutability a concern?
            return list(self._circuits.circuits.values())

        # calculate the answer
        for npts in range(2, self.dim + 2):
            if verbosity >= 1:
                print(f"Trying to find circuits with #{npts} points...")

            # iterate over all subsets
            for subconfig in itertools.combinations(self.labels, r=npts):
                if verbosity >= 2:
                    print(f"Checking if {subconfig} is a new circuit... ",
                                                                        end="")

                # check if we already contained relevant part of this circuit
                is_known = (self._circuits[subconfig] != 0)
                if is_known:
                    continue

                # compute, save the circuit
                circ = self.circuit(subconfig, set_non_dependencies=True)

        # return
        self._computed_all_circuits = True
        self._circuits.know_all_circuits = True
        return self.circuits(verbosity=0)

    def flip_graph(
        self,
        max_flips: int = None,
        only_fine: bool = False,
        only_regular: bool = True,
        only_pc_triang: bool = False,
        compute_node_labels: bool = False,
        verbosity: int = 0,
    ) -> (nx.Graph, list["Fan"], list[str]):
        """
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
        """
        # lazily compute the flip graph
        args = (
            max_flips,
            only_fine,
            only_regular,
            only_pc_triang,
            compute_node_labels,
        )

        if args not in self._flip_graphs:
            self._flip_graphs[args] = fan.flip_subgraph(
                self,
                max_flips=max_flips,
                only_fine=only_fine,
                only_regular=only_regular,
                only_pc_triang=only_pc_triang,
                compute_node_labels=compute_node_labels,
                verbosity=verbosity,
            )

        # return the output
        return [copy.copy(x) for x in self._flip_graphs[args]]

    def secondary_fan(self,
                      only_fine: bool=False,
                      project_lineality: bool=False,
                      formal_fan: bool=False,
                      verbosity: int=0):
        """
        **Description:**
        Compute the secondary fan of the vector configuration.

        **Arguments:**
        - `only_fine`:         Restrict to fine triangulations.
        - `project_lineality`: Project out lineality space with the gale
                               transform, mapping to the chamber complex.
        - `formal_fan`:        Save as a formal Fan object.
        - `verbosity`:         The verbosity level. Higher is more verbose

        **Returns:**
        The secondary fan triangulations.
        """
        # want the entire fan
        triangs = self.all_triangulations(only_regular=True,
                                          only_fine=only_fine,
                                          verbosity=verbosity)

        # compute all of the secondary cones
        fan   = [fan.secondary_cone(project_lineality=project_lineality) for \
                                                                fan in triangs]

        # map to a formal fan
        if formal_fan:
            rays = np.array(list(
                {tuple(r) for cone in fan for r in cone.rays()}
            ))
            vc = VectorConfiguration(rays)

            cones_as_labels = sorted([sorted(vc.vectors_to_labels(cone.rays()))\
                                                            for cone in fan])
            fan = vc.subdivide(cells=cones_as_labels)

        return fan, triangs

    # misc
    # ----
    def central_fan(self) -> "Fan":
        """
        **Description:**
        Generate the central fan of the vector configuration. Can be defined
        as lifting each vector by a height of 1.

        **Arguments:**
        None.

        **Returns:**
        The central fan.
        """
        return self.subdivide(heights=[1 for _ in self.labels])
