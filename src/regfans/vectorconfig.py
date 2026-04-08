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
from __future__ import annotations

import copy
import itertools
import warnings
from collections.abc import Generator, Iterable
from typing import TYPE_CHECKING

import flint
import networkx as nx
import numpy as np
import scipy as sp
import triangulumancer

if TYPE_CHECKING:
    from numpy.typing import ArrayLike

    from .circuits import Circuit
    from .fan import Fan

# local imports
from . import circuits, fan, util


class VectorConfiguration:
    """
    This class handles definition/operations on vector configurations. It is
    analogous to CYTools' Polytope class. This object can be triangulated,
    making a simplicial fan.

    Constructs a `VectorConfiguration` object describing a lattice vector
    configuration.

    Parameters
    ----------
    vectors : ArrayLike
        The vectors defining the VC in row format. I.e., vectors[i], the ith row, is treated as the ith vector in the config.
    labels : Iterable[int] | None
        A list of labels for the vectors. Only integral labels are allowed.
    eps : float
        Threshold for checking for non-integral vectors.
    gale_basis : Iterable[int] | None
        An optional basis for the gale transform. If provided, then the gale transform will be put in a basis such that the submatrix given by these labels equals the identity.
    """
    def __init__(
        self,
        vectors: ArrayLike,
        labels: Iterable[int] | None = None,
        eps: float = 1e-4,
        gale_basis: Iterable[int] | None = None,
    ) -> None:
        """
        Initializes a `VectorConfiguration` object.

        Parameters
        ----------
        vectors : ArrayLike
            The vectors defining the VC in row format. I.e., vectors[i], the ith row, is treated as the ith vector in the config.
        labels : Iterable[int] | None, optional
            A list of integer labels for the vectors. Only integral labels are allowed. Defaults to None.
        eps : float, optional
            Threshold for checking for non-integral vectors. Defaults to 0.0001.
        gale_basis : Iterable[int] | None, optional
            An optional basis for the gale transform. If provided, then the gale transform will be put in a basis such that the submatrix given by these labels equals the identity. Defaults to None.
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

        self._labels = tuple(label for label in labels)
        if not all([isinstance(lbl,int) for lbl in self._labels]):
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
        self._refinements = {}

        self._poly = {}

        self._flip_graphs = {}

        # allow setting of a particular basis of the Gale transform
        self._gale_basis    = gale_basis

        # allow caching of the Gale transform
        self._gale_in_basis = None
        self._gale          = None

    # defaults
    # ========
    def __repr__(self) -> str:
        """
        String representation of the VectorConfiguration.
        (more detailed than __str__)

        Returns
        -------
        out : str
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
        String description of the VectorConfiguration.
        (less detailed than __repr__ but more readable)

        Returns
        -------
        out : str
            String description of the object.
        """
        return (
            f"A {self.dim}-dimensional vector configuration consisting of "
            f"#{self.size} vectors"
        )

    def __hash__(self) -> int:
        """
        Hash for the vector configuration. Defined by hashing a dictionary from
        labels to vectors.

        Returns
        -------
        out : int
            The hash.
        """
        # immutable dictionary-like object mapping labels to vectors
        l2v_immut = [
            (label, tuple(self.vector(label))) for label in sorted(self._labels)
        ]
        l2v_immut = tuple(l2v_immut)

        return hash(l2v_immut)

    def __eq__(self, o: VectorConfiguration) -> bool:
        """
        Equality checking between two VectorConfiguration objects.

        Parameters
        ----------
        o : VectorConfiguration
            The other VectorConfiguration to compare against.

        Returns
        -------
        out : bool
            True if self==o. False if self!=o.
        """
        return not self.__ne__(o)

    def __ne__(self, o: VectorConfiguration) -> bool:
        """
        Inequality checking between two VectorConfiguration objects.

        Parameters
        ----------
        o : VectorConfiguration
            The other VectorConfiguration to compare against.

        Returns
        -------
        out : bool
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

    def copy(self) -> VectorConfiguration:
        """
        Copy method.

        Returns
        -------
        out : VectorConfiguration
            A copy of the vector configuration.
        """
        return copy.deepcopy(self)

    # getters
    # =======
    @property
    def labels(self) -> tuple[int]:
        """
        Returns the labels of the vectors in the VC.

        Returns
        -------
        out : tuple[int]
            The labels of the vectors in the VC.
        """
        return self._labels

    @property
    def labels_to_inds_dict(self) -> dict[int, int]:
        """
        Returns a dictionary mapping vector labels to their indices in the
        vector configuration.

        Returns
        -------
        out : dict[int, int]
            The mapping from labels to indices.
        """
        if self._labels_to_inds is None:
            self._labels_to_inds = {
                label: ind for ind, label in enumerate(self.labels)
            }

        return self._labels_to_inds

    @property
    def size(self) -> int:
        """
        Returns the number of the vectors in the VC.

        Returns
        -------
        out : int
            The number of the vectors in the VC.
        """
        return self._vectors.shape[0]

    @property
    def ambient_dim(self) -> int:
        """
        Returns the ambient dimension of the VC.

        Returns
        -------
        out : int
            The ambient dimension of the VC.
        """
        return self._vectors.shape[1]

    @property
    def dim(self) -> int:
        """
        Returns the dimension of the VC. I.e., the dimension of the subspace
        spanned by the vectors.

        Returns
        -------
        out : int
            The dimension of the VC.
        """
        if self._dim is None:
            self._dim = util.cone_dim(R=self.vectors())

        return self._dim

    # vectors
    # -------
    def vectors(self, which: int | Iterable[int] | None = None) -> ArrayLike:
        """
        Returns the vectors, optionally only those with given labels.

        Parameters
        ----------
        which : int | Iterable[int] | None, optional
            Either a single label, for which the single corresponding vector will be returned, or a list of labels. If not provided, then all vectors are returned. Defaults to None.

        Returns
        -------
        out : ArrayLike
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

    def vectors_to_labels(self, vectors: ArrayLike) -> int | list[int]:
        """
        Maps the vectors to their corresponding labels

        Parameters
        ----------
        vectors : ArrayLike
            Either a single vector, for which the single corresponding label will be returned, or a list of vectors.

        Returns
        -------
        out : int | list[int]
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
                      ambient_labels: Iterable[int] | None = None,
                      offset: int = 0) -> int | Iterable[int]:
        """
        Maps the labels to their indices in ambient_labels, optionally with a
        fixed offset.

        Parameters
        ----------
        labels : Iterable[int]
            The labels of interest.
        ambient_labels : Iterable[int] | None, optional
            The ambient labels to get the indices in. If None, use all labels of the VectorConfiguration. Defaults to None.
        offset : int, optional
            Return i+offset for i the index of a label in ambient_labels. Defaults to 0.

        Returns
        -------
        out : int | Iterable[int]
            The indices of the labels.
        """
        # optimization for standard labels 1, ..., N
        if (ambient_labels is None) and self._standard_labels:
            if not isinstance(labels, Iterable):
                return labels - 1 + offset
            else:
                return tuple(i - 1 + offset for i in labels)

        # get default labels
        # construct dict mapping label to index
        if ambient_labels is None:
            _labels_to_inds = self.labels_to_inds_dict
        else:
            _labels_to_inds = {label: ind for ind, label in\
                                                    enumerate(ambient_labels)}

        # either return a single index, or a tuple of indices
        if not isinstance(labels, Iterable):
            return _labels_to_inds[labels] + offset
        else:
            return tuple(_labels_to_inds[i] + offset for i in labels)

    # aliases
    label_to_ind = labels_to_inds

    # basic properties
    # ================
    def is_solid(self) -> bool:
        """
        Return whether or not the VC is full-dimensional.

        Returns
        -------
        out : bool
            True if the VC is full-dimensional. False otherwise.
        """
        return self.ambient_dim == self.dim

    # aliases
    is_full_dim = is_solid

    def is_totally_cyclic(self) -> bool:
        """
        Return whether or not the VC is totally cyclic. That is, whether
        self.conv() equals the subspace containing it (the supporting
        hyperplane).

        Returns
        -------
        out : bool
            True if the VC is totally cyclic. False otherwise.
        """
        if not self.is_solid():
            # could definitely be generalized to non-solid
            # likely just check if
            # len(dual_cone(self.vectors())) == 2*(ambient-dim)
            raise NotImplementedError("Not implemented for non-solid VCs")

        return len(util.dual_cone(self.vectors())) == 0

    def is_acyclic(self) -> bool:
        """
        Return whether or not the VC is acyclic. That is, whether there exists
        some direction psi such that
            psi.vi > 0 for all vi.

        This is equivalent to defining the cone {x: vi.x >= 0} and checking if
        it is full-dimensional.

        Returns
        -------
        out : bool
            True if the VC is acyclic. False otherwise.
        """
        return util.is_solid(H=self.vectors())

    def support(self) -> ArrayLike:
        """
        Get the support of the vector configuration as a hyperplane
        representation.

        Returns
        -------
        out : ArrayLike
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
            Whether to check if the vector is in the strict interior. Defaults to False.

        Returns
        -------
        out : bool
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
    def gale(self, set_basis: bool = False) -> ArrayLike:
        """
        Compute the gale transform of the config.

        I.e., a basis of the null-space of the vectors.

        Parameters
        ----------
        set_basis : bool, optional
            Whether to set a particular basis of the Gale transform. Defaults to False.

        Returns
        -------
        out : ArrayLike
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
        B = B//np.gcd.reduce(B, axis=1).reshape(-1, 1)

        if set_basis:
            # change basis
            P = np.zeros(shape=B.shape, dtype=int)

            gale_basis_inds = self.labels_to_inds(self._gale_basis)
            for i, j in enumerate(gale_basis_inds):
                P[i, j] = 1

            C = np.linalg.inv(P @ B.T)
            B = (B.T @ C).T

            # map back to integral
            Bint = np.rint(B).astype(int)
            if np.allclose(Bint, B):
                B = Bint
            else:
                raise ValueError

        # save/return
        if set_basis:
            self._gale_in_basis = B.T
            return self._gale_in_basis
        else:
            self._gale = B.T
            return self._gale

    def project(self, vec: ArrayLike) -> ArrayLike:
        """
        Project down a vector from height-space to chamber-space.

        Parameters
        ----------
        vec : ArrayLike
            The height-space vector.

        Returns
        -------
        out : ArrayLike
            The chamber-space vector.
        """
        return self.gale().T@vec

    # aliases
    proj = project

    def jorp(self, vec: ArrayLike) -> ArrayLike:
        """
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
        """
        return np.linalg.lstsq(self.gale().T, vec, rcond=None)[0]

    # generating triangulations
    # =========================
    def triangulate(
        self,
        heights: ArrayLike | None = None,
        cells: ArrayLike | None = None,
        tol: float = 1e-14,
        backend: str | None = None,
        make_fine: bool = None,
        check_heights: bool = True,
        cure_heights: bool = True,
        verbosity: int = 0,
    ) -> Fan:
        """
        Subdivide the vector configuration either by specified cells/simplices
        or by heights.

        Parameters
        ----------
        heights : ArrayLike | None, optional
            The heights to lift the vectors by. Defaults to None.
        cells : ArrayLike | None, optional
            The cells to use in the triangulation. Defaults to None.
        tol : float, optional
            Numerical tolerance used for curing negative heights. Defaults to 1e-14.
        backend : str | None, optional
            The lifting backend. Currently allowed to be "cgal" or "ppl". Defaults to None.
        make_fine : bool, optional
            Convert the triangulation to a fine triangulation, if not already fine. Defaults to None.
        check_heights : bool, optional
            Whether to check that the heights land in the secondary cone of the output triangulation. Defaults to True.
        cure_heights : bool, optional
            If the heights do not land in the secondary cone, try to cure them by linear flipping towards the heights. Defaults to True.
        verbosity : int, optional
            The verbosity level. Higher is more verbose. Defaults to 0.

        Returns
        -------
        out : Fan
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
        # flag as to whether the user didn't input heights
        default_heights = (heights is None)

        if make_fine is None:
            make_fine = default_heights

        # check the backend
        if backend is None:
            if default_heights:
                backend = "cgal"
            else:
                backend = "ppl"
        else:
            backend = backend.lower()
            if backend not in ["cgal", "ppl"]:
                raise ValueError(f"Unrecognized backend '{backend}'...")

        # warning in case the user request PPL backend but gave no heights
        if default_heights and (backend != "cgal"):
            msg = "Non-cgal backends are not trustworthy for Delaunay... "
            msg += f"changing from '{backend}' to cgal..."
            with warnings.catch_warnings():
                warnings.simplefilter("always")
                warnings.warn(msg)
            backend = "cgal"

        # allow small perturbations to the height of the origin to enforce that
        # cgal gives a star triangulation
        if backend == "cgal":
            make_cgal_star = default_heights

        # (if no heights are provided, compute Delaunay triangulation)
        # (need to add noise to ensure it is a *triangulation* and not a
        #  subdivision)
        # (allow retrying in case the noise brings the heights outside the
        #  secondary fan... exceedingly unlikely though)
        if default_heights:
            heights = np.sum(self.vectors()*self.vectors(), axis=1)
        else:
            heights = np.array(heights)

        # ensure the heights are non-negative
        already_nonneg = all([h_i >= 0 for h_i in heights])
        if not already_nonneg:
            if verbosity >= 1:
                print("Heights must be cured from negative components...")

            # more detailed check... see DRS 4.1.39
            B = self.gale(set_basis=False)
            Bh = B.T@heights
            heights_new, res = sp.optimize.nnls(B.T, Bh)
            if res > tol:
                print(f"Residuals {res} > tol {tol}...")
                raise ValueError("Invalid heights")

            # do the check
            if False: # too slow... instead just check if we found coeffs...
                H = sigma.hyperplanes()  # noqa: F821

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
                print('differ by linear eval of A?', max(res) < tol)

            if heights_new is None:
                raise ValueError("Heights outside support of secondary fan!")
            heights = heights_new

        # lift & compute simplices
        # ------------------------
        # check for heights=0
        if np.max(heights) == 0:
            return self.subdivide(cells=[self.labels])

        if verbosity >= 1:
            print("Constructing the triangulation via lifting...", flush=True)

        # nonzero heights -> lift via a point configuration
        if backend == "cgal":
            orig = np.zeros((1, self.ambient_dim), dtype=int)
            pts  = np.vstack([orig, self.vectors()])
            pc   = triangulumancer.PointConfiguration(pts)

            if make_cgal_star:
                # adjust heights for PC such that the triangulation is star...
                # just ensure that 0 is in all simplices
                #
                # this should always be true by construction...
                # but CGAL definitely perturbs heights a bit (e.g., lifting by
                # heights = 0 doesn't lead to trivial subdivision)
                #
                # maybe this causes errors leading to non-star triangulations...
                height_scale = max(1, np.min(heights))
                height_orig  = -1e-6*height_scale

                Niter = 0
                while True:
                    if verbosity >= 2:
                        msg  = f"Iteration {Niter}, trying the origin with "
                        msg += f"height {height_orig}..."
                        print(msg, end='\r')
                        Niter += 1
                    h_pc        = np.concatenate(([height_orig], heights))
                    simp_pcinds = pc.triangulate_with_heights(h_pc).simplices

                    # lower the height of the origin if not star
                    if not all([0 in simp for simp in simp_pcinds]):
                        height_orig -= 10*height_scale
                        continue

                    # star :)
                    if verbosity >= 2:
                        print()
                    break

                # check that we didn't lower the height of origin a crazy amount
                if (verbosity >= 1) and (height_orig < -np.min(heights)):
                    msg = "Significantly lowered the height of the origin... "
                    msg += "maybe something went wrong..."
                    with warnings.catch_warnings():
                        warnings.simplefilter("always")
                        warnings.warn(msg)
            else:
                h_pc        = np.concatenate(([0], heights))
                simp_pcinds = pc.triangulate_with_heights(h_pc).simplices

            # read the simplices as indices in the VC
            if make_cgal_star and (not all([0 in s for s in simp_pcinds])):
                msg = "cgal didn't produce a star triangulation... cells = "
                msg += f"{sorted([sorted(s) for s in simp_pcinds.tolist()])} "
                msg += "(0 corresponds to origin). maybe try PPL..."
                raise ValueError(msg)
            else:
                simp_pcinds = [s for s in simp_pcinds if 0 in s]
            simp_vcinds = [[pti-1 for pti in s if pti != 0] for s in simp_pcinds]

        elif backend == "ppl":
            # construct the rays of the lifted cone
            lifted = np.hstack([heights.reshape(-1,1),self.vectors()])
            lifted = np.array([util.primitive(v) for v in lifted]) # as integers

            H      = np.array(util.dual_cone(lifted))
            satd   = H@lifted.T

            # read the simplices as indices in the VC
            simp_vcinds = [np.where(facet == 0)[0].tolist() for facet in satd]

        # read the simplices as labels
        simp_labels = [[self.labels[vci] for vci in s] for s in simp_vcinds]
        simp_labels = [sorted(simp) for simp in simp_labels]

        # construct the fan
        f = self.triangulate(cells=simp_labels)

        # optionally, make the fan fine
        if make_fine and (not f.is_fine()):
            f = fan.make_fine(f)

        # some sanity checks
        if verbosity >= 1:
            print("Doing sanity checks on the triangulation...", flush=True)

        if not f.is_triangulation():
            if verbosity >= 1:
                msg = "Upon lifting, a non-triangulation subdivision was "
                msg += f"output...  (cells = {f.simplices()}) "
                msg += "double check with another backend (PPL is preferable) "
                msg += "OR perturb heights..."

                with warnings.catch_warnings():
                    warnings.simplefilter("always")
                    warnings.warn(msg)
        else:
            # yes a triangulation...
            # verify that the secondary cone contains the heights
            if (not default_heights) and check_heights:
                H = np.array(f.secondary_cone_hyperplanes())
                dists = H@heights

                if np.any((dists) <= 0):
                    if cure_heights:
                        _, __, f, *___ = f.flip_linear(h_target=heights)
                    else:
                        msg = "Heights not contained in secondary cone... "
                        msg += f"distances = {dists}..."
                        raise RuntimeError(msg)

        return f

    # aliases
    subdivide = triangulate

    def all_triangulations(
        self,
        only_fine: bool = False,
        only_regular: bool = True,
        verbosity: int = 0
    ) -> list[Fan]:
        """
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
        """
        G, triangs, labs = self.flip_graph(
            only_fine=only_fine, only_regular=only_regular, verbosity=verbosity
        )

        return triangs

    def random_triangulations_fast(
        self,
        method: str = "delaunay",
        h0: ArrayLike | None = None,
        sigma: float = 0.1,  # for delaunay
        N: int | None = None,
        as_list: bool = False,
        attempts_per_triang: int = 1000,
        backend: str | None = None,
        seed: int = 0,
        verbosity: int = 0,
    ) -> Generator[Fan] | list[Fan]:
        """
        Generate random regular triangulations by picking random heights.

        Parameters
        ----------
        method : str, optional
            Either "delaunay" or "isotropic". The former picks heights around some input height (e.g., the Delaunay heights). The latter picks heights isotropically. Defaults to "delaunay".
        h0 : ArrayLike | None, optional
            The reference heights, for Delaunay method. Defaults to None.
        sigma : float, optional
            How big of a distribution to study around h0. Defaults to 0.1.
        N : int | None, optional
            The number of triangulations to generate. If as_list, then code will keep track of all triangulations, retrying at most attempts_per_triang tries to get a new triangulation until the list has N triangs. O/w, then the first N height vectors are used (regardless of duplicates). Defaults to None.
        as_list : bool, optional
            Whether to return the triangulations as a list, or as a generator. Defaults to False.
        attempts_per_triang : int, optional
            Quit if we can't generate a new triangulation after this many tries. Defaults to 1000.
        backend : str | None, optional
            The lifting backend. See `VectorConfiguration.triangulate`. Defaults to None.
        seed : int, optional
            A random number seed. Defaults to 0.
        verbosity : int, optional
            The verbosity level. Higher is more verbose. Defaults to 0.

        Returns
        -------
        out : Generator[Fan] | list[Fan]
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
            np.random.seed(seed)
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
                    t = self.triangulate(heights=h, backend=backend)
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
                lmbda: Iterable | None = None,
                set_non_dependencies: bool = True,
                save_circuits: bool = True) -> Circuit:
        """
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
            Whether to save circuits... best to keep True for most circumstances. Defaults to True.

        Returns
        -------
        out : Circuit
            Circuit object containing
                - the support of the circuit as property 'Z',
                - the signed circuit as property 'Zpos' and 'Zneg',
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
            rel_dependence = tuple(-coeff for coeff in rel_dependence)
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

    def circuits(self, verbosity: int = 0) -> list[Circuit]:
        """
        Compute all possible circuits of this vector configuration.

        Parameters
        ----------
        verbosity : int, optional
            The verbosity level. Higher is more verbose. Defaults to 0.

        Returns
        -------
        out : list[Circuit]
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
                self.circuit(subconfig, set_non_dependencies=True)

        # return
        self._computed_all_circuits = True
        self._circuits.know_all_circuits = True
        return self.circuits(verbosity=0)

    def flip_graph(
        self,
        max_flips: int | None = None,
        only_fine: bool = False,
        only_regular: bool = True,
        only_pc_triang: bool = False,
        compute_node_labels: bool = False,
        verbosity: int = 0,
    ) -> (nx.Graph, list[Fan], list[dict]):
        """
        Compute the flip graph. Wrapper for flip_subgraph.

        Parameters
        ----------
        max_flips : int | None, optional
            The maximum number of flips to take from the seed. If none is provided, then the entire flip graph is calculated. Defaults to None.
        only_fine : bool, optional
            Whether to only compute fine triangulations. Defaults to False.
        only_regular : bool, optional
            Whether to only compute regular triangulations. Note, we never will see irregular triangulations that are not connected to regular ones. Defaults to True.
        only_pc_triang : bool, optional
            Whether to only compute triangulations that also correspond to star triangulations of the underlying point config. Defaults to False.
        compute_node_labels : bool, optional
            Whether to check whether each node is fine, regular, and a PC triangulation. Defaults to False.
        verbosity : int, optional
            The verbosity level. Higher is more verbose. Defaults to 0.

        Returns
        -------
        out : (nx.Graph, list[Fan], list[dict])
            - The flip graph as a networkx.Graph object.
            - A list of the triangulations
            - A list of the labels for each triangulation (labels are a
              dictionary from the property to a bool)
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
                      only_fine: bool = False,
                      formal_fan: bool = False,
                      verbosity: int = 0):
        """
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
        """
        # want the entire fan
        triangs = self.all_triangulations(only_regular=True,
                                          only_fine=only_fine,
                                          verbosity=verbosity)

        # compute the (hyperplanes of the) secondary cones
        secondary_cones = [t.secondary_cone_hyperplanes() for t in triangs]

        # map to a formal fan
        if formal_fan:
            fan_R    = [util.dual_cone(H) for H in secondary_cones]
            all_rays = np.array(list({
                 {tuple(r) for cone_R in fan_R for r in cone_R}
            }))

            # construct the VC
            vc = VectorConfiguration(all_rays)

            cones_as_labels = sorted([sorted(vc.vectors_to_labels(cone_R))\
                                                        for cone_R in fan_R])
            secondary_cones = vc.subdivide(cells=cones_as_labels)

        return secondary_cones, triangs

    # misc
    # ----
    def central_fan(self) -> Fan:
        """
        Generate the central fan of the vector configuration. Can be defined
        as lifting each vector by a height of 1.

        Returns
        -------
        out : Fan
            The central fan.
        """
        return self.subdivide(
            heights=[1 for _ in self.labels],
            check_heights=False)
