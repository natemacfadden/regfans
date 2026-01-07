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
#               on triangulations of vector configurations (i.e., fans).
# -----------------------------------------------------------------------------

# external imports
from collections.abc import Callable, Iterable
import flint
import itertools
import networkx as nx
import numpy as np
import scipy as sp
from typing import Union
import warnings

# local imports
from . import util, circuits


class Fan:
    """
    This class handles definition/operations on fans. It is analogous to
    CYTools' Triangulation class.

    **Description:**
    Constructs a `Fan` object describing a triangulation of a lattice vector
    configuration. This is handled by the hidden [`__init__`](#__init__)
    function.

    This class is *not* intended to be called directly. Instead, it is meant to
    be called through VectorConfiguration.subdivide.

    **Arguments:**
    - `vc`:      The ambient vector configuration that this fan is over.
    - `cones`:   The cones defining the fan. Each cone is a collection of
                 integer labels.
    - `heights`: The heights defining the fan, if it is regular. Can be
                 computed later.

    **Returns:**
    Nothing.
    """
    def __init__(self,
        vc: "VectorConfiguration",
        cones: list[list[int]],
        heights: list[float]=None):
        """
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
        """
        # read inputs
        # -----------
        self._vc = vc

        N_input_cones = len(cones)
        self._cones = {tuple(sorted([int(l) for l in c])) for c in cones}
        self._cones = tuple(sorted(self._cones))
        if len(self._cones) < N_input_cones:
            msg = "Fan: Input `cones` had duplicates..."
            
            if heights is not None:
                msg += f" heights={np.array(heights).tolist()}..."
            else:
                msg +=  " heights=None..."
            msg += " trimming..."
            warnings.warn(msg)
            
            assert self.is_valid()

        if (heights is None) or any([hi<0 for hi in heights]):
            self._heights = None
        else:
            self._heights = np.array(heights, copy=True)

        # initialize other attributes
        # ---------------------------
        self._used_labels = None
        self._labels_to_cones = None

        self._is_regular = None

        self._circuits = circuits.Circuits()
        self._computed_all_circuits = False

    # defaults
    # ========
    def __repr__(self) -> str:
        """
        **Description:**
        String representation of the Fan.
        (more detailed than __str__)

        **Arguments:**
        None.

        **Returns:**
        String representation of the object.
        """
        fine_str = "fine" if self.is_fine() else "non-fine"

        regular_str = ""
        if self._is_regular is not None:
            regular_str = ", "
            regular_str += "regular" if self._is_regular else "irregular"

        if self.is_triang():
            subdivision_str = "triangulation"
        else:
            subdivision_str = "subdivision"

        return (
            f"A "
            + fine_str
            + regular_str
            + f" {subdivision_str} of "
            + repr(self._vc)
        )

    def __str__(self) -> str:
        """
        **Description:**
        String description of the Fan.
        (less detailed than __repr__ but more readable)

        **Arguments:**
        None.

        **Returns:**
        String description of the object.
        """
        fine_str = "fine" if self.is_fine() else "non-fine"

        regular_str = ""
        if self._is_regular is not None:
            regular_str = ", "
            regular_str += "regular" if self._is_regular else "irregular"

        if self.is_triang():
            subdivision_str = "triangulation"
        else:
            subdivision_str = "subdivision"

        return (
            f"A "
            + fine_str
            + regular_str
            + f" {subdivision_str} of "
            + str(self._vc)
        )

    def __hash__(self) -> int:
        """
        **Description:**
        Hash for the fan. Defined by hashing vector configuration and the cones.

        **Arguments:**
        None.

        **Returns:**
        The hash.
        """
        return hash((hash(self.vc), tuple(sorted(self._cones))))

    def __eq__(self, o: "Fan") -> bool:
        """
        **Description:**
        Equality checking between two Fan objects.

        **Arguments:**
        - `o`: The other Fan to compare against.

        **Returns:**
        True if self==o. False if self!=o.
        """
        return (self.vc == o.vc) and set(self._cones) == set(o._cones)

    def __ne__(self, o: "Fan") -> bool:
        """
        **Description:**
        Inequality checking between two Fan objects.

        **Arguments:**
        - `o`: The other Fan to compare against.

        **Returns:**
        True if self!=o. False if self==o.
        """
        return not (self==o)

    # getters
    # =======
    @property
    def vector_config(self) -> "VectorConfiguration":
        """
        **Description:**
        Returns the associated vector configuration.

        **Arguments:**
        None.

        **Returns:**
        The associated vector configuration.
        """
        return self._vc

    # aliases
    vc = vector_config

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
        return self.vc.labels

    @property
    def used_labels(self) -> tuple[int]:
        """
        **Description:**
        Returns the labels of the vectors in the VC used by cones in the Fan.

        **Arguments:**
        None.

        **Returns:**
        The labels of the vectors in the VC used by cones in the Fan.
        """
        if self._used_labels is None:
            self._used_labels = [l for simp in self._cones for l in simp]
            self._used_labels = tuple(sorted(set(self._used_labels)))
        return self._used_labels

    @property
    def labels_to_cones(self) -> dict[ int, set[tuple[int]] ]:
        """
        **Description:**
        Returns a map from vector labels to the cones the vector appears in.

        **Arguments:**
        None.

        **Returns:**
        A map from vector label to a set of cones (as tuples of indices) that
        the vector appears in.
        """
        # lazily compute
        if self._labels_to_cones is None:
            self._labels_to_cones = {l:set() for l in self.labels}
            for c in self._cones:
                for l in c:
                    self._labels_to_cones[l].add(c)

        # return
        return self._labels_to_cones

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
        return self.vc.ambient_dim

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
        return self.vc.dim

    # less-trivial getters
    # --------------------
    def vectors(self,
        which: int | Iterable[int] = None,
        lifted: bool = False) -> "ArrayLike":
        """
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
        """
        # default labels
        if which is None:
            which = self.used_labels

        # get the optinally-lifted vectors
        vecs = self.vc.vectors(which=which)
        if lifted:
            inds = list(self.labels_to_inds(which))
            relevant_heights = self.heights()[inds]
            return np.hstack((vecs, relevant_heights.reshape(-1,1)))
        else:
            return vecs

    def cones(self,
        as_rays: bool = False,
        as_hyps: bool = False,
        as_inds: bool = False,
        ind_offset: int=0) -> Union[ tuple[tuple[int]], list["ArrayLike"] ]:
        """
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
        """
        # check that at most one of the format flags is set
        if (as_inds + as_rays + as_hyps) > 1:
            msg = "At most 1 of `as_inds`, `as_rays`, and `as_hyps` can be set."
            raise ValueError(msg)

        # format case-by-case
        if as_inds:
            # as indices
            cones = tuple([
                tuple([
                    self.vc.label_to_ind(i)+ind_offset for i in simp
                ])
                for simp in self._cones
            ])
        elif as_rays:
            # as rays
            cones = [self.vectors(which=simp).tolist() for simp in self._cones]
        elif as_hyps:
            # as hyperplanes
            cones = [util.dual_cone(self.vectors(which=simp)) for simp in\
                                                                    self._cones]
        else:
            # as labels
            cones = self._cones

        return cones

    cells = cones
    simps = cones
    simplices = cones

    def facets(self) -> dict[tuple[int], list[tuple[int]]]:
        """
        **Description:**
        Returns the facets of the cones. Save them as a dictionary from facet
        labels to a list of containing cones, stored by their labels.

        **Arguments:**
        None.

        **Returns:**
        A dictionary from facet labels to a list of containing cones.
        """
        if not self.is_triangulation():
            # the following assumes simplicial cones
            raise NotImplementedError

        # compute the facets as a map from facet labels to containing cones
        facets = dict()
        for cone_labels in self.cones():
            # any subset of #(dim-1) rays defines a facet
            for cc in itertools.combinations(cone_labels, r=self.dim - 1):
                # store each facet as a key in a dictionary
                # store containing solid cones as the associated value(s)
                facets[cc] = facets.get(cc, []) + [cone_labels]

        return facets

    # basic properties
    # ================
    def is_valid(self, verbosity: int = 0) -> bool:
        """
        **Description:**
        Return whether or not the cones define a valid polyhedral fan.

        Follows cor. 4.5.13 of "Triangulations" by De Loera, Rambau, Santos.

        **Arguments:**
        - `verbosity`: The verbosity level. Higher is more verbose.

        **Returns:**
        True if the cones define a valid fan. False otherwise.
        """
        # helpers
        # -------
        # get the cones
        _cones = {c:hyps for c,hyps in zip(self.cones(), self.cones(as_hyps=1))}

        # map from facet labels to containing cone labels
        _cone_facets = self.facets()

        # labels of vectors laying in each facet of A
        _A_facets = sorted([
            tuple(sorted([
                l for l,v in self.vc._labels_to_vectors.items()
                if np.dot(n,v)==0
                ]))\
            for n in self.vc.support()])

        # MaxMP
        if verbosity>=1:
            print("Checking MaxMP...")
        for f, containing in _cone_facets.items():
            in_A_facet = any([set(f).issubset(ff) for ff in _A_facets])

            if in_A_facet:
                if len(containing)!=1:
                    if verbosity>=2:
                        print("exterior wall not contained in facet")
                    return False
            elif len(containing) != 2:
                if verbosity>=2:
                    print(f"interior wall {f} contained in {len(containing)}!=2 simplices. {containing}")
                return False

        # MaxAdjHP
        if verbosity>=1:
            print("Checking MaxAdjHP...")
        tmp = (len(_cones) * (len(_cones)-1))//2
        for i, (c1, c2) in enumerate(itertools.combinations(_cones, 2)):
            if verbosity>=2:
                print(f"{i+1}/{tmp}", end='\r')

            # only required for adjacent cells
            label_inter = set(c1).intersection(c2)
            if len(label_inter)==0:
                continue

            # the intersected cone as in an abstract simplicial complex
            R_rhs = self.vectors(which=sorted(label_inter))
            H_rhs = np.array(util.dual_cone(R_rhs))

            # the intersected cone as in a geometric simplicial complex
            H1 = _cones[c1]
            H2 = _cones[c2]
            H_lhs = np.vstack([H1, H2])
            R_lhs = np.array(util.dual_cone(H_lhs))

            # need these two cones to be the same
            if np.all(H_rhs@R_lhs.T >= 0) and np.all(H_lhs@R_rhs.T >= 0):
                # cones are equal
                pass
            else:
                if verbosity>=2:
                    print(f"{c1},{c2} failed")
                return False

        # MaxAdjLP
        if verbosity>=1:
            print("Checking MaxAdjLP...")
        if not self.is_triangulation():
            msg = "MaxAdjLP not yet implemented for subdivisions..."
            raise NotImplementedError(msg)

        # IPP
        if verbosity>=1:
            print("Checking IPP...")
        x = None
        for simp in _cones:
            if x is None:
                x = util.find_interior_point(R = self.vectors(which=simp))
            else:
                # ensure other cones don't include this point
                if util.contains(p=x, R = self.vectors(which=simp)):
                    return False

        return True

    def respects_ptconfig(self) -> bool:
        """
        **Description:**
        Return whether or not the fan also defines a (star) subdivision of the
        underlying point configuration.

        **Arguments:**
        None.

        **Returns:**
        True if the fan defines a subdivision of the point configuration. False
        otherwise.
        """
        if not self.is_regular():
            # could be checked by MaxMP but I have't implemented that...
            raise NotImplementedError

        # just check if central subdivision is a refinement
        H = self.secondary_cone_hyperplanes()
        return all((H @ np.ones(self.vc.size)) >= 0)

    def is_triangulation(self) -> bool:
        """
        **Description:**
        Return whether or not the fan is a triangulation (not a subdivision).

        **Arguments:**
        None.

        **Returns:**
        True if the fan is a triangulation. False otherwise.
        """
        return all([len(c) == self.dim for c in self.cones()])

    is_triang = is_triangulation

    def is_fine(self) -> bool:
        """
        **Description:**
        Return whether or not the fan is fine.

        **Arguments:**
        None.

        **Returns:**
        True if the fan is fine. False otherwise.
        """
        return self.used_labels == self.labels

    def is_regular(self) -> bool:
        """
        **Description:**
        Return whether or not the fan is regular.

        **Arguments:**
        None.

        **Returns:**
        True if the fan is regular. False otherwise.
        """
        if self._is_regular is None:
            H = self.secondary_cone_hyperplanes()
            self._is_regular = util.is_solid(H=H)

        return self._is_regular

    def heights(self) -> list[float] | None:
        """
        **Description:**
        Return some heights defining the cone, if it is regular. Else, return
        None.

        **Arguments:**
        None.

        **Returns:**
        True if the fan is regular. False otherwise.
        """
        if self._heights is not None:
            return self._heights

        if self.is_regular():
            # get a non-negative point since those heights are easier
            H = self.secondary_cone_hyperplanes()
            self._heights = util.find_interior_point(H=H, nonneg=True)
            return self._heights
        else:
            return None

    # cones
    # -----
    def contains(self, c: Iterable[int] | Iterable[Iterable[int]]) -> bool:
        """
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
        """
        # recursively call for all cones
        if isinstance(c[0], Iterable):
            for cc in c:
                if not self.contains(cc):
                    return False
            return True

        # check a single cone
        l2c = self.labels_to_cones
        return len( l2c[c[0]].intersection(*[l2c[l] for l in c[1:]]) )>0

    # flip preliminaries
    # ------------------
    def circuit(
        self,
        labels: Iterable[int] = None,
        enforce_positive: int = None,
        lmbda: Iterable[float] = None,
        check_containment: bool = True,
        save_circuits_in_vc: bool = False,
        verbosity: int = 0,
    ) -> "Circuit":
        """
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
        """
        # check for sufficient inputs
        if lmbda is None:
            assert labels is not None
            assert enforce_positive is not None
        else:
            assert labels is None
            lmbda = np.array(lmbda).reshape(-1)
            try:
                labels = [self.labels[i] for i in np.where(lmbda!=0)[0]]
            except:
                print()
                print('lmbda',lmbda.tolist())
                print('nonzero',np.where(lmbda!=0)[0])
                print('labels',self.labels)
                print('num labels',len(self.labels))
                raise ValueError()
            enforce_positive = self.labels[np.where(lmbda>0)[0][0]]

        # check if it's a known circuits
        if labels in self._circuits:
            return self._circuits[labels]

        # not obviously known... see if we can fetch from VC
        circ = self.vc.circuit(labels,
                               lmbda=lmbda,
                               save_circuits=save_circuits_in_vc)

        # return None if not a circuit
        if circ is None:
            return None

        # got a circuit from VC... embed into Fan
        Z = circ.Z
        if Z in self._circuits:
            return self._circuits[Z]

        if verbosity >= 1:
            print(f"Z = {Z}")

        # re-orient if necessary
        n = circ.lmbda
        for label, coeff in zip(Z, n):
            if label == enforce_positive:
                # enforce `other` has a positive coefficient
                if coeff > 0:
                    pass
                elif coeff < 0:
                    n = tuple([-nn for nn in n])
                else:
                    return None

                break

        if verbosity >= 1:
            print(f"lambda = {n}")

        # compute Z_+, Z_-, Z_0
        Zpos = []
        Zneg = []
        Zzero = []
        for label, coeff in zip(Z, n):
            if coeff > 0:
                Zpos.append(label)
            elif coeff < 0:
                Zneg.append(label)
            else:
                Zzero.append(label)

        # check if circuit is flippable
        if verbosity >= 1:
            flippable = len(Zneg) > 0

            if flippable:
                print("Circuit is flippable!")
            else:
                print("Circuit is NOT flippable...")
                print("(boundary of secondary fan...)")

        # compute T_+, T_-
        Tpos = [tuple(sorted(set(Z) - {x})) for x in Zpos]
        Tneg = [tuple(sorted(set(Z) - {x})) for x in Zneg]

        if verbosity >= 1:
            print(f"T_+ = {Tpos}")
            print(f"T_- = {Tneg}")

        # check that T_+, T_- can be embedded
        # (if so, modify Tpos and Tneg to be the embedded version)
        pos_stars = [self.star(s) for s in Tpos]
        pos_links = [set(self.link(s)) for s in Tpos]
        
        if verbosity >= 1:
            print(f"Tpos has stars = {pos_stars}")
            print(f"     and links = {pos_links}")

        for the_star, the_link in zip(pos_stars, pos_links):
            # check that c is contained in the current fan:
            if len(the_star) == 0:
                if verbosity >= 1:
                    print(f"fan doesn't contain a cell in {Tpos}...")
                    print()
                return None

            # check that the link is constant
            if the_link != pos_links[0]:
                if verbosity >= 1:
                    print("non-constant link...")
                    print(f"links = {pos_links}...")
                    print()
                return None

        # we can flip :) Compute/save the neighbor
        Tpos_embedded = [cc for c in pos_stars for cc in c]
        Tneg_embedded = [cc for c in Tneg for cc in self.embed(c, the_link)]

        if check_containment and (not all([self.contains(c) for c in Tpos])):
            # this can fail
            # E.g., consider cell = {(-1,0,1), (1,0,1), (0,1,1)}
            #       of a 2D fan
            #
            #       consider point p = (0,3,1)
            #
            #       cell+{p} is linearly dependent, and corresponds
            #       to deletion of (0,1,1)
            #
            #       say there is also point q = (0,2,1). Then the
            #       circuit for  cell+{p} is not flippable
            #print('SKIPPING BAD EMBEDDINGS...')
            return None

        # save it
        Ztype = [
            sum(coeff > 0 for coeff in n),
            sum(coeff < 0 for coeff in n),
        ]

        out = circuits.Circuit(self.vc,
                               Z=Z,
                               Zpos=tuple(Zpos),
                               Zneg=tuple(Zneg),
                               lmbda=tuple(n),
                               signature=Ztype)

        # also save Tpos, Tneg
        out.Tpos = Tpos_embedded
        out.Tneg = Tneg_embedded

        # return
        self._circuits.set_circuit(out, verbosity=verbosity-1)
        return out

    def circuits(self,
        facets: dict[Iterable[int], Iterable[Iterable[int]]] = None,
        verbosity: int = 0) -> list["Circuit"]:
        """
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
        """
        # return answer if known
        if self._computed_all_circuits:
            return list(self._circuits.values())

        # output object in case complete_computation=False
        new_circuits = []

        # calculate the answer
        # ====================
        # (can skip combinatorial factors by just iterating over insertion
        # flips and solid cones)

        # setup
        to_insert = set(self.labels).difference(self.used_labels)
        if facets is None:
            facets = self.facets()
            complete_computation = True
        else:
            complete_computation = False

        # insertion/deletion
        # ------------------
        if (verbosity >= 1) and len(to_insert):
            print("Computing the insertion circuits")

        cones = self.cones()
        for missing in to_insert:
            for c in cones:
                # skip if missing isn't in the cone
                if not self.vc.cone_contains(c, missing):
                    continue

                if verbosity >= 2:
                    msg = f"Seeing if we can insert {missing} via {c}... "
                    print(msg, end='')

                # compute the circuit
                circ = self.circuit(c + (missing,),
                                    enforce_positive=missing,
                                    check_containment=False,
                                    verbosity=verbosity-1
                )

                if circ is not None:
                    if verbosity >= 2:
                        print("indeed!")
                        print(f"(circuit = {circ})")
                    break
                else:
                    if verbosity >= 2:
                        print("NOPE")

        # local folding
        # -------------
        if verbosity >= 1:
            print("Computing circuits associated to local folding...")

        for f in facets:
            if verbosity >= 1:
                print(f"Computing circuit associated to facet {f}... ", end="")
                print(f"(containing = {facets[f]}...) ")

            if len(facets[f]) == 1:
                # boundary of the fan. No associated local folding
                if verbosity >= 1:
                    print("only contained in 1x cone...")
                continue
            try:
                c1, c2 = facets[f]
            except:
                raise Exception(f"Facet {f} had !=2 cones {facets[f]}")

            only_c1 = [i for i in c1 if i not in f][0]
            only_c2 = [i for i in c2 if i not in f][0]

            circ = self.circuit(f + (only_c1, only_c2),
                                enforce_positive=only_c1,
                                check_containment=True,
                                verbosity=verbosity-1
            )

            if circ is not None:
                if verbosity >= 1:
                    print(f"-> It's {circ}...")
                new_circuits.append(circ)

        # return
        if complete_computation:
            self._computed_all_circuits = True
            return self.circuits(verbosity=0)
        else:
            return new_circuits

    def star(self,
        cell: Iterable[int],
        old_way: bool = False) -> Iterable[tuple[int]]:
        """
        **Description:**
        Compute the star of some cell. This is the subcomplex of all cones
        containing the cell (and their faces)

        **Arguments:**
        - `cell`:    The cell of interest.
        - `old_way`: Whether to do the computation in an old/slow manner.

        **Returns:**
        A list of all solid cones (as tuples of ints) containing the cell.
        """
        if old_way:
            # old way... bit slow...
            cell = set(cell)
            return [c for c in self.cones() if cell.issubset(c)]
        else:
            l2c = self.labels_to_cones
            return list(l2c[cell[0]].intersection(*[l2c[l] for l in cell[1:]]))

    def link(self, cell: Iterable[int]) -> list[tuple[int]]:
        """
        **Description:**
        Compute the link of some cell. This is the subcomplex of all cones in
        the star that don't intersect the cell.

        **Arguments:**
        - `cell`: The cell of interest.

        **Returns:**
        A list of all solid cones (as tuples of ints) containing the cell.
        """
        return [tuple(sorted(set(c) - set(cell))) for c in self.star(cell)]

    def embed(self,
        cell: Iterable[int],
        link: Iterable[Iterable[int]]) -> list[tuple[int]]:
        """
        **Description:**
        Embed some cell into the Fan bu combining it with each cell in the link.

        **Arguments:**
        - `cell`: The cell of interest.
        - `link`: The link of said cell.

        **Returns:**
        A list of solid cones representing the embedding of the cell into the
        Fan via the link.
        """
        return [tuple(sorted(cell + link_cell)) for link_cell in link]

    # flips
    # -----
    def flip(self,
        circ: "Circuit",
        formal: bool = True,
        verbosity: int = 0) -> Union["Fan", tuple[tuple[int]]]:
        """
        **Description:**
        Make a flip across a circuit.

        **Arguments:**
        - `circ`:      The circuit to flip through.
        - `formal`:    Whether to return a formal Fan (otherwise, just a tuple
                       of cones).
        - `verbosity`: The verbosity level. Higher is more verbose.

        **Returns:**
        The flipped Fan.
        """
        if verbosity >= 1:
            print(f"Flipping circuit {circ}...")
            if verbosity >= 2:
                print("out of all circuits:")
                print(self._circuits)
                print("with associated cone-to-circuit maps")
                print(self._circuits.cone_to_circuit)

        # get the positive, negative refinements
        Tpos = circ.Tpos
        Tneg = circ.Tneg

        # check if the circuit is even flippable
        if len(Tneg) == 0:
            if verbosity >= 1:
                print("not flippable (len(Tneg)=0)...")
            return None

        # compute the neighbor...
        neighb = tuple(
            sorted([c for c in self.cones() if c not in Tpos] + Tneg)
        )

        # return
        if formal:
            # update neighb from being a tuple of cones to a formal Fan object
            neighb = self.vc.subdivide(cells=neighb)

            # pass along circuit info to the neighbor
            neighb._circuits = self._circuits.copy()
            neighb._computed_all_circuits = False

            # copy cone information into neighbor
            neighb._labels_to_cones = {l:cones.copy() for l,cones in \
                                                self.labels_to_cones.items()}

            for c in Tpos:
                for l in c:
                    neighb._labels_to_cones[l].remove(c)
            for c in Tneg:
                for l in c:
                    neighb._labels_to_cones[l].add(c)

            if False:
                # THIS PREMISE SEEMS TO BE WRONG...
                # delete unclear non-circuits
                # (i.e., any which share a point with the any cone in Tpos)
                star_pts = {l for c in Tpos for l in c}
                for l in star_pts:
                    # iterate over non-dependencies that l touches
                    for non_depend in self._circuits.label_to_non_dependency[l]:

                        # iterate over all other points
                        for ll,l_nondepends in self._circuits.label_to_non_dependency.items():
                            # skip self._circuits.label_to_non_dependency[l]
                            if l==ll:
                                continue

                            # remove non_depend from
                            # self._circuits.label_to_non_dependency[ll], if it
                            # exists
                            if non_depend in l_nondepends:
                                msg =  f"removing non-dependency {non_depend} "
                                msg += f"(since {l} in star)"
                                print()
                                l_nondepends.remove(non_depend)

                        # empty the entry for l
                        self._circuits.label_to_non_dependency[l] = set()

            # delete now-irrelevant circuits
            # (i.e., any using a cone in Tpos)
            for c in Tpos:
                # for each deleted cone, find each circuit using it
                for Z in neighb._circuits.cone_to_circuit[c]:
                    # delete this circuit from all other cones
                    for cc in self._circuits[Z].Tpos:
                        if (cc!=c) and (len(neighb._circuits.cone_to_circuit.get(cc, []))!=0):
                            neighb._circuits.cone_to_circuit[cc].remove(Z)

                    # delete this circuit from the main list
                    neighb._circuits.pop(Z, None)
                    self.vc._circuits.pop(Z, None)
                    self.vc._computed_all_circuits = False
                    self.vc._circuits.know_all_circuits = False

                # delete the cone
                # (from cone_to_circuit)
                del neighb._circuits.cone_to_circuit[c]

            # compute new circuits
            facets = dict()
            for c in Tneg:
                for f in itertools.combinations(c, self.ambient_dim-1):
                    facets[f] = facets.get(f,[]) + [c]

            for f in facets:
                f_set = set(f)

                if len(facets[f])==1:
                    for c in neighb.cones():
                        if f_set.issubset(c) and c != facets[f][0]:
                            facets[f].append(c)
            new_circuits = neighb.circuits(facets = facets)
            neighb._computed_all_circuits = True
            neighb.vc._circuits.clear_cache()

            if verbosity >= 1:
                print(f"Neighbor, correspondingly, has circuits:")
                print(neighb._circuits)
                print("and associated cone-to-circuit maps")
                print(neighb._circuits.cone_to_circuit)

        return neighb

    def flip_linear(self,
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
        verbosity: int = 0) -> list[int|Exception, "ArrayLike", "Fan", "ArrayLike", int]:
        """
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
        """
        if check_regularity and (not self.is_regular()):
            raise ValueError("Assumes regular triangulation...")

        # warnings
        if not stop_at_deletion:
            warnings.warn("flip_linear struggles when points are deleted...")

        # initial data
        # ------------
        # get the initial triangulation + heights
        T_curr   = self
        sc_curr  = T_curr.secondary_cone_hyperplanes(via_circuits=True,
                                                     verbosity=-1)
        sc_curr  = np.array(sc_curr)
        
        if h_init is None:
            h_init = util.find_interior_point(H=sc_curr)
            h_curr = h_init
        else:
            h_curr   = np.array(h_init)
            if not util.contains(p=h_init, H=sc_curr):
                raise Exception

        # initialization hooks
        if hook_init is not None:
            hook_init(T_curr)

        # target data
        # -----------
        # split based off of a target point vs. a target direction
        seeking_target = (h_target is not None)

        if seeking_target:
            assert direction is None
            h_target = np.array(h_target)
            direction = h_target-h_curr
        else:
            assert h_target is None
            if verbosity >= 0:
                print("DIRECTIONS ARE HANDLED IN A VERY BAD WAY... ")
                print("h_target = h_init+1000*direction...")
            direction = 1000*np.array(direction)
            h_target = h_init+direction
        direction_norm2 = np.dot(direction,direction)
        
        # flip until we reach the target
        # ------------------------------
        num_flips = 0
        if record_fans:
            history_triangs = [self]
        if record_circs:
            history_circs   = []
        
        # main loop
        status = 1
        while True:
            # stop at a fan respecting the point configuration
            if stop_at_pct:
                is_pct = all(sc_curr @ np.ones(T_curr.size) >= 0)
                if is_pct:
                    break

            # reached our target!
            if seeking_target and util.contains(p=h_target, H=sc_curr):
                h_curr = h_target
                if verbosity >= 1:
                    print(f"Reached target in {num_flips} flip(s) :)")
                break

            # print progress
            progress = np.dot(direction, h_curr-h_init)/direction_norm2
            if verbosity >= 1:
                print(f"progress = {progress:.3f}; num_flips = {num_flips}",\
                                                                        end='')

            # took too many flips.. quit!
            if (max_N_flips is not None) and (num_flips >= max_N_flips):
                break

            # stepped too far
            if progress > 1:
                raise ValueError(f"progress={progress}>1")

            # find the first wall that we hit
            first_hit_ind, first_hit_dist = util.first_hit(h_curr, h_target, sc_curr)
            if first_hit_ind is None:
                msg = "first_hit_ind=None... should've been caught earlier... "
                msg += f"min(H@h_curr)={np.min(sc_curr@h_curr)}; "
                msg += f"min(H@h_target)={np.min(sc_curr@h_target)}...; "
                msg += f"H={sc_curr.tolist()}, h_curr={h_curr.tolist()}, "
                msg += f"h_target={h_target.tolist()}"
                raise ValueError(msg)
        
            # get the corresponding circuit
            try:
                circ = T_curr.circuit(lmbda = sc_curr[first_hit_ind],
                                      verbosity=verbosity-1)
            except:
                print(sc_curr.shape, first_hit_ind)
                print('curr dists', sc_curr@h_curr)
                print('target dists', sc_curr@h_target)
                print('h_target',h_target)
                raise ValueError()

            # check the circuit type
            if 0 in circ.signature:
                if verbosity >= 1:
                    print(f"; {circ.Z} is not flippable..." + 20*" ",flush=True)
                status = Exception("Hit non-flippable wall...")
                h_curr = util.lerp(h_curr, h_target, 0.99*first_hit_dist)

                # check that we can compute a next step
                dists = sc_curr@h_tmp
                n_i   = np.argmin(dists)
                n     = sc_curr[n_i]
                if np.dot(n, h_tmp)>0:
                    h_curr = h_tmp
                break
            if 1 in circ.signature and stop_at_deletion:
                if verbosity >= 1:
                    print(f"; {circ.Z} is deletion..." + 20*" ", flush=True)
                status = Exception("Hit deletion wall...")
                h_tmp = util.lerp(h_curr, h_target, 0.99*first_hit_dist)
                
                # check that we can compute a next step
                dists = sc_curr@h_tmp
                n_i   = np.argmin(dists)
                n     = sc_curr[n_i]
                if np.dot(n, h_tmp)>0:
                    h_curr = h_tmp
                break
            if verbosity >= 1:
                print(f"; flipping {circ.Z}..." + 20*" ", flush=True)

            # flip, update heights
            num_flips += 1

            T_new = T_curr.flip(circ)
            if hook_flip is not None:
                hook_flip(T_curr, T_new, circ)
            
            T_curr    = T_new
            sc_curr   = T_curr.secondary_cone_hyperplanes()
            sc_curr   = np.array(sc_curr)

            # save to history
            if record_circs:
                history_circs.append(circ)
            
            # compute the distance to next hyperplane
            # place self halfway across
            _, next_hit_dist = util.first_hit(h_curr, h_target, sc_curr,
                                         max_dist=float('inf'),
                                         verbosity=0)

            try:
                if (next_hit_dist is None) or (next_hit_dist >= 1):
                    h_curr = h_target
                    assert util.contains(p=h_curr, H=sc_curr)
                else:
                    tmp_dist = 0.5*(first_hit_dist+next_hit_dist)
                    h_curr = util.lerp(h_curr, h_target, tmp_dist)
            except:
                print("FAIL")
                print(f"Outputs of first_hit: {_, next_hit_dist}")
                print("Inputs  to first_hit: " + \
                            f"{h_curr.tolist(), h_target.tolist(), sc_curr}")

                dists = sc_curr@util.lerp(h_curr, h_target, first_hit_dist)
                i = np.argmin(dists)
                print(f"dists[argmin] = {dists[i]}")
                print(f"H[argmin]     = {sc_curr[i].tolist()}")
                print(f"circ          = {circ}")

                sc_hyps = {tuple(n) for n in sc_curr}
                n       = tuple(-sc_curr[i]) 
                print(f"n in sc_curr? = {n in sc_hyps}")
                
                raise ValueError()

            # check that the heights make sense
            if util.contains(p=h_curr, H=sc_curr):
                pass
            else:
                # maybe we landed just outside the SC since h_curr is basically
                # on a wall of the SC try to fix by nudging h_curr back in
                dists = sc_curr@h_curr
                if sum(dists<0) == 1:
                    # assume we have n.h_curr = -eps
                    # update h_curr->h_curr + 2*eps*n/(n.n)
                    # then n.(h_curr + 2*eps*n/(n.n)) = eps
                    v_i = np.where(dists<0)[0][0]
                    n   = sc_curr[v_i]
                    violation = -dists[v_i]

                    if violation > eps:
                        msg = f"Violated hyperplanes by {violation}>{eps}..."
                        raise Exception(msg)

                    dh = 2*violation*n/np.dot(n,n)
                    if np.all(sc_curr@(h_curr+dh)>0):
                        # nudging worked!
                        h_curr = h_curr + dh
                    else:
                        # add a bit of noise and retry
                        msg =  "Looks like we landed on a codim-2+ facet of "
                        msg += "the SC... forcing point to be in SC..."
                        warnings.warn(msg)
                        h_curr  = util.find_interior_point(H=sc_curr)
                        #raise Exception(f"h_curr looks like it landed on a wall of a secondary cone ({eps} off of wall {v_i}) but naive curing didn't fix it... nudged dists = {(sc_curr.hyperplanes()@(h_curr+dh)).tolist()}")
                else:
                    msg =   "Just flipped but new heights not in new secondary "
                    msg += f"cone... min(H@h)={min(dists)}"
                    raise Exception()

            # add to history
            if record_fans:
                history_triangs.append(T_curr)

        # return
        output = [status, h_curr, T_curr, sc_curr, num_flips]

        if record_fans:
            output.append(history_triangs)
        if record_circs:
            output.append(history_circs)
        
        return output

    def neighbors(self,
        only_fine: bool = False,
        formal: bool = True,
        verbosity: int = 0) -> tuple[ list[Union["Fan", tuple[tuple[int]]]], list["Circuit"] ]:
        """
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
        """
        neighbs = []
        neighb_circs = []

        circs = self.circuits(verbosity=verbosity-2)
        for circ in circs:
            # get the positive, negative refinements
            Tpos = circ.Tpos
            Tneg = circ.Tneg

            if verbosity >= 1:
                if verbosity >= 2:
                    print()
                msg =  f"Studying circuit {circ.Z} with Tpos = {Tpos}; "
                msg += f"Tneg={Tneg}..."
                print(msg, end=" ")

            # check if the circuit is even flippable
            if len(Tneg) == 0:
                if verbosity >= 1:
                    msg =  "not flippable (len(Tneg)=0)... bdry of secondary "
                    msg += "fan..."
                    print(msg)
                continue

            # check if the circuit is fine
            if only_fine and (1 in circ.signature):
                if verbosity >= 1:
                    print("insertion/deletion flip")
                continue

            # compute the flip
            neighbs.append(self.flip(circ, formal=formal,verbosity=verbosity-1))
            neighb_circs.append(circ)

            # print info
            if verbosity >= 1:
                print(f"fan contains Tpos ({Tpos}) and is flippable!")
                print(f"(neighbor = {neighbs[-1]}")

            continue

        # return
        return neighbs, neighb_circs

    neighbor_triangulations = neighbors

    def secondary_cone_hyperplanes(self,
        via_circuits: bool = False,
        verbosity: int = 0) -> "ArrayLike":
        """
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
        """
        # compute the normals via circuits
        if via_circuits:
            if verbosity>=0:
                msg =  "The construction with `via_circuits=True` should NOT "
                msg += "be used for checking regularity... it can give rise "
                msg += "to a SOLID cone output for an irregular fan "
                msg += "(presumably if the cone lays fully outside the support "
                msg += "of the secondary fan... i.e., has no valid heights).\n"
                msg += "For checking regularity, use local folding instead "
                msg += "by setting `via_circuits=False`.\n"
                msg += "Disable this warning with `verbosity=-1`."
                print(msg)

            # very simple construction: collect all the (appropriately
            # oriented) circuits as the hyperplanes
            H = []
            for circ in self.circuits(verbosity=verbosity-1):
                n = [0]*self.vc.size
                for l,coeff in zip(circ.Z, circ.lmbda):
                    n[l-1]=coeff
                H.append(n)

        # do so via local folding
        else:
            if verbosity >= 1:
                print("Computing facets...")
            facets = self.facets()

            if verbosity >= 1:
                print("Computing hyperplanes...")

            H = []
            # insertion/deletion
            # ------------------
            for missing in set(self.labels).difference(self.used_labels):
                # check inclusion of l in each cone
                for c in self.cones():
                    circ = self.circuit((missing,) + c, enforce_positive=missing)
                    
                    # check if actually a circuit
                    if circ is None:
                        continue

                    # compute/save normal
                    n = [0] * self.vc.size
                    for i, k in enumerate(circ.Z):
                        n[self.vc.label_to_ind(k)] = circ.lmbda[i]

                    H.append(n)

            # local folding
            # -------------
            for f in facets:
                if verbosity >= 2:
                    print(f"Computing hyperplanes associated to facet {f}...")
                # "MaxMP" - non-bdry facet must be shared by two maximal cells
                if len(facets[f]) != 2:
                    # bdry facet
                    continue

                # define objects s.t. non-shared vectors are separated
                c1, c2 = facets[f]

                only_c1 = [i for i in c1 if i not in f][0]
                only_c2 = [i for i in c2 if i not in f][0]

                circ = (only_c1, only_c2) + f

                if verbosity >= 2:
                    print(f"('circuit' = {circ})")

                # compute the normal
                A = self.vectors(circ).T.tolist()
                X, nullity = flint.fmpz_mat(A).nullspace()
                if nullity != 1:
                    continue
                normal = np.array([int(X[i, 0]) for i in range(X.nrows())])
                normal = normal//np.gcd.reduce(normal)

                # set the sign of the normal
                if normal[0] < 0:
                    normal *= -1

                normal = tuple(normal.tolist())

                if verbosity >= 2:
                    print(f"(n_small = {normal})")

                # compute/save normal
                n = [0] * self.vc.size
                for ind, label in enumerate(circ):
                    n[self.vc.label_to_ind(label)] = normal[ind]

                if verbosity >= 2:
                    print(f"(n_full = {n})")

                H.append(n)

        # return
        return H



# flip graph
# ----------
def flip_subgraph(
    seed,
    max_flips: int = None,
    only_fine: bool = False,
    only_regular: bool = True,
    only_pc_triang: bool = False,
    compute_node_labels: bool = False,
    verbosity: int = 0,
) -> tuple[ "networkx.Graph", list["Fan"], list[dict] ]:
    """
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
    """
    if verbosity >= 1:
        print(f"Computing the flip subgraph centered at {seed}...")
    if verbosity >= 2:
        print(f"(max_flips={max_flips})")
        print(f"(only_fine={only_fine})")
        print(f"(only_regular={only_regular})")
        print(f"(only_pc_triang={only_pc_triang})")

    # ensure that the input is a triangulation
    if isinstance(seed, Fan):
        tri_init = seed
        labels = [
            {
                "reg": tri_init.is_regular(),
                "fine": tri_init.is_fine(),
                "triang": tri_init.respects_ptconfig(),
            }
        ]
    else:
        tri_init = seed.subdivide()
        if not tri_init.is_triangulation():
            print("ERROR: seed isn't triangulation! Quitting!")
            return
        elif only_fine and (not tri_init.is_fine()):
            print("ERROR: seed isn't fine! Quitting!")
            return
        elif only_pc_triang and (not tri_init.respects_ptconfig()):
            print("ERROR: seed isn't PC triang! Quitting!")
            return

        labels = [
            {
                "reg": True,
                "fine": tri_init.is_fine(),
                "triang": True,
            }
        ]

    # hold the subgraph in a custom data structure
    triangs = [tri_init]

    G   = nx.Graph()
    n_i = 0
    val = tri_init.cones()
    fan_to_ind = {val: n_i}
    G.add_node(n_i, value=val)
    n_i += 1

    # for a given stage (e.g., built 2 layers out from the core), we know ALL
    # of the neighbors in the inner layers (e.g., layer 0 and 1). The outermost
    # layer will have more neighbors (unless we already have the entire graph).
    # Define the following variable to count the triangulations in inner layers
    num_fully_checked = 0

    # for each 'layer' in the subgraph (the `onion') of the triangulation graph
    layer = 0
    while (max_flips is None) or (layer < max_flips):
        if num_fully_checked == G.number_of_nodes():
            break  # we've already checked the outer layer for triangulations

        # increment the layer
        layer += 1
        if verbosity >= 1:
            print("Checking layer #", layer, " in the onion...")
            if verbosity >= 2:
                print(f"(this is {triangs[num_fully_checked:]}...)")

        # after adding this layer, how many FRSTs will be fully checked
        num_fully_checked_post_layer = G.number_of_nodes()

        # for every triangulation in the outer layer
        for i in range(num_fully_checked, num_fully_checked_post_layer):
            triang = triangs[i]

            # for each neighbor triangulation
            for neighb, circ in zip(*triang.neighbors(verbosity=verbosity-1)):
                if verbosity >= 2:
                    print(f"Studying neighbor = {neighb}...")
                # see if this is new
                if neighb.cones() in fan_to_ind:
                    G.add_edge(i, fan_to_ind[neighb.cones()], label=circ.data.copy())
                    if verbosity >= 2:
                        print(f"Adding edge {i, fan_to_ind[neighb.cones()]}")
                    continue

                # check various user-imposed restrictions
                new_label = {"reg": None, "fine": None, "triang": None}
                if compute_node_labels:
                    # compute regularity (/maybe check it)
                    # ------------------
                    if verbosity >= 2:
                        print("checking regularity...", end=" ")
                    new_label["reg"] = neighb.is_regular()
                    if verbosity >= 2:
                        print(new_label["reg"])

                    if only_regular and not new_label["reg"]:
                        continue

                    # compute fineness (/maybe check it)
                    # ----------------
                    if verbosity >= 2:
                        print("checking fineness...", end=" ")
                    new_label["fine"] = neighb.is_fine()
                    if verbosity >= 2:
                        print(new_label["fine"])

                    if only_fine and (not new_label["fine"]):
                        continue

                    # compute PC-triang (/maybe check it)
                    # -----------------
                    if verbosity >= 2:
                        print("checking PC triangulation...", end=" ")
                    new_label["triang"] = neighb.respects_ptconfig()
                    if verbosity >= 2:
                        print(new_label["triang"])

                    if only_pc_triang and (not new_label["triang"]):
                        continue
                else:
                    # lazily compute the labels
                    if only_regular:
                        new_label["reg"] = neighb.is_regular()
                        if not new_label["reg"]:
                            continue
                    if only_fine:
                        new_label["fine"] = neighb.is_fine()
                        if not new_label["fine"]:
                            continue
                    if only_pc_triang:
                        new_label["triang"] = neighb.respects_ptconfig()
                        if not new_label["triang"]:
                            continue

                # this is a new triangulation!
                triangs.append(neighb)
                labels.append(new_label)

                val = neighb.cones()
                fan_to_ind[val] = n_i
                G.add_node(n_i, value=val)
                G.add_edge(i, n_i, label=circ.data.copy())
                n_i += 1

        # now, we have fully checked the neighbors of triang
        num_fully_checked = num_fully_checked_post_layer

    return G, triangs, labels
