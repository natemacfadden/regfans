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
# Description:  This module contains a class containing information about
#               circuits
# -----------------------------------------------------------------------------

# external imports
from collections.abc import Iterable
from typing import Union

class Circuit():
    """
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
    """
    def __init__(self, vc, Z, Zpos, Zneg, lmbda, signature):
        """
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
        """
        self.vc        = vc
        self.Z         = Z
        self.Zpos      = Zpos
        self.Zneg      = Zneg
        self.lmbda     = lmbda
        self.signature = tuple(signature)

        # positive/negative triangulations
        self.Tpos      = []
        self.Tneg      = []

        self.normal = [0] * self.vc.size
        for i,j in enumerate([self.vc.label_to_ind(z) for z in self.Z]):
            self.normal[j] = self.lmbda[i]
        self.normal = tuple(self.normal)

    def __repr__(self):
        out = f"A circuit with (Z+,Z-)= ({self.Zpos}, {self.Zneg})"
        out += f"; lambda = {self.lmbda}"
        return out

    @property
    def data(self):
        # spoof dictionaries
        return {
            'Z':self.Z,
            'Z+':self.Zpos,
            'Z-':self.Zneg,
            'lambda':self.lmbda,
            'type':self.signature,
            }

    def reorient(self):
        # allow swapping Zpos <-> Zneg
        reoriented = Circuit(
            self.vc,
            self.Z,
            self.Zneg,
            self.Zpos,
            tuple([-x for x in self.lmbda]),
            (self.signature[1], self.signature[0])
            )
        return reoriented

class Circuits():
    """
    This class is a helper data structure to contain the circuits of some
    vector configuration.

    **Description:**
    Constructs a `Circuits` object describing all circuits of some VC. This
    is handled by the hidden [`__init__`](#__init__) function.

    **Arguments:**
    None.

    **Returns:**
    Nothing.
    """
    def __init__(self):
        """
        **Description:**
        Initializes a `Circuits` object.

        **Arguments:**
        None.

        **Returns:**
        Nothing.
        """
        self.clear_cache() # set attributes here

    # clear cache
    # -----------
    def clear_cache(self):
        # main data type - map from the (encoded) unsigned circuit to a Circuit
        # object
        # **Don't directly access... use getters/setters instead...***
        self.circuits = dict()

        # map from cone to the circuits it is involved in
        self.cone_to_circuit = dict()

        # non-dependencies
        self.non_dependencies = set()

        # whether we have computed/saved all circuits
        self.know_all_circuits = False

    # default methods
    # ---------------
    def __repr__(self):
        # hijack `dict` and `Circuit`
        return self.circuits.__repr__()

    def __str__(self):
        # hijack `dict` and `Circuit`
        return self.circuits.__str__()

    def __len__(self):
        # hijack `dict`
        return len(self.circuits)

    def __contains__(self, label_inds: Iterable[int]):
        encoding = self.encode(label_inds)
        return encoding in self.circuits

    def __getitem__(self, label_inds: Iterable[int]) -> Union["Circuit", int]:
        """
        **Description:**
        Get the circuit corresponding to the indicated indices.

        **Arguments:**
        - `label_inds`: The iterable of vector/label indices.

        **Returns:**
        Cases
            - if indices correspond to known circuit -> the `Circuit`
            - if indices correspond to non-circuit   -> -1
            - if indices aren't known                -> 0
        """
        encoding = self.encode(label_inds)
        
        # check if these labels are known to contain a circuit
        if encoding in self.circuits:
            return self.circuits[encoding]

        # if we know all circuits, then indices cannot correspond to a circuit
        if self.know_all_circuits:
            return -1
        
        # if we don't know all circuits, we check against known non-circuits
        for non_dependency in self.non_dependencies:
            if self.is_subset(encoding, non_dependency):
                return -1

        # unclear if this is a circuit
        return 0

    def set_circuit(self,
                    circuit: "Circuit",
                    verbosity: int = 0) -> None:
        """
        **Description:**
        Set the circuit properties corresponding to the indicated indices.

        **Arguments:**
        - `circuit`:   Dict describing the circuit.
        - `verbosity`: The verbosity level.

        **Returns:**
        Nothing.
        """
        encoding = self.encode(circuit.Z)

        # setting a circuit
        self.circuits[encoding] = circuit

        # keep a map from cones to the circuits they have
        for c in circuit.Tpos:
            self.cone_to_circuit[c] = self.cone_to_circuit.get(c,set())
            self.cone_to_circuit[c].add(encoding)

    def set_non_dependency(self,
                           label_inds: Iterable[int],
                           verbosity: int = 0) -> None:
        """
        **Description:**
        Record a set of points that is not dependent

        **Arguments:**
        - `label_inds`: The iterable of vector/label indices.
        - `verbosity`:  The verbosity level.

        **Returns:**
        Nothing.
        """
        encoding = self.encode(label_inds)

        new_non_dependencies= set()
        for non_dependency in self.non_dependencies:
            if self.is_subset(non_dependency, encoding):
                # non_dependency is weaker than encoding...
                # don't save it in our new list
                if verbosity >= 1:
                    print(f"Outdated non-dependency = {non_dependency}...")

                pass
            else:
                # not a subset of encoding -> not trivial
                new_non_dependencies.add(non_dependency)

        new_non_dependencies.add(encoding)
        self.non_dependencies = new_non_dependencies

    # dictionary methods
    # ------------------
    def values(self) -> Iterable["Circuit"]:
        """
        **Description:**
        Get the values (the actual circuits)

        **Arguments:**
        None

        **Returns:**
        The circuits.
        """
        return self.circuits.values()

    def copy(self) -> "Circuits":
        """
        **Description:**
        Copy the circuits object

        **Arguments:**
        None

        **Returns:**
        A copy of the circuits.
        """
        copied = Circuits()
        copied.circuits = {Z:circ for Z,circ in self.circuits.items()}
        copied.cone_to_circuit = {c:Zs.copy() for c,Zs in \
                                                self.cone_to_circuit.items()}
        copied.non_dependencies = self.non_dependencies.copy()
        copied.know_all_circuits = self.know_all_circuits

        return copied

    def pop(self, *args, **kwargs):
        """
        Pop an element from the circuits dict
        """
        out = self.circuits.pop(*args, **kwargs)

    # basic bit helpers
    # -----------------
    def encode(self, label_inds: Iterable[int]) -> int:
        """
        **Description:**
        Convert an iterable of integers to a binary vector, b, such that
            b_i = 1 <=> i in label_inds

        **Arguments:**
        - `label_inds`: The iterable of integers.

        **Returns:**
        The encoding
        """
        # as bitvector
        if isinstance(label_inds, int):
            return label_inds

        encoding = 0
        for label_ind in label_inds:
            encoding |= (1 << int(label_ind))
        return encoding

    def decode(self, encoding) -> list[int]:
        """
        **Description:**
        Convert a binary vector b to a list of of integers such that
            b_i = 1 <=> i in label_inds

        **Arguments:**
        - `encoding`: The encoding to map to label indices

        **Returns:**
        The label indices
        """
        # as bitvector
        label_inds = []

        for shift in range(len(bin(encoding))-2):
            if 1&(encoding>>shift):
                label_inds.append(shift)

        return label_inds
    
    def is_superset(self, setA, setB) -> bool:
        """
        **Description:**
        Check if the set encoded by setA is a superset of setB.

        **Arguments:**
        - `setA`: The candidate-superset encoding.
        - `setB`: The candidate-subset encoding.

        **Returns:**
        Whether setA is a superset of setB.
        """
        # as bitvector
        return (setA & setB) == setB
    
    def is_subset(self, setA: int, setB: int) -> bool:
        """
        **Description:**
        Check if the set encoded by setA is a subset of setB.

        **Arguments:**
        - `setA`: The candidate-superset encoding.
        - `setB`: The candidate-subset encoding.

        **Returns:**
        Whether setA is a subset of setB.
        """
        return self.is_superset(setB, setA)
