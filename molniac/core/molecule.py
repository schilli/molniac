from __future__ import print_function, division


class Molecule(object):
    """
    A molecule consisting of residues connected by covalent bonds.
    Atoms in this molecule can still be bound to other molecules by, e.g. sulfide bridges
    """

    def __init__(self, ID="", residues=None, selection=None):
        self.ID        = ID         # corresponding to the chain identifier in the PDB
        self.selection = selection  # a reference to the atom selection this molecule belongs to

        if residues is None:
            self.residues = []
        else:
            self.residues = residues   # a list of residues in this molecule

