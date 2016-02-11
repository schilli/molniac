from __future__ import print_function, division


class Molecule(object):
    """
    A molecule consisting of residues connected by covalent bonds.
    Atoms in this molecule can still be bound to other molecules by, e.g. sulfide bridges
    """

    def __init__(self, ID="", selection=None):
        self.ID        = ID         # corresponding to the chain identifier in the PDB
        self.selection = selection  # a reference to the atom selection this molecule belongs to
        self.residues  = []         # a list of residues in this molecule


    def __str__(self):
        return "Mol {:5s}".format(self.ID)


    def add(self, residue):
        """
        Add a residue to this molecule
        """
        self.residues.append(residue)


    def residues(self):
        """
        Generator for iteration through residues
        """
        for res in mol.residues:
            yield res


    def atoms(self):
        """
        Generator for iteration through atoms
        """
        for res in self.residues():
            for atm in res.atoms:
                yield atm 
 
