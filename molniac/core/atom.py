from __future__ import print_function, division

import numpy as np

categories = ["", "backbone", "sidechain"]

class Atom(object):
    """
    Contains all information on a classical atom
    """

    def __init__(self, ID=0, name="", insertion="", occupancy=1.0, bfac=0.0, element="", charge=0.0, xyz=None, residue=None, category=""):
        """
        Initialize atomic data
        """

        self.ID        = ID          # an integer identifying the atom uniquely
        self.name      = name        # name of the atom
        self.insertion = insertion   # insertion code (as in pdb files)
        self.occupancy = occupancy   # occupancy of the atom (as in X-ray crystallography files)
        self.bfac      = bfac        # B-factor or temperature factor
        self.element   = element     # chemical element corresponding to the abbreviations in the periodic table
        self.charge    = charge      # partial charge on the atom
        self.category  = category    # sidechain or backbone, etc.
        self.residue   = residue     # reference to the residue the atom belongs to
        self.neighbors = []          # list of atoms to which a bond exists

        if xyz is None:
            self.xyz = np.zeros(3)
        else:
            self.xyz = xyz           # atomic position in angstrom

    def __str__(self):
        s  = "{:6d}".format(self.ID)
        s += " " + self.name
        return s
