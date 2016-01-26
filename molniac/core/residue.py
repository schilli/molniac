from __future__ import print_function, division

categories = ["", "protein", "dna", "rna", "water", "ion"]
secondaryStructures = ["", "helix", "beta", "loop"]

class Residue(object):
    """
    A residue is a logical unit of atoms.
    """

    def __init__(self, ID=-1, name="", atoms=None, molecule=None, category="", secstruct=""):

        self.ID        = ID         # residue identifier (resid in the PDB)
        self.name      = name       # name of the residue
        self.molecule  = molecule   # reference to the molecule this residue belongs to
        self.category  = category   # the biochemical category
        self.secstruct = secstruct  # secondary structure (only applicable to protein residues)
        self.molecule  = molecule   # reference to the molecule this residue belongs to

        if atoms is None:
            self.atoms = []
        else:
            self.atoms     = atoms      # a list of atoms belonging to this residue
        
