from __future__ import print_function, division

secondaryStructures = ["", "helix", "beta", "loop"]

class Residue(object):
    """
    A residue is a logical unit of atoms.
    """

    def __init__(self, ID=-1, name="", molecule=None, category="", secstruct=""):

        self.ID        = ID         # residue identifier (resid in the PDB)
        self.name      = name       # name of the residue
        self.molecule  = molecule   # reference to the molecule this residue belongs to
        self.category  = category   # the biochemical category
        self.secstruct = secstruct  # secondary structure (only applicable to protein residues)
        self.molecule  = molecule   # reference to the molecule this residue belongs to
        self.atoms     = []         # a list of atoms belonging to this residue


    def __str__(self):
        return "{:3s} {:6d} {:10s} {:10s}".format(self.name, self.ID, self.category, self.secstruct)


    def add(self, atom):
        """
        Add an atom to this residue
        """
        self.atoms.append(atom)

        # if this residue is part of a selection, add the atom there too
        if self.molecule is not None:
            if self.molecule.selection is not None:
                self.molecule.selection.add_atom(atom)


        
