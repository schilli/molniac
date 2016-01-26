from __future__ import print_function, division

class Selection(object):
    """
    A selection contains a set of atoms organized as a hyrarchy of:
    molecules -> residues -> atoms
    Subselections can be created.
    """

    def __init__(self, molecules=None, parent=None):

        if molecules is None:
            self.molecules = []
        else:
            self.molecules = molecules  # list of molecules in the selection

        self.parent = parent            # if the selection was created as a subselection, it contains a reference to its parent selection


