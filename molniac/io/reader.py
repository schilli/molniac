from __future__ import print_function, division

import sys
from . import pdb_reader as pdb

def read(filename, **kwargs):
    """
    Load structure data from filename.
    The filetype is deduced from the filename extension.
    If applicable for the filetype, only the given model is loaded, or if None, the first model
    Supported file formats:
        PDB
    """

    extension = filename.split('.')[-1]
    if extension in ["pdb", "PDB"]:
        selection = pdb.read(filename, **kwargs)

#        for mol in selection.molecules:
#            print("mol", mol.ID)
#            for res in mol.residues:
#                print("{:3s} {:3d}".format(res.name, res.ID))
#                for atm in res.atoms:
#                    print("\t{:3s} {:5d} {:4.2f} {:2s} {:8.3f} {:8.3f} {:8.3f}".format(atm.name, atm.ID, atm.occupancy, atm.element, *atm.xyz))
    else:
        print("Filetype cannot be deduced from extension: {}".format(extension))
        sys.exit(1)

    return selection
 
