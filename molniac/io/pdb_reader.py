from __future__ import print_function, division

from ..core import *
import numpy as np

def read(filename, model=None):
    """
    Load structural data from a PDB file.
    If model is None, the first model ist read.
    """

    sel = selection.Selection()
    mol = molecule.Molecule()
    res = residue.Residue()

    if model is None:
        model = 1

    with open(filename, 'r') as f:
        lines = f.readlines() 

    currentModel = 0
    for line in lines:
        if line[0:5] == "MODEL":
            currentModel = int(line[5:])

        if line[0:3] == "TER":
            # reset molelcule and residue
            mol = molecule.Molecule()
            res = residue.Residue() 

        elif line[0:6].strip() in ["ATOM", "HETATM"]:

            if currentModel == 0:
                # there is only one model in the PDB file
                model = 0

            if model == currentModel:
                alt_loc = line[16:17]
                # onlye read the first of multiple alternate locations
                if alt_loc in [" ", "A"]: # actually other ways of specifying an alternate location are possible
                    atomid     =   int( line[ 6:11])         
                    atomname   =        line[12:16].strip()  
                    resname    =        line[17:20].strip()  
                    chain      =        line[21:22].strip()  
                    resid      =   int( line[22:26])         
                    insertion  =        line[26:27].strip() 
                    xyz        = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                    occupancy  = float( line[54:60])         
                    bfac       = float( line[60:66])        
                    element    =        line[76:78].strip()

                    try:
                        charge = float( line[78:80])          
                    except:
                        charge = 0.0

                    xyz = np.array(xyz)

                    if chain != mol.ID:
                        # start new molecule
#                        if mol.ID != "":
#                            sel.molecules.append(mol)
                        mol = molecule.Molecule(ID=chain, selection=sel)
                        sel.molecules.append(mol)

                    if resid != res.ID:
                        # start new residue
#                        if res.ID >= 0:
#                            mol.residues.append(res)
                        res = residue.Residue(ID=resid, name=resname, molecule=mol)
                        mol.residues.append(res)

                    atm = atom.Atom(ID=atomid, name=atomname, insertion=insertion, occupancy=occupancy, bfac=bfac, element=element, charge=charge, xyz=xyz, residue=res, category="")
                    res.atoms.append(atm)
#                    print(atm)

    return sel

