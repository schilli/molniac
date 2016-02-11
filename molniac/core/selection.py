from __future__ import print_function, division

import sys
import numpy as np
import molniac.data as data

class Selection(object):
    """
    A selection contains a set of atoms organized as a hyrarchy of:
    molecules -> residues -> atoms
    Subselections can be created.
    """

    def __init__(self, parent=None):
        self.molecules = []     # list of molecules in the selection
        self.parent    = parent # if the selection was created as a subselection, it contains a reference to its parent selection
        self.atomlist  = []     # references to atoms, corresponding to xyz array

        # memory management
        self._xyz      = np.zeros((0,3))   # internal memory to store atom positions. The xyz array should always be a view of this private array
        self._Rvdw     = np.zeros(0)       # same as before for Rvdw
        self._atom_ndx = np.zeros(0, dtype=np.int)  # and for atom indices
        self._used     = np.zeros(0, dtype=np.bool) # keep track of memory elements that are in use

        # arrays for coordinates, etc.
        self.xyz       = self._xyz      # contains all atom positions in an array for fast access
        self.Rvdw      = self._Rvdw     # array of VdW radii for each atom for fast access
        self.atom_ndx  = self._atom_ndx # indices linking coordinates and Rvdw to self.atomlist list 


    def __str__(self):
        nmols = len(self.molecules)
        nres  = 0
        natms = 0
        for mol in self.molecules:
            nres += len(mol.residues)
            for res in mol.residues:
                natms += len(res.atoms)
        return "selection of {} molecules, {} residues, {} atoms".format(nmols, nres, natms)


    def show(self):
        """
        Print all molecules, residues and atoms to stdout
        """
        print(self)
        for mol in self.molecules:
            print('\t', mol)
            for res in mol.residues:
                print(2*'\t', res)
                for atm in res.atoms:
                    print(3*'\t', atm)



    def add(self, molecule):
        """
        Add a molecule
        """
        self.molecules.append(molecule)
        for atom in molecule:
            self._add_atom(atom)


    def add_atom(self, atom):
        """
        Add atom coordinates, etc. to the self managed xyz array
        """
        self.atomlist.append(atom)

        # check for sufficient free allocated memory
        if self._used.sum() >= self._used.shape[0]:
            self._allocate_memory()

        # find first free slot in self._xyz
        try:
            free_ndx = np.nonzero(np.logical_not(self._used))[0].min()
        except ValueError:
            free_ndx = 0

        # insert atom
        self._xyz     [free_ndx,:] = atom.xyz
        self._Rvdw    [free_ndx  ] = data.elements[atom.element].Rvdw
        self._atom_ndx[free_ndx  ] = len(self.atomlist) - 1
        self._used    [free_ndx  ] = True

        # set views
        self.xyz       = self._xyz     [self._used,:]
        self.Rvdw      = self._Rvdw    [self._used  ]
        self.atom_ndx  = self._atom_ndx[self._used  ]
 


    def residues(self):
        """
        Generator for iteration through residues
        """
        for mol in self.molecules:
            for res in mol.residues:
                yield res


    def atoms(self):
        """
        Generator for iteration through atoms
        """
        for res in self.residues():
            for atm in res.atoms:
                yield atm
    


    def detect_bonds(self):
        """
        Walk through the selection and autodetect bonds between atoms.
        """ 
        for atm in self.atoms():
            atm.bonds = [] 

        for atm_idx in range(1, self.xyz.shape[0]):
            coords = self.xyz[atm_idx,:]
            Rvdw1  = self.Rvdw[atm_idx]
            dist   = ((self.xyz[:atm_idx] - coords)**2).sum(1)**0.5
            Rvdw   = 0.5 * (self.Rvdw[:atm_idx] + Rvdw1)
            
            atm1 = self.atomlist[self.atom_ndx[atm_idx]]
#            for atm2_idx in self.atom_ndx[:atm_idx][Rvdw < 1.1]:
#                atm2 = self.atomlist[atm2_idx]
#                atm1.bonds.append(atm2)
#                atm2.bonds.append(atm1)






    def detect_bonds_old(self):
        """
        Walk through the selection and autodetect bonds between atoms.
        """
        for atm1 in self.atoms():
            atm1.bonds = []

        for i1, atm1 in enumerate(self.atoms()):
            for i2, atm2 in enumerate(self.atoms()):
                if i2 >= i1:
                    break
                dist  = np.linalg.norm(atm1.xyz - atm2.xyz)
                Rvdw1 = data.elements[atm1.element].Rvdw
                Rvdw2 = data.elements[atm2.element].Rvdw
                Rvdw  = 0.5 * (Rvdw1 + Rvdw2)

                # detect if there is a bond
                if dist <= 1.1*Rvdw:
                    atm1.bonds.append(atm2)
                    atm2.bonds.append(atm1)



    def _allocate_memory(self):
        """
        Allocate more memory
        """
        # determine new array length
        oldarraylength = self._xyz.shape[0]
        newarraylength = 2 * oldarraylength
        if newarraylength < 1:
                newarraylength = 1

        # allocate new memory
        _xyz      = np.zeros((newarraylength,3))
        _Rvdw     = np.zeros( newarraylength)    
        _atom_ndx = np.zeros( newarraylength, dtype=np.int) 
        _used     = np.zeros( newarraylength, dtype=np.bool)

        # copy data to new memory
        _xyz     [:oldarraylength,:] = self._xyz     
        _Rvdw    [:oldarraylength  ] = self._Rvdw    
        _atom_ndx[:oldarraylength  ] = self._atom_ndx
        _used    [:oldarraylength  ] = self._used    

        # replace old memory
        self._xyz      = _xyz     
        self._Rvdw     = _Rvdw    
        self._atom_ndx = _atom_ndx
        self._used     = _used    

        # set correct array views
        self.xyz      = self._xyz     [self._used,:]
        self.Rvdw     = self._Rvdw    [self._used  ]
        self.atom_ndx = self._atom_ndx[self._used  ]
