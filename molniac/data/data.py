from __future__ import print_function, division

import os


def _load_list(filename):
    """
    Load list of strings from file
    """
    List = []
    with open(filename, 'r') as f:
        lines = f.readlines()
    for line in lines:
        if line[0] != '#' and len(line.strip()) > 0:
            List.append(line.strip())
    return List


class Element(object):
    """
    A container for informatin on an element
    """

    def __init__(self, name="", number=0, mass=1.0, Rvdw=3.0):
        self.name    = name    # name of the element
        self.number  = number  # number in the periodic table
        self.mass    = mass    # mass in u
        self.Rvdw    = Rvdw    # VdW radius


 
def _load_elements(filename):
    """
    Load information on chemical elements from a file
    """
    elements = {}
    with open(filename, 'r') as f:
        lines = f.readlines()
    for line in lines:
        if line[0] != '#' and len(line.strip()) > 0: 
            fields  = line.split()
            element = Element(name=fields[0].title(), number=int(fields[1]),
                    mass=float(fields[2]), Rvdw=float(fields[3]))
            elements[element.name] = element
    return elements



class Restype(object):
    """
    Holds information on residue types
    """

    def __init__(self, name, category, aliases):
        self.name     = name     # residue name
        self.category = category # category (protein, dna, rna, etc.
        self.aliases  = aliases  # list of alias names for this residue
        self.upper()

    def upper(self):
        """
        Convert all strings to uppercase
        """
        self.name     = self.name.upper()    
        self.category = self.category.upper()
        for i in range(len(self.aliases)):
            self.aliases[i] = self.aliases[i].upper() 


class RestypeDB(object):
    """
    Residue database
    """

    def __init__(self):
        self.restypes   = []
        self.categories = {}
        self.aliases    = {}

    def add(self, restype):
        """
        Add a residue type to the database
        """
        restype.upper()
        self.restypes.append(restype)

        try:
            self.categories[restype.category].append(restype.name)
        except KeyError:
            self.categories[restype.category] = [restype.name]
        for alias in restype.aliases:
            self.categories[restype.category].append(alias)

        for alias in restype.aliases:
            self.aliases[alias] = restype


    def belongs_to_category(self, resname, category):
        """
        Check if residue belongs to category
        Faster than function category
        """
        if resname.upper() in self.categories[category.upper()]:
            return True
        else:
            return False


    def category(self, resname):
        """
        Find category that resname belongs to.
        """
        for category in self.categories:
            if resname.upper() in self.categories[category]:
                return category
        return "LIGAND"





def _load_residue_database(filename):
    """
    Load database of residues
    """
    restypeDB = RestypeDB()
    with open(filename, 'r') as f:
        lines = f.readlines()
    for line in lines:
        if line[0] != '#' and len(line.strip()) > 0: 
            fields   = line.split(',')
            resname  = fields[0].strip()
            category = fields[1].strip()
            aliases  = fields[2].split()
            restype  = Restype(resname, category, aliases)
            restypeDB.add(restype)
    return restypeDB




def is_protein(resname):
    """
    Check if the reside is a protein residue
    """
    return restypeDB.belongs_to_category(resname, 'PROTEIN') 


def is_water(resname):
    """
    Check if the reside is a water molecule
    """
    return restypeDB.belongs_to_category(resname, 'WATER')


def is_solvent(resname):
    """
    Check if the given resname denotes a solvent molecule (including water)
    """
    iswater   = restypeDB.belongs_to_category(resname, 'WATER')
    issolvent = restypeDB.belongs_to_category(resname, 'SOLVENT')
    return iswater or issolvent
 

def is_ion(resname):
    """
    Check if the given atom is an ion
    """
    return restypeDB.belongs_to_category(resname, 'ION')


def is_dna(resname):
    """
    Check if the given atom is DNA
    """
    return restypeDB.belongs_to_category(resname, 'DNA')


def is_rna(resname):
    """
    Check if the given atom is RNA
    """
    return restypeDB.belongs_to_category(resname, 'RNA')
 

def is_nucleic(resname):
    """
    Check if the given atom is a nucleic acid
    """
    isdna = restypeDB.belongs_to_category(resname, 'DNA')
    isrna = restypeDB.belongs_to_category(resname, 'RNA')
    return isdna or isrna


def is_polymer(resname):
    """
    Check if resname belongs to a polymeric species (i.e. protein or nucleic acid)
    """
    isprotein = is_protein(resname)
    isnucleic = is_nucleic(resname)
    return isprotein or isnucleic


def is_ligand(resname):
    """
    Check if the given atom is a ligand (i.e. does not belong to any other category)
    """
    return restypeDB.belongs_to_category(resname, 'LIGAND')
 

def get_category(resname):
    """
    Determine residue category of resname
    e.g.: protein, water, ion, ligand, etc.
    """
    return restypeDB.category(resname)



datapath          = os.path.dirname(os.path.realpath(__file__))
restypeDB         = _load_residue_database(datapath + "/residues.dat")
elements          = _load_elements(datapath + "/elements.dat")
 
