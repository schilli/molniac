#!/usr/bin/env python

from __future__ import print_function, division

import molniac as mol
import numpy as np

trj = mol.load("1GNUh.pdb")
trj.mutate_protein_residue(oldres="resSeq 13", newres="s488")
trj.save('1GNUh_K13S488.pdb')

trj = mol.load("1GNUh_K13S488.pdb")
trj.mutate_protein_residue(oldres="resSeq 62", newres="s594")
trj.save('1GNUh_K13S488_F62S594.pdb') 


