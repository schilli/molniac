#!/usr/bin/env python

from __future__ import print_function, division

import molniac as mol
import numpy as np

trj     = mol.load("1CRNh.pdb")
#origtrj = mol.load("1CRNh.pdb")

newrestrj = trj.mutate_protein_residue(oldres="resSeq 46", newres="trp")


trj.save('1CRNh_T46S488.pdb')
#newrestrj.save('ala.pdb')
