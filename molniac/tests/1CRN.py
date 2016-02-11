#!/usr/bin/env python

from __future__ import print_function, division

import molniac as mol
import numpy as np

trj     = mol.load("1CRNh.pdb")
origtrj = mol.load("1CRNh.pdb")

newrestrj = trj.mutate_protein_residue(oldres="resSeq 24", newres="ala")

#for a1, a2 in zip(origtrj.top.atoms, trj.top.atoms):
#    if a1.index != a2.index or a1.serial != a2.serial or a1.name != a2.name:
#        print(a1, a2)
#    else:
#        p1   = origtrj.xyz[0,a1.index,:]
#        p2   =     trj.xyz[0,a2.index,:]
#        dist = np.linalg.norm(p1 - p2)
#        if dist > 1e-5:
#            print("dist = {:.3e}, ({:10s}, {:10s}), [{:.2e}, {:.2e}, {:.2e}], [{:.2e}, {:.2e}, {:.2e}]".format(dist, a1, a2, p1[0], p1[1], p1[2], p2[0], p2[1], p2[2]))



trj.save('1CRNh_A24A.pdb')
newrestrj.save('ala.pdb')
