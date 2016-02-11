#!/usr/bin/env python2
# generate amino acid pdb files with pymol

from __future__ import print_function, division

import sys, pymol
from pymol.exporting import _resn_to_aa as one_letter

# reclaim stdout from pymol
stdout = sys.stdout
pymol.finish_launching()
sys.stdout = stdout

for resname in one_letter.keys():
    resname = resname.lower()
    pymol.cmd.reinitialize()
    pymol.finish_launching()
    pymol.editor.attach_amino_acid('pk1', resname)
    pymol.cmd.save(filename = resname + '.pdb', selection=resname)
    print("Generated {}.pdb".format(resname))

