# Test PDB file reading functionality

from __future__ import print_function, division

import sys, os
import molniac



def test():
    """
    Run all PDB file tests
    """
    status = "passed"
    print("Running PDB tests ...", file=sys.stderr)

    if read_1crn() != "passed": status = "failed"

    print("\t*PDB tests " + status +'*', file=sys.stderr)
    return status




def read_1crn():
    """
    Tests basic PDB file reading
    """
    status = "passed"

    testpath = os.path.realpath(__file__) 
    filename = os.path.dirname(testpath) + "/1CRN.pdb"
    molniac.read(filename)


    print("\tReading 1CRN.pdb \t\t" + status, file=sys.stderr)

    return status
