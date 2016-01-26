from __future__ import print_function, division

import sys
from . import pdb_tests

def test(name=None):
    """
    Run the specified test or all tests.
    The test name is everything before the "_test.py" in the filename.
    """

    if name is "pdb":
        pdb_tests.test()

    elif name is None:
        pdb_tests.test()

    else:
        print("Test unknown: {}".format(test), file=sys.stderr)
        sys.exit(1)




