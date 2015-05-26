#!/usr/bin/env python
# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
##

"""
    test_conversion
    ~~~~~~~~~~~~~~

    Test NWChem-to-others conversion routines
"""

from collections import OrderedDict
import sys
import unittest
from src.EMSL_local import EMSL_local
from src import conversion

class ConversionTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def parse_nwchem(self, basis_name, element_symbol):
        #Test conversion from NWChem format to neutral intermediate
        c = conversion.Converter()
        el = EMSL_local("db/NWChem.db", fmt="nwchem", debug=False)
        basis = "\n".join(el.get_basis(basis_name, [element_symbol]))
        parsed = c.parse_one_nwchem(basis)
        return parsed

    def test_functions_per_shell_basic(self):
        #check functions per shell on parsed data
        reference = OrderedDict([("S", 5), ("P", 4), ("D", 2), ("F", 1)])
        parsed = self.parse_nwchem("cc-pVTZ", "Cl")
        self.assertEqual(reference, parsed.functions_per_shell)

    def test_functions_per_shell_sp(self):
        #check function ordering for basis set with SP functions
        reference = OrderedDict([("S", 1), ("SP", 6), ("D", 1)])
        parsed = self.parse_nwchem("6-31G*", "Cl")
        self.assertEqual(reference, parsed.functions_per_shell)


    def xtest_find_limits(self):
        c = conversion.Converter()
        counts = []
        shells = set()
        elements = [s[0] for s in c.elements[1:100]]
        el = EMSL_local("db/NWChem.db", fmt="nwchem", debug=False)
        for element in elements:
            for name, description in el.get_list_basis_available([element]):
                basis_raw = el.get_basis(name, [element])
                basis = "\n".join(basis_raw)
                parsed = c.parse_one_nwchem(basis)
                functions = parsed.functions_per_shell
                try:
                    value = max(functions.values())
                    for k in functions.keys():
                        shells.add(k)
                except ValueError:
                    continue
                key = (element, name)
                counts.append((value, key))
                print counts[-1], functions

        counts.sort()
        print counts[-1]
        print shells


 

def runSuite(cls, verbosity=2, name=None):
    """Run a unit test suite and return status code.

    @param cls: class that the suite should be constructed from
    @type cls : class
    @param verbosity: verbosity level to pass to test runner
    @type verbosity : int
    @param name: name of a specific test in the suite to run
    @type name : str
    @return: unit test run status code
    @rtype : int
    """
    try: 
        if name:
            suite = unittest.makeSuite(cls, name)
        else:
            suite = unittest.makeSuite(cls)
            
        return unittest.TextTestRunner(verbosity=verbosity).run(suite)
    
    except SystemExit:
        pass

def runTests():
    try:
        test_name = sys.argv[1]
        
    except IndexError:
        test_name = None

    if test_name:
        result = runSuite(ConversionTestCase, name = test_name)

    else:
        result = runSuite(ConversionTestCase)

    return result

if __name__ == '__main__':
    runTests()
