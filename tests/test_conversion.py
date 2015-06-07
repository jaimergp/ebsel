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
        #conversion from NWChem format to neutral intermediate
        c = conversion.Converter()
        el = EMSL_local(fmt="nwchem", debug=False)
        basis = "\n".join(el.get_basis(basis_name, [element_symbol]))
        parsed = c.parse_one_nwchem(basis)
        return parsed

    def parse_g94(self, basis_name, element_symbol):
        #conversion from Gaussian 94 format to neutral intermediate
        c = conversion.Converter()
        el = EMSL_local(fmt="g94", debug=False)
        basis = "\n".join(el.get_basis(basis_name, [element_symbol]))
        parsed = c.parse_one_g94(basis)
        return parsed

    def test_functions_per_shell_basic_nwchem(self):
        #check functions per shell on parsed data
        reference = OrderedDict([("S", 5), ("P", 4), ("D", 2), ("F", 1)])
        parsed = self.parse_nwchem("cc-pVTZ", "Cl")
        self.assertEqual(reference, parsed.functions_per_shell)

    def test_functions_per_shell_sp_nwchem(self):
        #check function ordering for basis set with SP functions
        reference = OrderedDict([("S", 1), ("SP", 6), ("D", 1)])
        parsed = self.parse_nwchem("6-31G*", "Cl")
        self.assertEqual(reference, parsed.functions_per_shell)

    def test_functions_per_shell_basic_g94(self):
        #check functions per shell on parsed data
        reference = OrderedDict([("S", 5), ("P", 4), ("D", 2), ("F", 1)])
        parsed = self.parse_g94("cc-pVTZ", "Cl")
        self.assertEqual(reference, parsed.functions_per_shell)

    def test_functions_per_shell_sp_g94(self):
        #check function ordering for basis set with SP functions
        reference = OrderedDict([("S", 1), ("SP", 6), ("D", 1)])
        parsed = self.parse_g94("6-31G*", "Cl")
        self.assertEqual(reference, parsed.functions_per_shell)

    def test_reformat_functions(self):
        #translate from "wide" form to "tall" form
        wide = ('P',
                [[663.3, 0.00240448, -0.000652145],
                 [156.8, 0.0192148, -0.00519445],
                 [49.98, 0.0885097, -0.0246938],
                 [18.42, 0.25602, -0.0728167],
                 [7.24, 0.436927, -0.13403],
                 [2.922, 0.350334, -0.0947742],
                 [0.3818, -0.00458423, 0.564667]])
        tall =  ('P',
                 [[[663.3, 0.00240448],
                   [156.8, 0.0192148],
                   [49.98, 0.0885097],
                   [18.42, 0.25602],
                   [7.24, 0.436927],
                   [2.922, 0.350334],
                   [0.3818, -0.00458423]],
                  [[663.3, -0.000652145],
                   [156.8, -0.00519445],
                   [49.98, -0.0246938],
                   [18.42, -0.0728167],
                   [7.24, -0.13403],
                   [2.922, -0.0947742],
                   [0.3818, 0.564667]]])
        parsed = self.parse_nwchem("cc-pVTZ", "Cl")
        self.assertEqual(wide, parsed.functions[3])
        converted = parsed._reformat_functions(parsed.functions)
        self.assertEqual(tall, converted[3])

    def test_compare_reformat(self):
        #reformatted data should have the same layout
        args = ("cc-pVTZ", "Li")
        p1 = self.parse_nwchem(*args)
        p2 = self.parse_g94(*args)
        self.assertEqual(p1.functions_per_shell, p2.functions_per_shell)
        c1 = p1._reformat_functions(p1.functions)
        c2 = p2._reformat_functions(p2.functions)
        self.assertEqual(c1, c2)

    def test_convert_from_nwchem_to_nwchem(self):
        #test 'internal' conversion from nwchem to nwchem
        basis_name = "cc-pVTZ"
        el = EMSL_local(fmt="nwchem", debug=False)
        data = el.get_basis("cc-pVTZ", ["Cl"], convert_from="")
        data2 = el.get_basis("cc-pVTZ", ["Cl"], convert_from="nwchem")
        self.assertTrue("BASIS SET reformatted" in data2[0])

    def test_convert_from_nwchem_to_gamess(self):
        #test 'internal' conversion from nwchem to gamess-us
        basis_name = "cc-pVTZ"
        el = EMSL_local(fmt="gamess-us", debug=False)
        data = el.get_basis(basis_name, ["Cl"], convert_from="")
        data2 = el.get_basis(basis_name, ["Cl"], convert_from="nwchem")
        self.assertTrue("BASIS SET reformatted" in data2[0])
        self.assertTrue("CHLORINE" in data2[0])

    def test_convert_from_nwchem_to_g94(self):
        #test 'internal' conversion from nwchem to Gaussian 94
        basis_name = "cc-pVTZ"
        el = EMSL_local(fmt="g94", debug=False)
        data = el.get_basis(basis_name, ["Cl"], convert_from="")
        data2 = el.get_basis(basis_name, ["Cl"], convert_from="nwchem")
        self.assertTrue("BASIS SET reformatted" in data2[0])

    def xtest_find_limits(self):
        c = conversion.Converter()
        counts = []
        shells = set()
        elements = [s[0] for s in c.elements[1:100]]
        el = EMSL_local(fmt="nwchem", debug=False)
        for element in elements:
            for name, description in el.get_available_basis_sets([element]):
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
