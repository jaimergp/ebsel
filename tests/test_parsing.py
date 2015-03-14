#!/usr/bin/env python
# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
##

"""
    test_parsing
    ~~~~~~~~~~~~~~

    Test parsing of basis set data from raw EMSL BSE downloads.
"""

import sys
import unittest
from src.EMSL_dump import EMSL_dump

class ParserTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_nwchem_basic(self):
        #extract basis set data from a popular Pople basis set
        helium = """#BASIS SET: (4s) -> [2s]
He    S
     38.4216340              0.0237660        
      5.7780300              0.1546790        
      1.2417740              0.4696300        
He    S
      0.2979640              1.0000000        """
        
        ed = EMSL_dump(None, format="NWChem")
        name = "6-31G*"
        description = "6-31G* Split Valence + Polarization Basis"
        elements = "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn".split()
        with open("tests/samples/nwchem-6-31Gs.html") as infile:
            text = infile.read()

        parser_method = ed.parser_map[ed.format]
        name, description, parsed = parser_method(text, name, description,
                                                  elements)
        self.assertEquals(len(elements), len(parsed))
        self.assertEquals("He", parsed[1][0])
        self.assertEquals(helium, parsed[1][1])

    def test_nwchem_single(self):
        #Code did not handle a single-element basis set before
        ed = EMSL_dump(None, format="NWChem")
        name = "B2 basis set for Zn"
        description = "N/A"
        elements = ["Zn"]
        with open("tests/samples/nwchem-B2_basis_set_for_Zn.html") as infile:
            text = infile.read()

        parser_method = ed.parser_map[ed.format]
        name, description, parsed = parser_method(text, name, description,
                                                  elements)
        self.assertEquals(1, len(parsed))
        self.assertEquals("Zn", parsed[0][0])
        
    def test_gamess_us_basic(self):
        #extract basis set data from a popular Pople basis set
        helium = """HELIUM\nS   3\n  1     38.4216340              0.0237660        \n  2      5.7780300              0.1546790        \n  3      1.2417740              0.4696300        \nS   1\n  1      0.2979640              1.0000000"""
        
        ed = EMSL_dump(None, format="GAMESS-US")
        name = "6-31G*"
        description = "6-31G* Split Valence + Polarization Basis"
        elements = "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn".split()
        with open("tests/samples/gamess-us-6-31Gs.html") as infile:
            text = infile.read()

        parser_method = ed.parser_map[ed.format]
        name, description, parsed = parser_method(text, name, description,
                                                  elements)
        self.assertEquals(len(elements), len(parsed))
        self.assertEquals("He", parsed[1][0])
        self.assertEquals(helium, parsed[1][1])

    def test_gaussian94_basic(self):
        #extract basis set data from a popular Pople basis set
        helium = """He     0 \nS   3   1.00\n     38.4216340              0.0237660        \n      5.7780300              0.1546790        \n      1.2417740              0.4696300        \nS   1   1.00\n      0.2979640              1.0000000        """
        
        ed = EMSL_dump(None, format="Gaussian94")
        name = "6-31G*"
        description = "6-31G* Split Valence + Polarization Basis"
        elements = "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn".split()
        with open("tests/samples/gaussian94-6-31Gs.html") as infile:
            text = infile.read()

        parser_method = ed.parser_map[ed.format]
        name, description, parsed = parser_method(text, name, description,
                                                  elements)
        self.assertEquals(len(elements), len(parsed))
        self.assertEquals("He", parsed[1][0])
        self.assertEquals(helium, parsed[1][1])

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
        result = runSuite(ParserTestCase, name = test_name)

    else:
        result = runSuite(ParserTestCase)

    return result

if __name__ == '__main__':
    runTests()
