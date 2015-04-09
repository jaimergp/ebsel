#!/usr/bin/env python
# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
##

"""
    test_local
    ~~~~~~~~~~~~~~

    Test data export functions
"""

import sys
import unittest
from src.EMSL_local import EMSL_local

class LocalTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_gamess_us_am_pass(self):
        #GAMESS-US angular momentum check passes for max am <= G
        el = EMSL_local("db/Gamess-us.db", fmt="gamess-us", debug=False)
        el.get_basis("pCs-3", ["Cl"])
        self.assertFalse(el.am_too_large)
        self.assertEquals("G", el.max_am)

    def test_gamess_us_am_fail(self):
        #GAMESS-US angular momentum check fails for max am <= I
        el = EMSL_local("db/Gamess-us.db", fmt="gamess-us", debug=False)
        el.get_basis("cc-pv6z", ["Cl"])
        self.assertFalse(el.am_too_large)
        self.assertEqual("I", el.max_am)

    def test_gamess_us_am_L(self):
        #GAMESS-US angular momentum check special case for SP "L" basis
        el = EMSL_local("db/Gamess-us.db", fmt="gamess-us", debug=False)
        el.get_basis("6-31G", ["Cl"])
        self.assertFalse(el.am_too_large)
        self.assertEqual("P", el.max_am)

    def test_nwchem_am_pass(self):
        #NWChem angular momentum check passes for max am <= I
        el = EMSL_local("db/NWChem.db", fmt="nwchem", debug=False)
        el.get_basis("cc-pv6z", ["Ne"])
        self.assertFalse(el.am_too_large)
        self.assertEqual("I", el.max_am)

    def test_nwchem_am_fail(self):
        #NWchem angular momentum check fails for max am > I
        el = EMSL_local("db/NWChem.db", fmt="nwchem", debug=False)
        el.get_basis("cc-pv8z", ["Ne"])
        self.assertTrue(el.am_too_large)
        self.assertEqual("L", el.max_am)

    def test_gaussian94_am(self):
        #There is no upper am limit for this format! But verify max_am
        el = EMSL_local("db/Gaussian94.db", fmt="g94", debug=False)
        el.get_basis("cc-pv8z", ["Ne"])
        self.assertFalse(el.am_too_large)
        self.assertEqual("L", el.max_am)

    def test_cartesian_or_spherical(self):
        #most basis sets treated as using spherical (pure) functions, while
        #a few older ones are treated as using cartesians
        expected_cartesian = [u'3-21G', u'4-31G', u'6-31G', u'6-31G*',
                              u'6-31G**'] 
        el = EMSL_local("db/NWChem.db", fmt="nwchem", debug=False)
        assigned = {}
        names = el.get_list_basis_available()
        names.sort()
        for name, description in names:
            fn_type = el.spherical_or_cartesian(name)
            try:
                assigned[fn_type].append(name)
            except KeyError:
                assigned[fn_type] = [name]

        self.assertEqual(expected_cartesian, assigned["cartesian"])
 

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
        result = runSuite(LocalTestCase, name = test_name)

    else:
        result = runSuite(LocalTestCase)

    return result

if __name__ == '__main__':
    runTests()
