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
        el = EMSL_local(fmt="gamess-us", debug=False)
        el.get_basis("pCs-3", ["Cl"])
        self.assertFalse(el.am_too_large)
        self.assertEquals("G", el.max_am)

    def test_gamess_us_am_fail(self):
        #GAMESS-US angular momentum check fails for max am <= I
        el = EMSL_local(fmt="gamess-us", debug=False)
        el.get_basis("cc-pv6z", ["Cl"])
        self.assertFalse(el.am_too_large)
        self.assertEqual("I", el.max_am)

    def test_gamess_us_am_L(self):
        #GAMESS-US angular momentum check special case for SP "L" basis
        el = EMSL_local(fmt="gamess-us", debug=False)
        el.get_basis("6-31G", ["Cl"])
        self.assertFalse(el.am_too_large)
        self.assertEqual("P", el.max_am)

    def test_nwchem_am_pass(self):
        #NWChem angular momentum check passes for max am <= I
        el = EMSL_local(fmt="nwchem", debug=False)
        el.get_basis("cc-pv6z", ["Ne"])
        self.assertFalse(el.am_too_large)
        self.assertEqual("I", el.max_am)

    def test_nwchem_am_fail(self):
        #NWchem angular momentum check fails for max am > I
        el = EMSL_local(fmt="nwchem", debug=False)
        el.get_basis("cc-pv8z", ["Ne"])
        self.assertTrue(el.am_too_large)
        self.assertEqual("L", el.max_am)

    def test_gaussian94_am(self):
        #There is no upper am limit for this format! But verify max_am
        el = EMSL_local(fmt="g94", debug=False)
        el.get_basis("cc-pv8z", ["Ne"])
        self.assertFalse(el.am_too_large)
        self.assertEqual("L", el.max_am)

    def test_cartesian_or_spherical(self):
        #most basis sets treated as using spherical (pure) functions, while
        #a few older ones are treated as using cartesians
        expected_cartesian = ["3-21G", "4-31G", "6-31G", "6-31G*",
                              "6-31G**", "DZ (Dunning)", "DZP (Dunning)"] 
        el = EMSL_local(fmt="nwchem", debug=False)
        assigned = {}
        names = el.get_available_basis_sets()
        names.sort()
        for name, description in names:
            fn_type = el.spherical_or_cartesian(name)
            try:
                assigned[fn_type].append(name)
            except KeyError:
                assigned[fn_type] = [name]

        self.assertEqual(expected_cartesian, assigned["cartesian"])

    def test_get_available_basis_sets_name_filter(self):
        #test get_available_basis_sets with basis set name filtering
        el = EMSL_local(fmt="nwchem")
        basis_names = ["cc-pVTZ", "fakename", "6-31G"]
        expected = [("6-31G", "VDZ Valence Double Zeta: 2 Funct.'s/Valence AO"),
                    ("cc-pVTZ", "VTZ2P Valence Triple Zeta + Polarization on All Atoms")]

        names = el.get_available_basis_sets(allowed_basis_names=basis_names)
        self.assertEqual(expected, names)

    def test_get_available_basis_sets_element_filter(self):
        #test get_available_basis_sets with element filtering
        el = EMSL_local(fmt="nwchem")
        elements = ["Hg", "Pb", "Th"]
        expected = [("ANO-RCC",
                     "full ANO-RCC basis, reduce to get MB, VDZP, VTZP and VQZP quality"),
                    ("CRENBL ECP", "N/A"),
                    ("CRENBL ECP-number2", "1D UNCONTR Uncontracted"),
                    ("SARC-DKH", "N/A"),
                    ("SARC-ZORA",
                     "Segmented all-electron relativistically contracted basis sets for ZORA"),
                    ("Stuttgart RLC ECP", "DZ Double Zeta Basis Set designed for an ECP"),
                    ("Stuttgart RLC ECP-number2", "N/A"),
                    ("UGBS", "UGBS basis by de Castro and Jorge")]

        names = el.get_available_basis_sets(elements=elements)
        self.assertEqual(expected, names)

    def test_get_available_basis_sets_combined_filter(self):
        #test get_available_basis_sets with element + basis set name filtering
        el = EMSL_local(fmt="nwchem")
        elements = ["Hg", "Pb", "Th"]
        basis_names = ["SARC-ZORA", "ANO-RCC"]
        expected = [("ANO-RCC",
                     "full ANO-RCC basis, reduce to get MB, VDZP, VTZP and VQZP quality"),
                    ("SARC-ZORA",
                     "Segmented all-electron relativistically contracted basis sets for ZORA")]

        names = el.get_available_basis_sets(elements=elements,
                                            allowed_basis_names=basis_names)
        self.assertEqual(expected, names)

    def test_get_available_basis_sets_fs_basic(self):
        #test that we can get the name of supplemental basis set stored on
        #the file system
        el = EMSL_local(fmt="nwchem")
        expected = [("g3mp2large", "db/nwchem/g3mp2large.nwbas"),
                    ("g3largexp", "db/nwchem/g3largexp.nwbas")]
        names = el.get_available_basis_sets_fs("nwchem")
        self.assertEqual(sorted(expected), sorted(names))

    def test_get_available_elements_fs(self):
        #test that we can get elements from supplemental basis set stored on
        #the file system
        el = EMSL_local(fmt="nwchem")
        expected = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
                    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca',
                    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']

        names = el.get_available_elements_fs("nwchem", "g3mp2large")
        self.assertEqual(expected, names)

    def test_get_available_elements(self):
        #verify element listing from standard db-stored basis set
        el = EMSL_local(fmt="nwchem")
        expected = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
                    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca',
                    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', "I"]

        names = el.get_available_elements("6-311G")
        self.assertEqual(expected, names)

    def test_get_available_elements_fused(self):
        #element data is automatically supplemented via basis sets stored
        #on the file system
        el = EMSL_local(fmt="nwchem")
        expected = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
                    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca',
                    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']

        names = el.get_available_elements("g3mp2large")
        self.assertEqual(expected, names)


    def test_get_available_basis_sets_fs_name_filter(self):
        #test that we can get the name of supplemental basis set stored on
        #the file system -- with name filtering
        el = EMSL_local(fmt="nwchem")
        expected = []
        basis_names = ["g3mp2gigantic"]
        names = el.get_available_basis_sets_fs("nwchem", allowed_basis_names=basis_names)
        self.assertEqual(expected, names)

    def test_get_available_basis_sets_fs_element_filter(self):
        #test that we can get the name of supplemental basis set stored on
        #the file system -- with element filtering
        el = EMSL_local(fmt="nwchem")
        expected1 = [("g3mp2large", "db/nwchem/g3mp2large.nwbas")]
        expected2 = []
        elements1 = ["Ar", "Kr"]
        elements2 = ["Kr", "Xe"]

        #g3mp2large has krypton parameters but not xenon, so second retrieval
        #should produce nothing
        names1 = el.get_available_basis_sets_fs("nwchem", elements=elements1)
        self.assertEqual(expected1, names1)

        names2 = el.get_available_basis_sets_fs("nwchem", elements=elements2)
        self.assertEqual(expected2, names2)

    def test_get_available_basis_sets_supplemented(self):
        #test get_available_basis_sets supplemented with basis data from
        #the file system
        el = EMSL_local(fmt="nwchem")
        basis_names = ["cc-pVTZ", "g3mp2large", "6-31G"]
        expected = [("6-31G", "VDZ Valence Double Zeta: 2 Funct.'s/Valence AO"),
                    ("cc-pVTZ", "VTZ2P Valence Triple Zeta + Polarization on All Atoms"),
                    ("g3mp2large", "db/nwchem/g3mp2large.nwbas")]

        names = el.get_available_basis_sets(allowed_basis_names=basis_names)
        self.assertEqual(expected, names)

    def test_get_basis_supplemented(self):
        #test that basis set data gets auto-translated when there
        #is no "native" version available but an .nwbas to convert
        el = EMSL_local(fmt="g94")
        elements = ["Li", "Cl"]
        result = el.get_basis("g3mp2large", elements=elements)
        self.assertTrue("BASIS SET reformatted" in result[0])


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
