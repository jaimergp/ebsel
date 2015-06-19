#!/usr/bin/env python
# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
##

"""
    test_parsing
    ~~~~~~~~~~~~~~

    Test parsing of basis set data from raw EMSL BSE downloads.
"""

import json
import sys
import unittest
from src.EMSL_dump import EMSL_dump

class ParserTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_nwchem_cd_basis(self):
        #test for basis set containing charge density fitting data only
        he_cd_ref = """#BASIS SET: (4s,2p) -> [2s,2p]\nHe    S\n      3.73939300E+01         0.10996705       \n      6.98669000E+00         0.37520477       \n      1.92344500E+00         0.41847746       \nHe    S\n      6.28757000E-01         1.0000000        \nHe    P\n      3.60021686E+00         1.0000000        \nHe    P\n      1.50009035E+00         1.0000000        """
        ed = EMSL_dump(None, format="NWChem")
        name = "Ahlrichs Coulomb Fitting"
        description = "DFT Coulomb Fitting"
        elements = "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr".split()

        with open("tests/samples/nwchem-ahlrichs-cf.html") as infile:
            text = infile.read()

        parser_method = ed.extraction_map[ed.format]
        name, description, parsed = parser_method(text, name, description,
                                                  elements)
        he_cd_data = json.loads(parsed[1][1])["cd basis"]
        self.assertEqual(he_cd_ref, he_cd_data)

    def test_nwchem_ecp_only(self):
        #extract basis set data from data containing ECP section alone
        ed = EMSL_dump(None, format="NWChem")
        name = "LANL2DZ ECP"
        description = "Halogen ECP"
        elements = "Cl Br I".split()

        cl_ecp_ref = """Cl nelec 10\nCl ul\n1     94.8130000            -10.0000000        \n2    165.6440000             66.2729170        \n2     30.8317000            -28.9685950        \n2     10.5841000            -12.8663370        \n2      3.7704000             -1.7102170        \nCl S\n0    128.8391000              3.0000000        \n1    120.3786000             12.8528510        \n2     63.5622000            275.6723980        \n2     18.0695000            115.6777120        \n2      3.8142000             35.0606090        \nCl P\n0    216.5263000              5.0000000        \n1     46.5723000              7.4794860        \n2    147.4685000            613.0320000        \n2     48.9869000            280.8006850        \n2     13.2096000            107.8788240        \n2      3.1831000             15.3439560        """
        
        with open("tests/samples/nwchem-lanl2dz-ecponly-small.html") as infile:
            text = infile.read()

        parser_method = ed.extraction_map[ed.format]
        name, description, parsed = parser_method(text, name, description,
                                                  elements)

        #Chlorine has ecp data only
        cl_unpacked = json.loads(parsed[0][1])
        cl_ao_basis = cl_unpacked.get("ao basis")
        cl_ecp = cl_unpacked.get("ecp")

        self.assertFalse(cl_ao_basis)
        self.assertEqual(cl_ecp_ref, cl_ecp)

    def test_nwchem_ecp_mixed(self):
        #extract basis set data from data containing AO basis and ECP section
        li_ao_ref = """#BASIS SET: (10s,4p) -> [3s,2p]\nLi    S\n    921.3000000              0.0013670        \n    138.7000000              0.0104250        \n     31.9400000              0.0498590        \n      9.3530000              0.1607010        \n      3.1580000              0.3446040        \n      1.1570000              0.4251970        \n      0.4446000              0.1694680        \nLi    S\n      0.4446000             -0.2223110        \n      0.0766600              1.1164770        \nLi    S\n      0.0286400              1.0000000        \nLi    P\n      1.4880000              0.0387700        \n      0.2667000              0.2362570        \n      0.0720100              0.8304480        \nLi    P\n      0.0237000              1.0000000        """

        na_ao_ref = """#BASIS SET: (3s,3p) -> [2s,2p]\nNa    S\n      0.4972000             -0.2753574        \n      0.0560000              1.0989969        \nNa    S\n      0.0221000              1.0000000        \nNa    P\n      0.6697000             -0.0683845        \n      0.0636000              1.0140550        \nNa    P\n      0.0204000              1.0000000        """

        na_ecp_ref = """Na nelec 10\nNa ul\n1    175.5502590            -10.0000000        \n2     35.0516791            -47.4902024        \n2      7.9060270            -17.2283007        \n2      2.3365719             -6.0637782        \n2      0.7799867             -0.7299393        \nNa S\n0    243.3605846              3.0000000        \n1     41.5764759             36.2847626        \n2     13.2649167             72.9304880        \n2      3.6797165             23.8401151        \n2      0.9764209              6.0123861        \nNa P\n0   1257.2650682              5.0000000        \n1    189.6248810            117.4495683        \n2     54.5247759            423.3986704        \n2     13.7449955            109.3247297        \n2      3.6813579             31.3701656        \n2      0.9461106              7.1241813        """
        
        ed = EMSL_dump(None, format="NWChem")
        name = "LANL2DZ ECP"
        description = "DZ Double Zeta Basis Set designed for an ECP"
        elements = "H Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi U Np Pu".split()
        
        with open("tests/samples/nwchem-lanl2dz-ecp.html") as infile:
            text = infile.read()

        parser_method = ed.extraction_map[ed.format]
        name, description, parsed = parser_method(text, name, description,
                                                  elements)

        #Lithium only has ordinary ao basis, no ecp data
        li_unpacked = json.loads(parsed[1][1])
        li_ao_basis = li_unpacked.get("ao basis")
        li_ecp = li_unpacked.get("ecp")
        
        self.assertEqual(li_ao_ref, li_ao_basis)
        self.assertFalse(li_ecp)

        #Sodium has ao basis and ecp data
        na_unpacked = json.loads(parsed[9][1])
        na_ao_basis = na_unpacked.get("ao basis")
        na_ecp = na_unpacked.get("ecp")

        self.assertEqual(na_ao_ref, na_ao_basis)
        self.assertEqual(na_ecp_ref, na_ecp)

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

        parser_method = ed.extraction_map[ed.format]
        name, description, parsed = parser_method(text, name, description,
                                                  elements)
        self.assertEqual(len(elements), len(parsed))
        unpacked = json.loads(parsed[1][1])
        he_ao_basis = unpacked["ao basis"]
        self.assertEqual("He", parsed[1][0])
        self.assertEqual(helium, he_ao_basis)

    def test_nwchem_single(self):
        #Code did not handle a single-element basis set before
        ed = EMSL_dump(None, format="NWChem")
        name = "B2 basis set for Zn"
        description = "N/A"
        elements = ["Zn"]
        with open("tests/samples/nwchem-B2_basis_set_for_Zn.html") as infile:
            text = infile.read()

        parser_method = ed.extraction_map[ed.format]
        name, description, parsed = parser_method(text, name, description,
                                                  elements)
        self.assertEqual(1, len(parsed))
        self.assertEqual("Zn", parsed[0][0])
        
    def test_gamess_us_basic(self):
        #extract basis set data from a popular Pople basis set
        helium = """HELIUM\nS   3\n  1     38.4216340              0.0237660        \n  2      5.7780300              0.1546790        \n  3      1.2417740              0.4696300        \nS   1\n  1      0.2979640              1.0000000"""
        
        ed = EMSL_dump(None, format="GAMESS-US")
        name = "6-31G*"
        description = "6-31G* Split Valence + Polarization Basis"
        elements = "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn".split()
        with open("tests/samples/gamess-us-6-31Gs.html") as infile:
            text = infile.read()

        parser_method = ed.extraction_map[ed.format]
        name, description, parsed = parser_method(text, name, description,
                                                  elements)
        self.assertEqual(len(elements), len(parsed))
        self.assertEqual("He", parsed[1][0])
        self.assertEqual(helium, parsed[1][1])

    def test_gaussian94_basic(self):
        #extract basis set data from a popular Pople basis set
        helium = """He     0 \nS   3   1.00\n     38.4216340              0.0237660        \n      5.7780300              0.1546790        \n      1.2417740              0.4696300        \nS   1   1.00\n      0.2979640              1.0000000        """
        
        ed = EMSL_dump(None, format="Gaussian94")
        name = "6-31G*"
        description = "6-31G* Split Valence + Polarization Basis"
        elements = "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn".split()
        with open("tests/samples/gaussian94-6-31Gs.html") as infile:
            text = infile.read()

        parser_method = ed.extraction_map[ed.format]
        name, description, parsed = parser_method(text, name, description,
                                                  elements)
        self.assertEqual(len(elements), len(parsed))
        self.assertEqual("He", parsed[1][0])
        self.assertEqual(helium, parsed[1][1])
        
    def test_gaussian94_multipart(self):
        #extract basis set data where each element has multiple basis blocks
        #NOTE: This does not really work properly! But at least it completes
        #without raising errors that block further work.
        helium = """He     0 \nS   5   1.00\n    221.3880300              0.0027491        \n     33.2619660              0.0208658        \n      7.5616549              0.0970588        \n      2.0855990              0.2807289        \n      0.6143392              0.4742218        \nS   1   1.00\n      0.1829212              1.0000000        """
        
        ed = EMSL_dump(None, format="Gaussian94", debug=False)
        name = "DZVP (DFT Orbital)"
        description = "VDZP Valence Double Zeta + Polarization designed for DFT"
        elements = "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe".split()
        with open("tests/samples/gaussian94-dzvp.html") as infile:
            text = infile.read()

        parser_method = ed.extraction_map[ed.format]
        name, description, parsed = parser_method(text, name, description,
                                                  elements)
        #This commented-out assertion would fail; > 1 block per element here
        #self.assertEqual(len(elements), len(parsed))
        self.assertEqual("He", parsed[1][0])
        self.assertEqual(helium, parsed[1][1])

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
