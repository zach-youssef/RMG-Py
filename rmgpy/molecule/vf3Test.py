#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2018 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

import unittest

from rmgpy.molecule import getElement
from rmgpy.molecule.vf2Test import TestVF2
from rmgpy.molecule.molecule import Molecule, Atom, Bond
from rmgpy.molecule.vf3 import VF3


###############################################################################

class TestVF3(TestVF2):
    """
    Applies the tests for the old VF2 algorithm to the
    updated VF3 implementation, along with some additional unit tests
    """

    def setUp(self):
        super(TestVF3, self).setUp()
        self.vf2 = VF3()
        self.vf3 = VF3()

    def testGraphIsomorphism(self):
        self.assertTrue(self.vf3.isIsomorphic(self.mol, self.mol2, None))

        examples = ["CC(C)=O", "C1=CC=CC=C1O[H]", "[S](=O)(=O)([O-])[O-]", "C1=CC=C2C(=C1)C=CC=C2"]

        for smiles in examples:
            mol = Molecule().fromSMILES(smiles)
            self.assertTrue(self.vf3.isIsomorphic(mol, mol.copy(True), None))
            self.assertTrue(self.vf3.isSubgraphIsomorphic(mol, mol.copy(True), None))
            for smiles2 in examples:
                if smiles is not smiles2:
                    self.assertFalse(self.vf3.isIsomorphic(mol, Molecule().fromSMILES(smiles2), None))
                    self.assertFalse(self.vf3.isSubgraphIsomorphic(mol, Molecule().fromSMILES(smiles2), None))

    def testSubgraphIsomorphism(self):
        mol = Molecule()

        c1 = Atom(getElement(6))
        c2 = c1.copy()

        mol.addAtom(c1)
        mol.addAtom(c2)
        mol.addBond(Bond(c1, c2))

        mol2 = Molecule()

        c1 = c1.copy()
        c2 = c1.copy()
        c3 = c1.copy()

        mol2.addAtom(c1)
        mol2.addAtom(c2)
        mol2.addAtom(c3)
        mol2.addBond(Bond(c1, c2))
        mol2.addBond(Bond(c2, c3))

        self.assertTrue(self.vf3.isSubgraphIsomorphic(mol2, mol, None))
        self.assertFalse(self.vf3.isSubgraphIsomorphic(mol, mol2, None))

    def testNumberOfMappings(self):
        examples = {"CC": 72,
                    "CCC": 144,
                    "CC[O]": 12,
                    "C1=CC=CC=C1": 6,
                    "C1(CCCC1)O[H]": 32,
                    "C1=CC=C2C(=C1)C=CC=C2": 1}

        for smiles in examples.keys():
            mol = Molecule().fromSMILES(smiles)
            self.assertEquals(len(self.vf3.findIsomorphism(mol, mol.copy(True), None)), examples[smiles])

    def testInitialMapping(self):
        mol1 = Molecule().fromSMILES("CCC")
        mol2 = mol1.copy(True)

        mol1.sortVertices()
        mol2.sortVertices()

        initMap = {mol1.vertices[0]: mol2.vertices[0], mol1.vertices[2]: mol2.vertices[2]}

        self.assertTrue(self.vf3.isIsomorphic(mol1, mol2, initMap))


###############################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
