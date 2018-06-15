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

from rmgpy.molecule.vf2Test import TestVF2
from rmgpy.molecule.molecule import Molecule
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

        for example in examples:
            mol = Molecule().fromSMILES(example)
            self.assertTrue(self.vf3.isIsomorphic(mol, mol.copy(True), None))



    def testSubgraphIsomorphism(self):
        pass#TODO


###############################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
