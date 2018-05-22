from rmgpy.molecule.vf2Test import TestVF2
from rmgpy.molecule.vf3 import VF3


class TestVF3(TestVF2):
    """
    Applies the tests for the old VF2 algorithm to the
    updated VF3 implementation
    """

    def setUp(self):
        super(TestVF3, self).setUp()
        self.vf2 = VF3()
