import unittest, os
import numpy as np

from use_generator import Generator

class UseGeneratorTest( unittest.TestCase ):

    def setUp(self):
        self.g = Generator()
    def test_get_hfqua_dal(self):
        assert ".PARALLELL" in self.g.get_hfqua_dal()
        assert ".DIPLEN" in self.g.get_hfqua_dal()
        assert "**WAVE FUNCTION" in self.g.get_hfqua_dal()

    #def test_gen_mols_param( self ):
    #    assert self.g.gen_mols_param() == 0

    def eq(self, a, b, decimal = 7):
        np.testing.assert_almost_equal( a, b, decimal = decimal )

if __name__ == '__main__':
    unittest.main()
