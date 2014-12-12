import unittest, subprocess, os
import numpy as np

from use_generator import Generator
from template import Template
from molecules import Water, Atom, Property


class WaterTest( unittest.TestCase ):

    def setUp(self):
        self.ut_alpha = np.random.random( (6, ) )
        self.ut_beat = np.random.random(  (10, ) )

        self.g = Generator()
        self.w = self.g.get_mol( center = np.random.uniform( -10, 10, [3] ), mol = "water" )

        self.t1 = np.random.uniform( 0, np.pi/2 )
        self.t2 = np.random.uniform( 0, np.pi   )
        self.t3 = np.random.uniform( 0, np.pi/2 )

    def test_dalton_output_same_as_loprop(self):
        w = self.g.get_mol( center = [0,0,0], mol = "water", AA = False )
        val = "/tmp"
        w.rotate( 1, 1, 1 )
### Add dalton output from a rotated water by 1,1,1 and assert all components identical to transform_ut_10 from loprop template

        assert 2==2

    def eq(self, a, b, decimal = 7):
        np.testing.assert_almost_equal( a, b, decimal = decimal )



if __name__ == '__main__':
    unittest.main()
