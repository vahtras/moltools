import unittest, mock, os
import numpy as np

from use_generator import Generator
from molecules import Water, Molecule, Methanol


FILE_XYZ =os.path.join( os.path.dirname( os.path.realpath( __file__ ) ), "tmp.xyz" )


class GeneratorTestCase( unittest.TestCase ):

    def setUp(self):
        self.g = Generator()

    def test_get_hfqua_dal(self):
        assert "**DALTON" in self.g.get_hfqua_dal()
        assert ".PARALLELL" in self.g.get_hfqua_dal()
        assert ".DIPLEN" in self.g.get_hfqua_dal()
        assert "**WAVE FUNCTION" in self.g.get_hfqua_dal()
        assert ".HF" in self.g.get_hfqua_dal()
        assert "**END" in self.g.get_hfqua_dal()

    def test_get_mol(self):
        w = self.g.get_mol( center = [0, 1, 1], mol = "water", AA = False )
        m = self.g.get_mol( center = [0, 1, 1], mol = "methanol", AA = False )
        self.assertIsInstance( w, Molecule )
        self.assertIsInstance( m, Molecule )
        self.assertNotIsInstance( w, Methanol )
        self.assertNotIsInstance( m, Water )

    def test_vary_parameters(self):
#Could modify in future this test todo
        a = self.g.vary_parameters( {"dog" : {"min":"small", "max":5, "points":33} } )
        self.assertIsNone(a)

    def test_polar_to_cartesian_5_0_0(self):
        v = self.g.polar_to_cartesian( 5, 0, 0 )
        v_ref = np.zeros(3)
        v_ref[2] = 5
        self.eq( v, v_ref )

    def test_polar_to_cartesian_3_pi_pi(self):
        v = self.g.polar_to_cartesian( 3, np.pi, np.pi )
        v_ref = np.zeros(3)
        v_ref[2] = -3
        self.almost_eq( v, v_ref, decimal = 14 )

#Add mock here todo for testing
    #def test_build_pna(self):
    #    d = FILE_XYZ
    #    res = self.g.build_pna( xyz = d )
    #    self.assertIsNone( res )

    def almost_eq(self, a, b, decimal = 7):
        np.testing.assert_almost_equal( a, b, decimal = decimal )

    def eq(self, a, b, ):
        np.testing.assert_equal( a, b )

    def tearDown(self):
        files = [f for f in os.listdir( os.path.dirname( os.path.realpath(__file__) ) ) if f.endswith(".pot") or f == "pna.mol" or f.startswith("5.00-") ]
        for fi in files:
            os.remove( os.path.join(os.path.dirname( __file__ ) , fi) )


if __name__ == '__main__':
    unittest.main()
