import unittest, os
import numpy as np

import warnings
warnings.simplefilter('error')
from nose.plugins.attrib import attr
from moltools import Generator, Water, Molecule


FILE_XYZ =os.path.join( os.path.dirname( os.path.realpath( __file__ ) ), "tmp.xyz" )


@attr(speed = 'fast' )
class GeneratorTestCase( unittest.TestCase ):

    def setUp(self):
        self.g = Generator()

    def test_get_hfqua_dal(self):
        string = self.g.get_hfqua_dal()
        assert "**DALTON" in string
        assert ".PARALLELL" in string
        assert ".DIPLEN" in string
        assert "**WAVE FUNCTION" in string
        assert ".HF" in string
        assert "**END" in string

    def test_get_mol(self):
        w = Water.get_standard() 
        self.assertIsInstance( w, Molecule )

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

    #@mock.patch( "use_generator.open" ,create = True)
    #def test_gen_mols_param(self, mock_open):
    #    mock_open.return_value = mock.MagicMock( spec = file )
    #    ret = self.g.gen_mols_param()
    #    assert ret == 0

#Mocking added to ensure that the file tmp.xyz isn't opened writable during tests
#    @mock.patch( "use_generator.write" ,create = True)
#    @mock.patch( "use_generator.open" ,create = True)
#    def test_build_pna(self, mock_open, mock_write):
#        mock_open.return_value = mock.MagicMock( spec = file )
#        mock_write.return_value = mock.MagicMock( spec = file )
#        d = FILE_XYZ
#        res = self.g.build_pna( xyz = d, waters = 1 )
#        self.assertIsNone( res )

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
