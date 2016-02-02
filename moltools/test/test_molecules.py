import unittest, os
import numpy as np

import warnings
warnings.simplefilter('error')
from nose.plugins.attrib import attr
from moltools import Water, Molecule, Property, Template, Generator, unique

@attr(speed = 'fast' )
class MoleculesTestCase( unittest.TestCase ):
    def setUp(self):
        wat = Water.get_standard()
        t1, t2, t3 = wat.get_euler()
        kw_dict = Template().get( *("TIP3P", "HF", "ANOPVDZ",
            True, "0.0"))
        for at in wat:
            at.p = Property.from_template( at.name, kw_dict )
            at.p = at.p.rotate( t1, t2, t3 )
# will only test this quadrupole
            at.Property["quadrupole"] = np.arange( 6 )
        self.wat = wat

#    @mock.patch( "use_generator.open" ,create = True)
#    def test_gen_mols_param(self, mock_open):
#        mock_open.return_value = mock.MagicMock( spec = file )
#        ret = self.g.gen_mols_param()
#        assert ret == 0

    def test_dist_to_mol(self):
        wat1, wat2 = Water.get_standard(), Water.get_standard()
        wat1.translate( [0,0,0] )
        wat2.translate( [0,0,1] )
        self.eq(  np.linalg.norm( wat1.dist_to_mol( wat2 )) , 1, decimal = 7 )

        wat1.translate( [0, 0, 0] )
        wat2.translate( [1, 1, 1] )
        self.eq(  np.linalg.norm( wat1.dist_to_mol( wat2 )) ,np.sqrt(3) , decimal = 7 )

# Ignore the opening of the file, just assert its a correct type returned
#    @mock.patch( "molecules.open", create = True )
#    def test_from_mol(self, mock_open):
#        mock_open.return_value = mock.MagicMock( spec = file )
#        w = Water.from_mol_file( "tip3p.mol" )
#        self.assertIsInstance( w, Molecule )


    def test_get_xyz_string(self):
        st = """3

O            0.000000   0.000000   0.000000
H            0.756950   0.000000   0.585882
H           -0.756950   0.000000   0.585882
"""
        self.wat.to_AA()
        self.assertEqual( st, self.wat.get_xyz_string() )
        self.wat.to_AU()

    def test_get_mol_string(self):
        st = """ATOMBASIS


Atomtypes=2 Charge=0 Nosymm Angstrom
Charge=1.0 Atoms=2 Basis=ano-1 2 0 0 0
0-MOL-HW2        -0.75695   0.00000   0.58588
0-MOL-HW1         0.75695   0.00000   0.58588
Charge=8.0 Atoms=1 Basis=ano-1 4 3 1 0
0-MOL-OW          0.00000   0.00000   0.00000
"""
        self.wat.to_AA()
        assert self.wat.get_mol_string() == st 
        self.wat.to_AU()


    def test_translate_coc(self):
        self.wat.translate_coc( [0, 0, 0] )
        self.eq( self.wat.coc , [0, 0, 0] )

        self.wat.translate_coc( [1, 2, 3] )
        self.eq( self.wat.coc , [1, 2, 3] )

    def tearDown(self):
        files = [f for f in os.listdir( os.path.dirname( os.path.realpath(__file__) ) ) if f.endswith(".pot") or f == "pna.mol" or f.startswith("5.00-") ]
        for fi in files:
            os.remove( os.path.join(os.path.dirname( __file__ ) , fi) )

    def eq(self, a, b, decimal = 7):
        np.testing.assert_almost_equal( a, b, decimal = decimal )


if __name__ == '__main__':
    unittest.main()
