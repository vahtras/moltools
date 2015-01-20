import unittest, mock, os
import numpy as np

from use_generator import Generator
from molecules import Water, Molecule, Methanol, Property
from template import Template

class MoleculesTestCase( unittest.TestCase ):
    def setUp(self):
        wat = Generator().get_mol( AA = False )
        t1, t2, t3 = wat.get_euler()
        kw_dict = Template().get( *("TIP3P", "HF", "ANOPVDZ",
            True, "0.0"))
        for at in wat:
            Property.add_prop_from_template( at, kw_dict )
            at.Property.transform_ut_properties( t1, t2, t3 )
# will only test this quadrupole
            at.Property["quadrupole"] = np.arange( 6 )
        self.wat = wat

#    @mock.patch( "use_generator.open" ,create = True)
#    def test_gen_mols_param(self, mock_open):
#        mock_open.return_value = mock.MagicMock( spec = file )
#        ret = self.g.gen_mols_param()
#        assert ret == 0

    def test_dist_to_mol(self):
        wat1 = Generator().get_mol( [0,0,0] )
        wat2 = Generator().get_mol( [0,0,1] )
        self.eq(  np.linalg.norm( wat1.dist_to_mol( wat2 )) , 1, decimal = 7 )

        wat1.translate( [0, 0, 0] )
        wat2.translate( [1, 1, 1] )
        self.eq(  np.linalg.norm( wat1.dist_to_mol( wat2 )) ,np.sqrt(3) , decimal = 7 )

# Ignore the opening of the file, just assert its a correct type returned
    @mock.patch( "molecules.open", create = True )
    def test_from_mol(self, mock_open):
        mock_open.return_value = mock.MagicMock( spec = file )
        w = Water.from_mol_file( "tip3p.mol" )
        self.assertIsInstance( w, Molecule )


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
Charge=8.0 Atoms=1 Basis=ano-1 3 2 1
O 0.00000 0.00000 0.00000
Charge=1.0 Atoms=2 Basis=ano-1 2 1
H 0.75695 0.00000 0.58588
H -0.75695 0.00000 0.58588
"""
        self.wat.to_AA()
        assert self.wat.get_mol_string() == st 
        self.wat.to_AU()


    def test_translate_coc(self):
        self.wat.translate_coc( [0, 0, 0] )
        self.eq( self.wat.coc , [0, 0, 0] )

        self.wat.translate_coc( [1, 2, 3] )
        self.eq( self.wat.coc , [1, 2, 3] )

    def test_transform_dist_quadrupole(self):
        print self.wat.o.Property["quadrupole"]

    def tearDown(self):
        files = [f for f in os.listdir( os.path.dirname( os.path.realpath(__file__) ) ) if f.endswith(".pot") or f == "pna.mol" or f.startswith("5.00-") ]
        for fi in files:
            os.remove( os.path.join(os.path.dirname( __file__ ) , fi) )

    def eq(self, a, b, decimal = 7):
        np.testing.assert_almost_equal( a, b, decimal = decimal )


if __name__ == '__main__':
    unittest.main()
