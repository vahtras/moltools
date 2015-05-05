import unittest, os
import numpy as np

from molecules import Cluster, Atom, Water
from use_generator import Generator

FILE_XYZ = os.path.join( os.path.dirname(__file__), 'pna_waters.xyz' )
FILE_MOL = os.path.join( os.path.dirname(__file__), 'tip3p44_10qm.mol' )
FILE_PDB = os.path.join( os.path.dirname(__file__), 'tip3p0.pdb' )

class WaterTest( unittest.TestCase ):
    def setUp(self):
        self.c = Cluster.get_water_cluster(
                FILE_MOL,
                in_AA = False,
                out_AA = False,
                N_waters = 10)

    def test_dipole(self):
        w1 = Water.get_standard()
        w2 = Water.get_standard()
        w2.o.p.q -= 1.0
        w2.translate_by_r( np.array([ 3, 3, 3] ) )

        c = Cluster( w1, w2 )
        d1 = c.p.d.copy()
        r = np.random.uniform(-10,10, (3,) ) 
        for i in c:
            i.translate_by_r( r )

        np.testing.assert_allclose( d1, c.p.d, atol =1e-7 )

    def test_size(self):
        assert len( self.c) == 10

    def test_min_len(self):
        assert len( self.c.min_dist_coo() ) == 9

    def test_dists(self):
        print len( self.c )
        print self.c.min_dist_coo()[:len(self.c)]

    def test_two_waters(self):
        c = Cluster.get_water_cluster(
                FILE_MOL,
                in_AA = False,
                out_AA = False,
                N_waters = 2)
        assert len(c) == 2
        print c.min_dist_coo()
        assert c.min_dist_coo()[0] > 4.5

    def test_get_all_molecules_from_file( self ):
        c = Cluster.get_all_molecules_from_file( FILE_XYZ,
                in_AA = False,
                out_AA = False,
                )

    def test_get_water_cluster(self):
        c = Cluster.get_water_cluster( FILE_PDB, in_AA = True, N_waters = 10 )
        assert len(c) == 10

    def test_add_atom(self):
        a = Atom( {'element':"H"} )
        self.assertNotIn( a, self.c ) 
        self.c.add_atom( a )
        self.assertIn( a, self.c )

    def test_add_mol(self):
        a = len(self.c)
        w = Generator().get_mol() 
        self.assertNotIn( w, self.c ) 
        self.c.add_mol( w )
        self.assertEqual( len( self.c ), (a + 1) )
        self.assertIn( w, self.c ) 

    def test_copy_cluster(self):
        c2 = self.c.copy_cluster()
        self.assertNotEqual( c2, self.c )
        self.assertEqual( len(self.c), len(c2) )

    def test_set_qm_mm(self):
        self.c.set_qm_mm( N_qm = 4, N_mm = 2 )
        self.assertEqual( len([i for i in self.c if i.in_qm ]), 4 )
        self.assertEqual( len([i for i in self.c if i.in_mm ]), 2 )

    def test_three_waters(self):
        c = Cluster.get_water_cluster(
                FILE_MOL,
                in_AA = False,
                out_AA = False,
                N_waters = 3)
        assert len(c) == 3
        assert c.min_dist_coo()[0] > 4.5

        dists = c.min_dist_coo()

        assert np.mean( [ dists[0], dists[1] ] ) == \
                (dists[0] + dists[1])/2


    def eq(self, a, b):
        np.testing.assert_almost_equal( a, b, decimal = 3)



if __name__ == '__main__':
    unittest.main()
