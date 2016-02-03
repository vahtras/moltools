import unittest, os
import numpy as np


import warnings
warnings.simplefilter('error')
from nose.plugins.attrib import attr

from moltools import Atom, Molecule, Water, Cluster, Property, Generator, Residue, System
from moltools import utilz 

FILE_XYZ = os.path.join( os.path.dirname(__file__), 'pna_waters.xyz' )
FILE_MOL = os.path.join( os.path.dirname(__file__), 'tip3p44_10qm.mol' )
FILE_PDB = os.path.join( os.path.dirname(__file__), 'tip3p0.pdb' )

@attr(speed = 'fast' )
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
        w2.LoProp = True
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
        w = Water.get_standard()
        self.assertNotIn( w, self.c ) 
        self.c.add_mol( w )
        self.assertEqual( len( self.c ), (a + 1) )
        self.assertIn( w, self.c ) 

    def test_copy_cluster(self):
        c2 = self.c.copy_cluster()
        #self.assertNotEqual( c2, self.c )
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

    def test_sum_property(self):
#Make sure that cluster total dipole moment is same weather molecular units are LoProp or Simple mode
        a1 = Atom(element='H', z = 1.   )
        a2 = Atom(element='H', z = 1.   )
        a3 = Atom(element='O', z = 0.0  )
        a1.p.q = 0.25
        a2.p.q = 0.25
        a3.p.q = -0.5
        a1.p.d = np.array( [ 0., 0., -0.2 ] )
        a2.p.d = np.array( [ 0., 0., -0.2 ] )
        a3.p.d = np.array( [ 0., 0., 0.5 ] )
        m1 = Molecule( a1, a2, a3 )
        np.testing.assert_allclose( m1.p.d, np.array([0.,0.,0.6]), atol = 1e-7 )

        m2 = m1.copy()
        m2.Property = Property()
        m2.p.d = np.array( [0., 0., 0.6] )
        c = Cluster(m1, m2 )
        np.testing.assert_allclose( c.p.d, np.array([0.,0.,1.2]), atol = 1e-7 )

        m1.to_Property()
        np.testing.assert_allclose( c.p.d, np.array([0.,0.,1.2]), atol = 1e-7 )

    def test_both_loprop_and_none(self):
        w1, w2, w3, w4 = [Water.get_standard() for i in range(4)]
        c1 = Cluster( w1, w2.t(0, 0, 5) )
        c2 = Cluster( w3, w4.t(0, 0, 5) )

        c1.attach_properties( loprop = False )
        c2.attach_properties( loprop = True )

        np.testing.assert_allclose( c1.p.d, c2.p.d )

    def eq(self, a, b):
        np.testing.assert_almost_equal( a, b, decimal = 3)

    def test_attach_both_ways(self):
        c1 = Cluster( [Water.get_standard() for i in range(2 )] )
        c2 = Cluster( [Water.get_standard() for i in range(2 )] )
        
        for w1, w2 in zip(c1, c2):
            t_vec, r_mat = np.random.random( 3 ), np.random.random( 3 )
            w1.rotate( *r_mat )
            w2.rotate( *r_mat )
            w1.t( t_vec )
            w2.t( t_vec )

        np.testing.assert_allclose( w1.o.r, w2.o.r )

        c1.attach_properties( model = 'spc', force_template = True )
        for wat in c2:
            wat.attach_properties( model = 'spc', force_template = True )

        np.testing.assert_allclose( c1.p.a, c2.p.a )


    def test_system_and_cluster( self ):
        sys = System.from_pdb_string( open(FILE_PDB).read() )
        clus = Cluster.get_water_cluster( FILE_PDB, in_AA = True,
                out_AA = True,
                N_waters = 1000)
        np.testing.assert_allclose( sys.coc , clus.coc )



if __name__ == '__main__':
    unittest.main()
