import unittest, warnings 
import numpy as np

warnings.simplefilter('error')
from nose.plugins.attrib import attr
from moltools import Molecule, Water, Atom, Property, Template, Generator

@attr(speed = 'fast' )
class WaterEulerTest( unittest.TestCase ):

    def setUp(self):
        self.ut_alpha = np.random.random( (6, ) )
        self.ut_beat = np.random.random(  (10, ) )

        self.g = Generator()
        self.w = Water.get_standard()
        self.w.translate_by_r( np.random.uniform( -10, 10, [3] ) )

        self.t1 = np.random.uniform( 0, np.pi/2 )
        self.t2 = np.random.uniform( 0, np.pi   )
        self.t3 = np.random.uniform( 0, np.pi/2 )

    def test_euler(self):
        a1 = Atom( pdb_name = 'ORIGO')
        a2 = Atom( z = 1, pdb_name = 'Z')
        a3 = Atom( x = 1, pdb_name = 'XZ')
        m = Molecule( a1, a2, a3 )
        m.t( 6, 6, 6 )
        m.rotate( 0, 0, 0 )
        f = lambda x: (x.get_atom_by_pdbname('ORIGO').r,
                x.get_atom_by_pdbname('Z').r,
                x.get_atom_by_pdbname('XZ').r)
        r1, r2, r3 = m.get_euler( key = f )
        np.testing.assert_allclose( [r1, r2, r3], 0.0, atol=1e-7 )


    def test_inv_rotation(self):
        """rotate two water molecules by random angles,
        translate them in random space, one by oxygen other by center of mass
        rotate them to standard orientation and center by mass"""
        w1 = Water.get_standard()
        w2 = Water.get_standard()

        t1, t2, t3 = np.random.uniform( -np.pi, np.pi, [3] )
        d1, d2, d3 = np.random.uniform( -np.pi, np.pi, [3] )

        r1, r2 = np.random.uniform( -100, 100, [2,3] )

        w1.translate_o( r1 )
        w2.translate( r1 )

        w1.rotate( t1, t2, t3 )
        w2.rotate( d1, d2, d3 )

        w1.inv_rotate( t3, t2, t1 )
        w2.inv_rotate( d3, d2, d1 )

        w1.center()
        w2.center()

        self.eq ( w1.com, w2.com )

        
    def test_center_get_euler(self):
        w = Water.get_standard()

        t1, t2, t3 = w.get_euler( lambda x: (x.o.r, (x.h1.r-x.h2.r)/2.0 + x.h2.r, x.h1.r ))
        self.eq( [t1, t2, t3] , np.zeros(3) )

    def test_moved_get_euler(self):
        w = Water.get_standard()
        w.translate([5,5,5])
        t1, t2 ,t3 =  w.get_euler( lambda x: (x.o.r, (x.h1.r-x.h2.r)/2.0 + x.h2.r, x.h1.r ))
        w.center()

        self.eq( [t1, t2, t3] , np.zeros(3) )

        t1, t2 ,t3 =  w.get_euler(lambda x: (x.o.r, (x.h1.r-x.h2.r)/2.0 + x.h2.r, x.h1.r ))
        w.center()

        self.eq( [t1, t2, t3] , np.zeros(3) )

    def eq(self, a, b):
        np.testing.assert_almost_equal( a, b, decimal = 3)



if __name__ == '__main__':
    unittest.main()
