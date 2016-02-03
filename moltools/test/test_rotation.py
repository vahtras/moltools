import unittest
import numpy as np
from numpy.linalg import norm

import warnings
warnings.simplefilter('error')
from nose.plugins.attrib import attr

from moltools import utilz
from moltools import Generator, Template, Water, Atom, Property

@attr(speed = 'fast' )
class WaterTest( unittest.TestCase ):

    def setUp(self):
        np.random.seed(111)
        self.ut_alpha = np.random.random( (6, ) )
        self.ut_beat = np.random.random(  (10, ) )

        self.g = Generator()
        r = np.random.uniform( -10, 10, [3] )
        self.w = Water.get_standard()
        self.w.translate_o( r )

        t1 = np.random.uniform( 0, np.pi )
        t2 = np.random.uniform( 0, np.pi )
        t3 = np.random.uniform( 0, np.pi )
        self.w.rotate( t1, t2, t3 )

        assert len( self.w ) == 3

    #def test_random_rotation(self):
    #    t1, t2, t3 = self.w.get_euler
    #    w = self.g.get_mol( center = np.random.uniform( -1,1,[3] ), mol= "water" )
    #    w.rotate( t1, t2, t3 )
    #    w.translate( self.w.r )

    #    self.eq( w.com, self.w.com )


    def test_translation(self):
        w = Water.get_standard()
        w.translate( self.w.o.r )
        self.eq( w.com , self.w.o.r )

        w.translate( [1,1,1] )
        w.translate( [0,0,0] )
        w.translate( [1,1,1] )
        w.translate( [1,1,1] )
        w.translate( [2333,1,1] )
        w.translate( [0,0,0] )

        self.eq( w.com, [0,0,0], decimal = 7 )


    def test_center_water(self):
        w = Water.get_standard()

#Move the water to random place
        w.translate( np.random.uniform(-100,100, [3]  ))
#Call center function

        w.center()
#Assert that oxygen is in origo

        self.eq( w.com, [0,0,0], decimal = 7 )


    def test_rotate_around_z(self):
        w = Water.get_standard()
        w.translate( self.w.com )
        self.eq( w.com, self.w.com )

        origin = w.o.r.copy()
        w.translate_by_r( -w.o.r )
        w.rotate( np.pi/2 , 0 , np.pi/2 )
        w.translate_by_r( origin )
        self.eq( w.com, self.w.com )

    def test_rotate_to_random(self):
        w = Water.get_standard( )
        w.translate( self.w.com )
        t1, t2, t3 = self.w.get_euler()
        origin = w.com.copy()
        w.translate_by_r( -w.com )
        w.rotate( t1, t2, t3 )
        w.translate_by_r( origin )
        self.eq( w.com, self.w.com )

    def test_rotate_inv_Z(self):
        w = Water.get_standard()

        self.eq( w.h1.y, 0 )
        self.eq( w.h2.y, 0 )
        self.eq( w.o.r, [0,0,0] )

#rotate around z axis by 90 degree, both hydrogens now in zy plane
        w.rotate( np.pi/2 ,0 ,0)
        self.eq( w.h1.x, 0 )
        self.eq( w.h2.x, 0 )
        self.eq( w.o.r, [0,0,0] )

#rotate around z axis by 90 degree
        #w.rotate( np.pi/2 ,0 ,0)
        #self.eq( w.o.r, [0,0,0] )
        #self.eq( w.h1.y, 0 )
        #self.eq( w.h2.y, 0 )
    def test_get_euler(self):
        pass

    def test_rotate_dip(self):
        w = Water.get_standard( )
        w.translate_o( np.zeros(3) )

        for at in w:
            at.p = Property.from_template( at.name, Template().get( model = 'TIP3p',
                method ='HF',
                basis = 'ANOPVDZ') )
        w.LoProp = True

        dip = 0.7870603

        np.testing.assert_allclose( w.p.d ,[0, 0, dip ], atol = 1e-7)

        w.rotate( np.pi/2, 0, 0 )

        self.eq( w.p.d ,[0, 0, dip])

        print "After 90 around Z-axis: counter-clock"
        print w.p
        w.rotate( 0, np.pi/2 , 0 )
        self.eq( w.p.d ,[-dip,  0, 0])

        print "After 90 around Y-axis: clock"
        print w.p.d

        w.rotate( 0, 0, np.pi/2 )
        print "After 90 around Z-axis: counter-clock"
        print w.p.d
        t, r1, r2, r3 = utilz.center_and_xz(w.o.r, w.o.r+(w.h1.r-w.h2.r)/2,w.h1.r )
        w.inv_rotate(r1, r2, r3)
        #for at in w:
        #    Property.add_prop_from_template( at, Template().get() )
        #w.rotate( np.pi/2, np.pi/2, np.pi/2 )
        #print w.p
        #self.eq( w.p ,[0, -0.78704747, 0])

        #w.rotate( 0, np.pi/2, 0 )
        #self.eq( w.p ,[0, 0,-0.41021007, ])

    def test_y_rotation(self):
        w = Water.get_standard( )
#rotate around y axis by 90 degree, assert

        w.translate( self.w.o.r )
        w.rotate( 0 , np.pi/2 ,0)
        w.rotate( 0 , np.pi/2 ,0)
        w.center()

        self.eq( w.h1.y, 0 )
        self.eq( w.h2.y, 0 )

        #self.eq( np.arccos( np.dot(w.p, np.array([0,0,1])) /(norm(w.p)*norm([0,0,1])) , 0  ))

    def test_ut_2_to_square(self):
        assert len( self.ut_alpha ) == 6

    def test_(self):
        assert len( self.ut_alpha ) == 6

    def eq(self, a, b, decimal = 4):
        np.testing.assert_almost_equal( a, b, decimal = decimal)

if __name__ == '__main__':
    unittest.main()
