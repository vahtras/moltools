import unittest
import numpy as np
from numpy.linalg import norm

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
        assert len( self.w ) == 3

    def test_translation(self):
        w = self.g.get_mol( center = [0,0,0], mol = "water" )

        w.translate( [1,1,1] )
        w.translate( [0,0,0] )
        w.translate( [1,1,1] )
        w.translate( [1,1,1] )
        w.translate( [2333,1,1] )
        w.translate( [0,0,0] )

        assert np.equal( w.o.r , [0,0,0] ).all()


    def test_center_water(self):
        w = self.g.get_mol( center = [0,0,0], mol = "water" )

#Move the water to random place
        w.translate( np.random.uniform(-100,100, [3]  ))
#Call center function
        w.center()
#Assert that oxygen is in origo

        assert np.equal( w.o.r, [0,0,0] ).all()

    def test_get_euler(self):
        w = self.g.get_mol( center = [0,0,0], mol = "water" )


    def test_rotate_inv_Z(self):
        w = self.g.get_mol( center = [0,0,0], mol = "water" )

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
        w = self.g.get_mol( center = [0,0,0], mol = "water" )

        self.eq( w.p ,[0, 0, 0.41655264])

        w.rotate( np.pi/7, 0, 0 )
        self.eq( w.p ,[0, 0, 0.41655264])
        print "After 90 around Z-axis: counter-clock"
        print w.p
        w.rotate( 0, np.pi/2 , 0 )
        self.eq( w.p ,[-0.41655264, 0, 0])

        print "After 90 around Y-axis: clock"
        print w.p

        w.rotate( 0, 0, np.pi/2 )
        print "After 90 around Z-axis: counter-clock"
        print w.p
        w.rotate( np.pi/2, np.pi/2, np.pi/2 )
        print w.p
        self.eq( w.p ,[0, -0.41655264, 0])
        #assert 2==3

        #w.rotate( 0, np.pi/2, 0 )
        #self.eq( w.p ,[0, 0,-0.41655264, ])

    def test_y_rotation(self):
        w = self.g.get_mol( center = [0,0,0], mol = "water" )
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

    def eq(self, a, b):
        np.testing.assert_almost_equal( a, b, decimal = 3)

    def eq(self, a, b):
        np.testing.assert_almost_equal( a, b, decimal = 3)

if __name__ == '__main__':
    unittest.main()
