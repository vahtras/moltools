import unittest
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
        assert len( self.w ) == 3

    def test_translation(self):
        w = self.g.get_mol( center = [0,0,0], mol = "water" )
        w.translate( self.w.o.r )
        assert np.equal( w.o.r , self.w.o.r ).all()



    def test_center_water(self):
        w = self.g.get_mol( center = [0,0,0], mol = "water" )

#Move the water to random place
        print w.o.r
        w.translate( np.random.uniform(-100,100, [3]  ))
        print w.o.r
#Call center function

        w.center()
        print w.o.r
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
        self.eq( w.o.r, [0,0,0] )
        #self.eq( w.h2.y, 0 )
        #self.eq( w.h1.y, 0 )
    

    def test_y_rotation(self):
        w = self.g.get_mol( center = [0,0,0], mol = "water" )
#rotate around y axis by 90 degree, assert
        w.rotate( 0 , np.pi/2 ,0)
        #self.eq( w.p, [-1, 0, 0] )

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
