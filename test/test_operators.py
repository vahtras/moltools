
import unittest
import numpy as np

from use_generator import Generator
from template import Template
from molecules import Water, Atom, Property, Molecule, Rotator

class WaterTest( unittest.TestCase ):

    def setUp(self):
        self.ut_alpha = np.random.random( (6, ) )
        self.ut_beat = np.random.random(  (10, ) )

        self.g = Generator()
        self.w = self.g.get_mol( center = np.random.uniform( -10, 10, [3] ), mol = "water" )

        self.t1 = np.random.uniform( 0, np.pi/2 )
        self.t2 = np.random.uniform( 0, np.pi   )
        self.t3 = np.random.uniform( 0, np.pi/2 )


    def test_negative_y_get_euler(self):

        w = self.g.get_mol( center = [0,0,0], mol = "water" )

        t1 = 0
        t2 = 0
        t3 = 0

        t1 = np.pi/2
        #t2 = np.pi/2
        #t3 = np.pi/2

        #t1, t2, t3 = w.get_euler()
        Rz1 = Rotator.get_Rz( t1 )
        #Ry = Molecule.get_Ry_inv( t2 )
        #Rz2 = Molecule.get_Rz( t3 )



        #assert isinstance( w,  )



    def eq(self, a, b):
        np.testing.assert_almost_equal( a, b, decimal = 3)



if __name__ == '__main__':
    unittest.main()
