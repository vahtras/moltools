
import unittest
import numpy as np

from molecules import Water

class ReflectionTest( unittest.TestCase ):

    def setUp(self):
        self.w = Water.get_standard()
        self.w.attach_properties()

    def test_reflection(self):
        q = self.w.p.q
        d = self.w.p.d.copy()
        Q = self.w.p.Q.copy()
        a = self.w.p.a.copy()
        b = self.w.p.b.copy()
        self.w.rotate( 0, np.pi/2, 0, )
        self.w.reflect( lambda x: map( np.array,((0,0,0),(0,1,0),(0,0,1))))
        self.w.rotate( 0, np.pi/2, 0, )
        np.testing.assert_allclose( q, self.w.p.q, atol = 1e-7 )
        np.testing.assert_allclose( d, self.w.p.d, atol = 1e-7 )
        np.testing.assert_allclose( Q, self.w.p.Q, atol = 1e-7 )
        np.testing.assert_allclose( a, self.w.p.a, atol = 1e-7 )
        np.testing.assert_allclose( b, self.w.p.b, atol = 1e-7 )

if __name__ == '__main__':
    unittest.main()
