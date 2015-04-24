import unittest
import numpy as np

from utilz import *

class RotatorTest( unittest.TestCase ):

    def setUp(self):
        pass

    def test_get_euler(self):
        v1, v2 = map(lambda x: np.array(x), [[0,0,1], [1,0,0]] )
        r1, r2, r3 = get_euler( v1, v2 )
        np.testing.assert_allclose( [r1,r2,r3], [0,0,0] )

    def test_rotate_point_by_two_points(self):
        p = np.array( [ 1,1,1] )
        p1 = np.array( [0, 0, 0] )
        p2 = np.array( [0, 0, 1] )
        theta = np.pi/2
        p_out = rotate_point_by_two_points(p, p1, p2, theta)
        np.testing.assert_allclose( p_out, np.array( [-1,1,1] ), atol=1e-10 )

    def test_rotate_point_by_two_points2(self):
        p = np.array( [ -1, -1, 0] )
        p1 = np.array( [ 0, 0,  1] )
        p2 = np.array( [ 0, 0, -1] )
        theta = np.pi * 3/2
        p_out = rotate_point_by_two_points(p, p1, p2, theta)
        np.testing.assert_allclose( p_out, np.array( [ 1, -1, 0] ), atol=1e-10 )



if __name__ == '__main__':
    unittest.main()
