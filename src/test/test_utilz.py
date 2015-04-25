import unittest
import numpy as np

from utilz import *

class RotatorTest( unittest.TestCase ):

    def setUp(self):
        pass

    def test_center_and_xz_1(self):
        p1 = np.array( [5, 5, 5] )
        p2 = np.array( [5, 5, 6] )
        p3 = np.array( [6, 5, 5] )
        t_v, r1, r2, r3 = center_and_xz( p1, p2, p3 )
        np.testing.assert_allclose( t_v, np.full( (3,), -5 ) , atol =1e-10 )
        np.testing.assert_allclose( r1, 0.0 , atol =1e-10 )
        np.testing.assert_allclose( r2, 0.0 , atol =1e-10 )
        np.testing.assert_allclose( r3, 0.0 , atol =1e-10 )


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

    def test_rotate_point_by_two_points_2(self):
        p = np.array( [ -1, -1, 0] )
        p1 = np.array( [ 0, 0,  1] )
        p2 = np.array( [ 0, 0, -1] )
        theta = np.pi * 3/2
        p_out = rotate_point_by_two_points(p, p1, p2, theta)
        np.testing.assert_allclose( p_out, np.array( [ 1, -1, 0] ), atol=1e-10 )

    def test_rotate_point_by_two_points_3(self):
        p = np.array( [ 0, 3, 1] )
        p1 = np.array( [ 3, 3,  0] )
        p2 = np.array( [ 3, 3,  1] )
        theta = np.pi * 3/2
        p_out = rotate_point_by_two_points(p, p1, p2, theta)
        np.testing.assert_allclose( p_out, np.array( [ 3, 6, 1] ), atol=1e-10 )




if __name__ == '__main__':
    unittest.main()
