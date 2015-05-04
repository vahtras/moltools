import unittest, os
import numpy as np

import utilz
from molecules import Property, Water

class WaterTest( unittest.TestCase ):

    def setUp(self):
        self.p1 = Property()
        self.p1['charge'] = 0.22
        self.p1['dipole'] = np.array( [0.0, 0.0, -0.5 ] )
        self.p1['quadrupole'] = np.arange( 6 )
        self.p1['alpha'] = np.arange( 6, 0, -1 )
        self.p1['beta'] = np.arange( 10 )

    def test_potline_0_0_0(self):
        assert self.p1.potline(max_l = 0, pol = 0, hyper = 0) == "0.22000 "

    def test_potline_1_0_0(self):
        assert self.p1.potline(max_l = 1, pol = 0, hyper = 0) == \
                "0.22000 0.00000 0.00000 -0.50000 "

    def test_potline_2_0_0(self):
        assert self.p1.potline(max_l = 2, pol = 0, hyper = 0) == \
                "0.22000 0.00000 0.00000 -0.50000 0.00000 1.00000 2.00000 3.00000 4.00000 5.00000 "

    def test_potline_0_1_0(self):
        print self.p1.potline(max_l = 0, pol = 1, hyper = 0)
        assert self.p1.potline(max_l = 0, pol = 1, hyper = 0) == \
                "0.22000 3.33333 "

    def test_potline_0_2_0(self):
        assert self.p1.potline(max_l = 0, pol = 2, hyper = 0) == \
                "0.22000 6.00000 5.00000 4.00000 3.00000 2.00000 1.00000 "
    def test_potline_0_22_1(self):
        print self.p1.potline(max_l = 0, pol = 22, hyper = 1) 
        assert self.p1.potline(max_l = 0, pol = 22, hyper = 1) == \
                "0.22000 6.00000 5.00000 4.00000 3.00000 2.00000 1.00000 " + \
                "0.00000 1.00000 2.00000 3.00000 4.00000 5.00000 6.00000 7.00000 8.00000 9.00000 "


    def test_property_beta_proj(self):
        w = Water.get_standard() 
        w.translate_by_r( np.random.uniform( -10, 10, (3,) ))
        w.rotate( *np.random.uniform( -10, 10, (3,) ))
        w[0].p.d = np.random.random( 3 )
        w[0].p.b = np.random.random( 10 )
        w[1].p.d = np.random.random( 3 )
        w[1].p.b = np.random.random( 10 )
        w[2].p.d = np.random.random( 3 )
        w[2].p.b = np.random.random( 10 )
        BP = w.p.b_proj
        t, r1, r2, r3 = utilz.center_and_xz( w[0].r, w[0].r + (w.h1.r-w.h2.r)/2, w.h1.r )
        w.translate_by_r( t )
        w.inv_rotate( r1, r2, r3)
        np.testing.assert_allclose( BP, w.p.b_proj, atol = 1e-7 )

    def test_property_add(self):
        p1 = Property()
        p1["charge"] = 0.3
        p2 = Property()
        p2["charge"] = 0.7
        assert (p1 + p2)["charge"] == 1.0

    def test_property_sub(self):
        p1 = Property()
        p1["charge"] = 0.3
        p2 = Property()
        p2["charge"] = 0.7
        self.eq( (p2 - p1)["charge"] , 0.4, decimal = 7 )

    def test_inv_rotate(self):
        p = Property()
        p.d = np.array( [ 0, 0, 1] )
        p.a = np.array( [ 3, 0, 0, 2, 0, 7 ] )
        p.Q = np.array( [ 3, 0, 0, 2, 0, 7 ] )
        p.b = np.array( [ 3, 0, 0, 0, 0, 0, 5, 0, 0, 10 ] )

        P = p.inv_rotate( 0, np.pi , 0 )

        np.testing.assert_allclose( P.d, np.array([0, 0, -1]) , atol = 1e-7 )
        np.testing.assert_allclose( P.a, np.array([3, 0, 0, 2, 0, 7]) , atol = 1e-7 )
        np.testing.assert_allclose( P.Q, np.array([3, 0, 0, 2, 0, 7]) , atol = 1e-7 )
        #np.testing.assert_allclose( P.b, np.array([3, 0, 0, 0, 0, 0, 5, 0, 0, 10 ]), atol = 1e-7 )



    def eq(self, a, b, decimal = 7):
        np.testing.assert_almost_equal( a, b, decimal = decimal )

if __name__ == '__main__':
    unittest.main()
