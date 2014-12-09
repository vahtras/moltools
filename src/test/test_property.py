import unittest, os
import numpy as np

from molecules import Property

class WaterTest( unittest.TestCase ):

    def setUp(self):
        self.p1 = Property()
        self.p1['charge'] = np.array( [0.22 ] )
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

    def eq(self, a, b, decimal = 7):
        np.testing.assert_almost_equal( a, b, decimal = decimal )

if __name__ == '__main__':
    unittest.main()
