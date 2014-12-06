import unittest, os
import numpy as np

from molecules import Cluster

FILE = os.path.join( os.path.dirname(__file__), 'tip3p44_10qm.mol' )

class WaterTest( unittest.TestCase ):
    def setUp(self):
        self.c = Cluster.get_water_cluster(
                FILE,
                in_AA = False,
                out_AA = False,
                N_waters = 10)

    def test_size(self):
        assert len( self.c) == 10

    def test_min_len(self):
        assert len( self.c.min_dist_coo() ) == 9

    def test_dists(self):
        print len( self.c )
        print self.c.min_dist_coo()[:len(self.c)]

    def test_two_waters(self):
        c = Cluster.get_water_cluster(
                FILE,
                in_AA = False,
                out_AA = False,
                N_waters = 2)
        assert len(c) == 2
        print c.min_dist_coo()
        assert c.min_dist_coo()[0] > 4.5

    def test_three_waters(self):
        c = Cluster.get_water_cluster(
                FILE,
                in_AA = False,
                out_AA = False,
                N_waters = 3)
        assert len(c) == 3
        assert c.min_dist_coo()[0] > 4.5

        dists = c.min_dist_coo()

        assert np.mean( [ dists[0], dists[1] ] ) == \
                (dists[0] + dists[1])/2


    def eq(self, a, b):
        np.testing.assert_almost_equal( a, b, decimal = 3)



if __name__ == '__main__':
    unittest.main()
