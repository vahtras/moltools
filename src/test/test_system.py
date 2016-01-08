import unittest, os
import numpy as np

import warnings
warnings.simplefilter('error')
from nose.plugins.attrib import attr

from moltools import Cluster, Atom, Water, System, Chain, Generator


@attr(speed = 'fast' )
class SystemText( unittest.TestCase ):

    def test_three_waters(self):
        s = System()

        ch1, ch2 = Chain(), Chain()
        ch1.chain_id = 'A'
        ch2.chain_id = 'B'

        w1, w2 = Water.get_standard( AA = True ), Water.get_standard( AA = True )

        ch1.add( w1 )
        ch2.add( w2 )

#After rotating and translating both waters oxygen is 2 AA apart
        w1.rotate( 0,  np.pi/2, 0 )
        w2.rotate( 0, -np.pi/2, 0 )

        w1.t( -1, 0, 0 )
        w2.t(  1, 0, 0 )

        s.add( ch1 )
        s.add( ch2 )

        assert len( s.min_dist_atoms_separate_res_chain( 2.1 ) ) == 2

        ch2.chain_id = 'A'

        assert len( s.min_dist_atoms_separate_res_chain( 2.1 ) ) == 0


if __name__ == '__main__':
    unittest.main()
