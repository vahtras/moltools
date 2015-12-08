import unittest, os
import numpy as np

import warnings
warnings.simplefilter('error')
from nose.plugins.attrib import attr
from molecules import Atom

@attr(speed = 'fast' )
class AtomTestCase( unittest.TestCase ):

    def test_translate(self):
        h = Atom( element = 'H' )
        h.t( np.array( [1, 2, 3,] ) )
        np.testing.assert_allclose( [h.x, h.y, h.z], np.array([1, 2, 3,]) )

if __name__ == '__main__':
    unittest.main()
