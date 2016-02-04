import unittest
import numpy as np

from ..utilz import o_filter
import warnings
warnings.simplefilter('error')
from nose.plugins.attrib import attr

@attr(speed = 'fast' )
class OFilterrTest( unittest.TestCase ):

    def setUp(self):
        pass
    def test_o_filter(self):
        files = ["hfqua_9.23-0.00-0.00-0.00-0.00-0.00.out",
                 "hfqua_3.00-0.00-0.00-0.00-0.00-0.00.out",
                 "hfqua_3.00-0.00-0.00-0.00-0.00-0.00.out",
                 "hfqua_3.00-1.57-0.00-0.00-0.00-0.00.out",
                 "hfqua_3.00-1.57-0.00-0.00-0.00-0.00.out",
                 "hfqua_3.78-0.00-0.55-0.00-0.00-0.00.out",
                 "hfqua_4.56-0.00-0.55-0.00-0.00-0.00.out",
                 "hfqua_5.33-0.00-0.55-0.00-0.00-0.00.out",
                 "hfqua_6.11-0.00-0.55-0.33-0.00-0.00.out",
                 "hfqua_6.89-0.00-0.00-0.33-0.00-0.00.out",
                 "hfqua_7.67-0.00-0.00-0.33-0.00-0.00.out",
                 "hfqua_8.44-0.00-0.00-0.00-1.00-0.00.out",
                 "hfqua_9.22-0.00-0.00-0.00-0.00-1.00.out",
                ]
#Default should take first 3 files only
        out = o_filter( files, vary = 'r' )
        assert len(out) == 3

        out = o_filter( files, vary = 'rho3' )
        assert len(out) == 0

        out = o_filter( files, r = 9.22, vary = 'rho3' )
        print out
        assert len(out) == 1

if __name__ == '__main__':
    unittest.main()
