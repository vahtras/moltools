import unittest, os
import numpy as np

import warnings
warnings.simplefilter('error')
from nose.plugins.attrib import attr
import utilz
from pdbreader import NewResidue, NewAtom
from property import Property
import time

class TimeException( Exception ):
    def __init__(self, val):
        self.val = val
        
def f( x, *args, **kwargs ):
    def wrapped(  *args, **kwargs ):
        d1 = time.time()
        val = x( *args, **kwargs )
        delta = time.time() - d1
        if delta > 0.01:
            raise TimeException( delta )
        return val
    return wrapped

@attr(speed = 'fast' )
class MfccTestCase( unittest.TestCase ):

    def test_mfcc(self):

#For res 1
        res1 = NewResidue( NewAtom(x=0),
               NewAtom(x=1),
               NewAtom(x=2),
               NewAtom(x=3),
               NewAtom(x=4),
               NewAtom(x=5),
               NewAtom(x=1, y = 1),
               NewAtom(x=3, y = 1),
               )

    @attr(speed = 'slow' )
    def test_props_from_targz(self):
        assert 0
        pass

if __name__ == '__main__':
    unittest.main()
