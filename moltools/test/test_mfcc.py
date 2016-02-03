import unittest, os, time
import numpy as np

import warnings
warnings.simplefilter('error')
from nose.plugins.attrib import attr
from moltools import utilz
from moltools import Residue, Atom, Property

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
        res1 = Residue( Atom(x=0),
               Atom(x=1),
               Atom(x=2),
               Atom(x=3),
               Atom(x=4),
               Atom(x=5),
               Atom(x=1, y = 1),
               Atom(x=3, y = 1),
               )

if __name__ == '__main__':
    unittest.main()
