
import unittest, os
import numpy as np

from nose.plugins.attrib import attr
import utilz
from pdbreader import NewResidue, NewAtom
from property import Property

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

if __name__ == '__main__':
    unittest.main()
