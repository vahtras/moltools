import unittest, os
import numpy as np
from molecules import Cluster, Atom
from use_generator import Generator

FILE_XYZ = os.path.join( os.path.dirname(__file__), 'pna_waters.xyz' )
FILE_MOL = os.path.join( os.path.dirname(__file__), 'tip3p44_10qm.mol' )
FILE_PDB = os.path.join( os.path.dirname(__file__), 'tip3p0.pdb' )

from dstruct import Cell
class CellTest( unittest.TestCase ):
    def setUp(self):
        self.c = Cluster.get_water_cluster(
                FILE_MOL,
                in_AA = False,
                out_AA = False,
                N_waters = 10)
        
    def test_init(self):
        c = Cell( my_min = map(float, [0, 0, 0]),
                    my_max = map(float, [1, 1, 1] ),
                    my_cutoff = 0.4)
        assert len(c) == 3

        c = Cell( my_min = map(float, [-10, 0, 0]),
                    my_max = map(float, [0, 1, 1] ),
                    my_cutoff = 12)
        assert len(c) == 1

        c = Cell( my_min = map(float, [-5, 0, 0]),
                    my_max = map(float, [10, 1, 1] ),
                    my_cutoff = 4.9)
        assert len(c) == 4

    def test_add(self):
        c = Cell( my_cutoff = 2.9 )
        a1 = Atom( element = 'H', x = 3 )
        a2 = Atom( element = 'H', x = 3, y = 3 )
        a3 = Atom( element = 'H', x = 3, y = 3, z= 3 )
        c.add( a1 )
        c.add( a2 )
        c.add( a3 )
        assert a1 in c[1][0][0]
        assert a2 in c[1][1][0]
        assert a3 in c[1][1][1]

if __name__ == '__main__':
    unittest.main()
