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
        pass
        
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

    def test_update(self):
        c = Cell( my_cutoff = 3 )
        a = Atom( z = 5 )
        c.add( a )
        assert a in c[0][0][1]
        a.z = 0 
        c.update()
        assert a in c[0][0][0]

    def test_get_closest(self):
        cell = Cell.from_xyz( FILE_XYZ )
#ensure at1 exists
        for at in cell:
            at1 = at
        x, y, z = cell.get_index( at1 )
        ats = 0
        tmp = []
        for i in range( x-1, x+2 ):
            for j in range( y-1, y+2 ):
                for k in range( z-1, z+2 ):
                    try:
                        for at in cell[i][j][k]:
                            if at in tmp:
                                continue
                            tmp.append(at)
                    except IndexError:
                        pass
        assert len(tmp) -1 == len(cell.get_closest( at1 ))


if __name__ == '__main__':
    unittest.main()
