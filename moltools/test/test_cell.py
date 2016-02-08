import unittest, os, warnings, warnings
import numpy as np
from applequist.particles import PointDipoleList

warnings.simplefilter('error')
from nose.plugins.attrib import attr

from moltools import Cluster, Atom, Generator
from moltools import dstruct
from moltools.dstruct import Cell

FILE_XYZ = os.path.join( os.path.dirname(__file__), 'pna_waters.xyz' )
FILE_MOL = os.path.join( os.path.dirname(__file__), 'tip3p44_10qm.mol' )
FILE_PDB = os.path.join( os.path.dirname(__file__), 'tip3p0.pdb' )
POTSTRING = """AU
6 1 22 1
1      0.000000   0.000000   0.000000 -0.66229 0.00000 0.00000 0.34276 4.10574 0.00000 0.00000 4.79229 0.00000 4.01912 0.00000 0.00000 -3.33162 0.00000 0.00000 0.00000 0.00000 -0.32216 0.00000 0.79137
1      1.430429   0.000000   1.107157 0.33114 -0.16617 0.00000 -0.11629 1.53802 0.00000 1.19765 0.90661 0.00000 1.37138 -4.52137 0.00000 -5.08061 -1.35494 0.00000 -4.83365 0.00000 -0.46317 0.00000 -3.47921
1     -1.430429   0.000000   1.107157 0.33114 0.16617 0.00000 -0.11629 1.53802 0.00000 -1.19765 0.90661 0.00000 1.37138 4.52137 0.00000 -5.08061 1.35494 0.00000 4.83365 0.00000 -0.46317 0.00000 -3.47921
2     15.000000  15.000000  15.000000 -0.66229 0.00000 0.00000 0.34276 4.10574 0.00000 0.00000 4.79229 0.00000 4.01912 0.00000 0.00000 -3.33162 0.00000 0.00000 0.00000 0.00000 -0.32216 0.00000 0.79137
2     16.430429  15.000000  16.107157 0.33114 -0.16617 0.00000 -0.11629 1.53802 0.00000 1.19765 0.90661 0.00000 1.37138 -4.52137 0.00000 -5.08061 -1.35494 0.00000 -4.83365 0.00000 -0.46317 0.00000 -3.47921
2     13.569571  15.000000  16.107157 0.33114 0.16617 0.00000 -0.11629 1.53802 0.00000 -1.19765 0.90661 0.00000 1.37138 4.52137 0.00000 -5.08061 1.35494 0.00000 4.83365 0.00000 -0.46317 0.00000 -3.47921"""


@attr(speed = 'fast' )
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
        c = c.update()
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

    def test_from_PointDipoleList(self, ):
        _str = POTSTRING
        pdl = PointDipoleList.from_string( _str )
        cell = dstruct.Cell.from_PointDipoleList( pdl, co = 5 )

        assert isinstance( cell, dstruct.Cell )



if __name__ == '__main__':
    unittest.main()
