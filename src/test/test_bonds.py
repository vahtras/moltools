import unittest, os
import numpy as np

from molecules import Cluster, Atom, Molecule, Water
from use_generator import Generator

FILE_MOL = os.path.join( os.path.dirname(__file__), 'tip3p44_10qm.mol' )
FILE_PDB = os.path.join( os.path.dirname(__file__), 'tip3p0.pdb' )

class BondTestCase( unittest.TestCase ):

    def test_populate_bonds(self):
        w = Water.get_standard()
        w.populate_bonds()
        assert len( w.o.bonds ) == 2
        assert len( w.h1.bonds ) == 1
        assert len( w.h2.bonds ) == 1


if __name__ == '__main__':
    unittest.main()
