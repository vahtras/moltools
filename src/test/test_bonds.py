import unittest, os
import numpy as np

from molecules import Cluster, Atom, Molecule
from use_generator import Generator

FILE_MOL = os.path.join( os.path.dirname(__file__), 'tip3p44_10qm.mol' )
FILE_PDB = os.path.join( os.path.dirname(__file__), 'tip3p0.pdb' )

class BondTestCase( unittest.TestCase ):

    def setUp(self):
        self.m = Molecule()

    def test_populate_bonds(self):

if __name__ == '__main__':
    unittest.main()
