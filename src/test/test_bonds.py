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

    def test_transfer_props(self):
        w = Water.get_standard()
        w.attach_properties()
        w.o.transfer_props()
        np.testing.assert_allclose( w.o.p.b, np.zeros(10), atol = 1e-7 )

    def test_transfer_props(self):
        w = Water.get_standard()
        w.attach_properties()
        w.populate_bonds()
        B = w.p.b.copy()

        w.h1.transfer_props()
        w.o.transfer_props()
        np.testing.assert_allclose( w.p.b, B, atol = 1e-7 )

    def test_transfer_props(self):
        w = Water.get_standard()
        w.attach_properties()
        w.populate_bonds()
        B = w.p.b.copy()

        transfer = { 'charge' : 1,
            'quadrupole' : 1,
            'dipole' : 1,
            'alpha' : 1,
            'beta' : 1 }

        w.transfer_props( [w.o, w.h1], transfer = transfer )
        np.testing.assert_allclose( w.p.b, B, atol = 1e-7 )






if __name__ == '__main__':
    unittest.main()
