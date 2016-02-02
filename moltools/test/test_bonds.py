import unittest, os, warnings
import numpy as np

from moltools import Cluster, Atom, Molecule, Water, Generator
warnings.simplefilter('error')
from nose.plugins.attrib import attr


FILE_MOL = os.path.join( os.path.dirname(__file__), 'tip3p44_10qm.mol' )
FILE_PDB = os.path.join( os.path.dirname(__file__), 'tip3p0.pdb' )

@attr(speed = 'fast' )
class BondTestCase( unittest.TestCase ):

    def test_populate_bonds(self):
        w = Water.get_standard()
        w.populate_bonds()
        assert len( w.o.bonds ) == 2
        assert len( w.h1.bonds ) == 1
        assert len( w.h2.bonds ) == 1

    #def test_transfer_props_1(self):
    #    w = Water.get_standard()
    #    w.attach_properties()
    #    w.populate_bonds()
    #    w.o.transfer_props( {'beta' : 1} )
    #    np.testing.assert_allclose( w.o.p.b, np.zeros(10), atol = 1e-7 )

    #def test_transfer_props_2(self):
    #    w = Water.get_standard()
    #    w.attach_properties()
    #    w.populate_bonds()
    #    B = w.p.b.copy()

    #    w.h1.transfer_props()
    #    w.o.transfer_props()
    #    np.testing.assert_allclose( w.p.b, B, atol = 1e-7 )

    #def test_transfer_props_3(self):
    #    w = Water.get_standard()
    #    w.attach_properties()
    #    w.populate_bonds()
    #    C = w.p.q.copy()
    #    #C = w.p.q.copy()
    #    #D = w.p.d.copy()
    #    #Q = w.p.Q.copy()
    #    A = w.p.a.copy()
    #    B = w.p.b.copy()

    #    transfer = { 'charge' : 1,
    #        'quadrupole' : 1,
    #        'dipole' : 1,
    #        'alpha' : 1,
    #        'beta' : 1 }

    #    w.transfer_props( [w.o, w.h1], transfer = transfer )
    #    #np.testing.assert_allclose( w.p.q, C, atol = 1e-7 )
    #    #np.testing.assert_allclose( w.p.d, D, atol = 1e-7 )
    #    #np.testing.assert_allclose( w.p.Q, Q, atol = 1e-7 )
    #    np.testing.assert_allclose( w.p.a, A, atol = 1e-7 )
    #    np.testing.assert_allclose( w.p.b, B, atol = 1e-7 )

    #def test_transfer_props_4(self):
    #    w = Water.get_standard()
    #    w.attach_properties()
    #    w.populate_bonds()

    #    q_h1 = w.h1.p.q
    #    d_h1 = w.h1.p.d.copy()
    #    Q_h1 = w.h1.p.Q.copy()
    #    a_h1 = w.h1.p.a.copy()
    #    b_h1 = w.h1.p.b.copy()
    #    
    #    q_o = w.o.p.q

    #    w.h1.transfer_props( {'charge' : 1 }  )
    #    np.testing.assert_allclose( w.h1.p.q, 0.0, atol = 1e-7 )
    #    np.testing.assert_allclose( w.h1.p.d, d_h1, atol = 1e-7 )
    #    np.testing.assert_allclose( w.h1.p.Q, Q_h1, atol = 1e-7 )
    #    np.testing.assert_allclose( w.h1.p.a, a_h1, atol = 1e-7 )
    #    np.testing.assert_allclose( w.h1.p.b, b_h1, atol = 1e-7 )

    #    np.testing.assert_allclose( w.o.p.q, (q_h1 + q_o), atol = 1e-7 )
#

    #def test_transfer_props_5(self):
    #    w = Water.get_standard()
    #    w.attach_properties()
    #    w.populate_bonds()

    #    q_h1 = w.h1.p.q
    #    d_h1 = w.h1.p.d.copy()
    #    Q_h1 = w.h1.p.Q.copy()
    #    a_h1 = w.h1.p.a.copy()
    #    b_h1 = w.h1.p.b.copy()
    #    w.h1.transfer_props( { 'charge':1, 'dipole' : 1 }  )
    #    np.testing.assert_allclose( w.h1.p.q, 0.0, atol = 1e-7 )
    #    np.testing.assert_allclose( w.h1.p.d, np.zeros(3), atol = 1e-7 )
    #    np.testing.assert_allclose( w.h1.p.Q, Q_h1, atol = 1e-7 )
    #    np.testing.assert_allclose( w.h1.p.a, a_h1, atol = 1e-7 )
    #    np.testing.assert_allclose( w.h1.p.b, b_h1, atol = 1e-7 )

    #def test_transfer_props_6(self):
    #    w = Water.get_standard()
    #    w.attach_properties()
    #    w.populate_bonds()

    #    q_h1 = w.h1.p.q
    #    d_h1 = w.h1.p.d.copy()
    #    Q_h1 = w.h1.p.Q.copy()
    #    a_h1 = w.h1.p.a.copy()
    #    b_h1 = w.h1.p.b.copy()
    #    w.h1.transfer_props( { 'charge':1, 'dipole' : 1, 'quadurpole':1 }  )
    #    np.testing.assert_allclose( w.h1.p.q, 0.0, atol = 1e-7 )
    #    np.testing.assert_allclose( w.h1.p.d, np.zeros(3), atol = 1e-7 )
    #    np.testing.assert_allclose( w.h1.p.Q, np.zeros(6), atol = 1e-7 )
    #    np.testing.assert_allclose( w.h1.p.a, a_h1, atol = 1e-7 )
    #    np.testing.assert_allclose( w.h1.p.b, b_h1, atol = 1e-7 )

    #def test_transfer_props_7(self):
    #    w = Water.get_standard()
    #    w.attach_properties()
    #    w.populate_bonds()

    #    q_h1 = w.h1.p.q
    #    d_h1 = w.h1.p.d.copy()
    #    Q_h1 = w.h1.p.Q.copy()
    #    a_h1 = w.h1.p.a.copy()
    #    b_h1 = w.h1.p.b.copy()
    #    w.h1.transfer_props( { 'charge':1,
    #        'alpha' : 1,
    #        'dipole' : 1, 'quadurpole':1 }  )
    #    np.testing.assert_allclose( w.h1.p.q, 0.0, atol = 1e-7 )
    #    np.testing.assert_allclose( w.h1.p.d, np.zeros(3), atol = 1e-7 )
    #    np.testing.assert_allclose( w.h1.p.Q, np.zeros(6), atol = 1e-7 )
    #    np.testing.assert_allclose( w.h1.p.a, np.zeros(6), atol = 1e-7 )
    #    np.testing.assert_allclose( w.h1.p.b, b_h1, atol = 1e-7 )

    #def test_transfer_props_8(self):
    #    w = Water.get_standard()
    #    w.attach_properties()
    #    w.populate_bonds()

    #    q_h1 = w.h1.p.q
    #    d_h1 = w.h1.p.d.copy()
    #    Q_h1 = w.h1.p.Q.copy()
    #    a_h1 = w.h1.p.a.copy()
    #    b_h1 = w.h1.p.b.copy()
    #    w.h1.transfer_props( { 'charge':1,
    #        'alpha' : 1, 'beta' : 1,
    #        'dipole' : 1, 'quadurpole':1 }  )
    #    np.testing.assert_allclose( w.h1.p.q, 0.0, atol = 1e-7 )
    #    np.testing.assert_allclose( w.h1.p.d, np.zeros(3), atol = 1e-7 )
    #    np.testing.assert_allclose( w.h1.p.Q, np.zeros(6), atol = 1e-7 )
    #    np.testing.assert_allclose( w.h1.p.a, np.zeros(6), atol = 1e-7 )
    #    np.testing.assert_allclose( w.h1.p.b, np.zeros(10), atol = 1e-7 )




if __name__ == '__main__':
    unittest.main()
