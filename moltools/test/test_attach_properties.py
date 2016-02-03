import unittest, os, warnings
import numpy as np

warnings.simplefilter('error')
from nose.plugins.attrib import attr
from moltools import Molecule, Property, Water
from moltools import utilz

WATER = """3

OW               22.96300  30.31700  45.73800
HW1              23.12300  30.72100  46.59000
HW2              22.51700  30.99500  45.23000"""

@attr(speed = 'fast' )
class AttachPropertyTest( unittest.TestCase ):

    def setUp(self):
        pass

    def test_attach_properties(self):
        w = Water.get_standard()
        t, r1, r2, r3 = utilz.center_and_xz( w.o.r, (w.h1.r - w.h2.r)/2.0 + w.h2.r, w.h1.r )
        w.o.pdb_name = 'OW'
        w.h1.pdb_name = 'HW1'
        w.h2.pdb_name = 'HW2'
        w.attach_properties( 
                model = 'TIP3P_PDB',
                method = 'B3LYP',
                basis ='ANO631',
                template_key = lambda x: x.pdb_name,
                force_template = True
                )
        np.testing.assert_allclose( w.p.q , 0.0, atol = 1e-4 )

    def test_attach_to_water(self):

        w = Molecule.from_xyz_string( WATER )
        p1, p2, p3 = w[0].r, (w[1].r - w[2].r)/2.0 + w[2].r, w[1].r 
        t, r1, r2, r3 = utilz.center_and_xz( p1, p2, p3 )
        w[0].pdb_name = 'OW'
        w[1].pdb_name = 'HW1'
        w[2].pdb_name = 'HW2'
        w.attach_properties( 
                model = 'TIP3P_PDB',
                method = 'B3LYP',
                basis ='ANO631',
                template_key = lambda x: x.pdb_name,
                force_template = True
                )
        np.testing.assert_allclose( w.p.q , 0.0, atol = 1e-4 )

    def test_centered(self):

        w = Molecule.from_xyz_string( WATER )
        p1, p2, p3 = w[0].r, (w[1].r - w[2].r)/2.0 + w[2].r, w[1].r 
        t, r1, r2, r3 = utilz.center_and_xz( p1, p2, p3 )
        w.attach_properties( 
                model = 'TIP3P',
                method = 'B3LYP',
                basis ='ANO631',
                loprop = 0, 
                template_key = lambda x: x.pdb_name,
                force_template = True
                )

        np.testing.assert_equal( w.is_Property, True )
        np.testing.assert_allclose( w.Property.q , 0.0, atol = 1e-4 )

    def test_centered_loprop(self):

        w = Water.get_standard()

        w.attach_properties()
        w.rotate( *np.random.random( (3,) ))
        w.t( np.random.uniform( 0, 5, (3,) ))
        p_ref = w.p.copy()

        w.attach_properties()
        np.testing.assert_allclose( w.p.d, p_ref['dipole'], atol = 1e-7 )
        np.testing.assert_allclose( w.p.b, p_ref['beta'], atol = 1e-7 )
        np.testing.assert_allclose( w.p.a, p_ref['alpha'], atol = 1e-7 )

    def test_centered_simple(self):

        w = Water.get_standard()

        w.attach_properties( loprop = 0 )
        w.rotate( *np.random.random( (3,) ))
        w.t( np.random.uniform( 0, 5, (3,) ))
        p_ref = w.p.copy()

        w.attach_properties( loprop = 0 )
        np.testing.assert_allclose( w.p.d, p_ref['dipole'], atol = 1e-7 )
        np.testing.assert_allclose( w.p.b, p_ref['beta'], atol = 1e-7 )
        np.testing.assert_allclose( w.p.a, p_ref['alpha'], atol = 1e-7 )



#
#        w.attach_properties2( 
#                model = 'TIP3P',
#                method = 'B3LYP',
#                basis ='ANO631',
#                loprop = 0, 
#                template_key = lambda x: x.pdb_name,
#                force_template = True
#                )
#
#        np.testing.assert_equal( w.is_Property, True )
#        np.testing.assert_allclose( w.Property.q , 0.0, atol = 1e-4 )
#
#






if __name__ == '__main__':
    unittest.main()
