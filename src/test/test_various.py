import unittest, os
import numpy as np

from moltools import Rotator, Property, Generator, Template, Cluster, Water

import warnings
warnings.simplefilter('error')
from nose.plugins.attrib import attr

FILE = os.path.join( os.path.dirname(__file__), 'tip3p44_10qm.mol' )

FILE_STR = """ATOMBASIS
QM: WAT1 WAT173 WAT199 WAT245 WAT47 WAT227 WAT202 WAT306 WAT217 WAT65
MM: 
Atomtypes=2 Charge=0 Nosymm 
Charge=8.0 Atoms=10 Basis=ano-1 3 2 1
O      31.18048  24.32078  29.08289
O      33.78830  22.63892  25.30343
O      28.78053  28.64825  28.57266
O      29.32855  25.34123  34.03397
O      31.25607  19.65315  31.86078
O      36.58510  24.18849  30.76474
O      26.81521  23.16804  24.62313
O      29.89547  27.24985  23.07356
O      26.32388  19.53977  30.40569
O      35.31898  28.91281  32.82454
Charge=1.0 Atoms=20 Basis=ano-1 2 1
H      30.89702  26.09712  28.87502
H      30.29231  23.88614  30.59467
H      32.99462  21.29721  24.37747
H      32.59778  23.05466  26.60734
H      27.17426  28.02464  29.10178
H      28.42148  30.17893  27.66559
H      27.55221  24.96328  34.10956
H      29.38524  27.13647  33.88279
H      31.23717  18.08468  32.74895
H      32.27652  20.73030  32.90013
H      36.18826  25.88925  31.29386
H      35.26229  23.75386  29.61201
H      28.55376  22.65782  24.60423
H      25.88925  21.65626  24.94438
H      28.11912  27.02308  23.28143
H      30.53797  27.30654  24.75541
H      25.47351  20.93817  31.18048
H      28.08133  19.86102  30.67026
H      36.83076  29.64980  33.52374
H      34.33632  28.44038  34.26073
"""

@attr(speed = 'fast' )
class WaterTest( unittest.TestCase ):

    def setUp(self):
        self.c = Cluster.get_water_cluster(
                FILE,
                in_AA = False,
                out_AA = False,
                N_waters = 10)

    def test_beta_rotation(self):
        """reference beta is for a tip3p rotated 1.57 radians round z axis"""
        bxxz = -2.2760
        byyz = -14.409
        bzzz = -5.851
        #w = Generator().get_mol( [0,0,0] )
        v = np.zeros( (3,3,3))

# This is the vectors when alligned in x-z plane, switched xxz and yyz components
        v[0,0,2] = byyz
        v[1,1,2] = bxxz
        v[2,2,2] = bzzz
        vnew = Rotator.transform_3( v, np.pi/2, 0,0  )
        self.eq( vnew[0,0,2], bxxz )
        self.eq( vnew[1,1,2], byyz )
        self.eq( vnew[2,2,2], bzzz )

    def test_size(self):
        c = Cluster.get_water_cluster_from_string(
                FILE_STR,
                in_AA = False,
                out_AA = False,
                N_waters = 10)

# Read in distributed properties, transform to atomic sites from waters euler angles
        loprop = True
        for wat in c:
            t1, t2, t3  = wat.get_euler( key = lambda x: (x.o.r, (x.h1.r-x.h2.r)/2 + x.h2.r, x.h1.r ))
            kwargs_dict = Template().get( *("TIP3P", "HF", "ANOPVDZ",
                loprop, "0.0" ))
            for at in wat:
                at.p = Property.from_template( at.name, kwargs_dict )
                at.p = at.p.rotate( t1, t2 ,t3)

# Read in the properties for the oxygen atom, the projected dipole vector should
# be the same
        for wat in c:
            q = np.zeros( (1, ) )
            d = np.zeros( (3, ) )
            a = np.zeros( (3,3, ) )
            b = np.zeros( (3,3,3, ) )

            for at in wat:
                q += at.Property["charge"]
                d += at.Property["dipole"] + at.Property["charge"] * (at.r-wat.coc)
                a += Rotator.ut_2_square( at.Property["alpha"] )
                b += Rotator.ut_3_square( at.Property["beta"] )

#False loprop means we just label the point by 'X'
            kwargs_dict = Template().get( *("TIP3P", "HF", "ANOPVDZ",
                False , "0.0" ))
            d_ref = kwargs_dict[('X','dipole')]
            b_ref = Rotator.ut_3_square( kwargs_dict[('X','beta')] )

            self.eq( q, 0.0 )
            self.eq( np.dot( np.einsum('iij->j',b), d)/np.linalg.norm(d), 
                    np.dot( np.einsum('iij->j',b_ref),d_ref)/np.linalg.norm(d_ref),
                    decimal = 4)

    def test_tensor_to_ut(self):
        t = np.zeros( (3, 6) )
        self.eq( t, np.zeros( ( 3, 6, )) )
#xxx component
        t[0,0] = -3.0
#yyy component
        t[1,3] = 0.25
#zzz component
        t[2,5] = 5.0

#xxy component
        t[0,1] = 3.5
        t[1,0] = 4.5

#xxz component
        t[0,2] = 2.0
        t[2,0] = 3.0

#xyy component
        t[0,3] = -2.0

#xyz component
        t[0,4] = 1.0
        t[1,2] = 3.0
        t[2,1] = 8.0
        t[1,1] = 2.0
#xzz component
        t[0,5] = 2.0
        t[2,2] = -1.0
#yyz component
        t[1,4] = 3.5
        t[2,3] = 7.5
#yzz component
        t[1,5] = -4.0
        t[2,4] = -6.0

# AFTER TRANSFORMAATION

        square = Rotator.tensor_to_ut( t )

#xxx, one permutations, -3.0
        self.eq( square[0], -3.0 )
#xxy, two permutations, 3.5 and 4.5 should be 3.0
        self.eq( square[1], 4.0 )
#xxz, two permutations, 2.0 and 3.0 should be 2.5
        self.eq( square[2], 2.5 )
#xyy, two permutations, -2.0 and 2.0 should be 0
        self.eq( square[3], 0.0 )
#xyz, three permutations, 1.0 and 3.0  and 8.0 should be 4.0
        self.eq( square[4], 4.0 )
#xzz, two permutations, 2.0 and -1.0 should be 0.5
        self.eq( square[5], 0.5 )
#yyy, one permutations, 0.25
        self.eq( square[6], 0.25 )
#yyz, two permutations, 3.5 and 7.5, should be 5.5
        self.eq( square[7], 5.5 )
#yzz, two permutations, -4.0 and -6.0, should be -5.0
        self.eq( square[8], -5.0 )
#zzz, one permutations, 5.0
        self.eq( square[9], 5.0 )


    def test_property_add(self):
        p1 = Property()
        p1["charge"] = 0.3
        p2 = Property()
        p2["charge"] = 0.7
        assert (p1 + p2)["charge"] == 1.0

    def test_property_sub(self):
        p1 = Property()
        p1["charge"] = 0.3
        p2 = Property()
        p2["charge"] = 0.7
        self.eq( (p2 - p1)["charge"] , 0.4, decimal = 7 )

    def test_reflection(self):
        """docstring for test_reflection"""
        w = Water.get_standard()
        w.attach_properties( force_template = True)

        w.rotate( 0, np.pi/3.0, 0 )

        b_ref = w.p.b_proj.copy()
        d_ref = w.p.d.copy()
        w.reflect( plane = 'zy' )

        np.testing.assert_allclose( w.p.q, 0.0, atol = 1e-7 )
        np.testing.assert_allclose( w.p.b_proj, b_ref, atol = 1e-7 )
        np.testing.assert_allclose( abs(w.p.d), abs( d_ref ), atol = 1e-5 )

    def test_reflection2(self):
        w = Water.get_standard()
        w.attach_properties( force_template = True )

        q = w.p.q
        d = w.p.d.copy()
        Q = w.p.Q.copy()
        a = w.p.a.copy()
        b = w.p.b.copy()
        w.rotate( 0, np.pi/2, 0, )
        w.reflect( 'zy' )
        w.rotate( 0, np.pi/2, 0, )
        np.testing.assert_allclose( q, w.p.q, atol = 1e-7 )
        np.testing.assert_allclose( d, w.p.d, atol = 1e-7 )
        np.testing.assert_allclose( Q, w.p.Q, atol = 1e-7 )
        np.testing.assert_allclose( a, w.p.a, atol = 1e-7 )
        np.testing.assert_allclose( b, w.p.b, atol = 1e-7 )

    def eq(self, a, b, decimal = 7):
        np.testing.assert_almost_equal( a, b, decimal = decimal )

if __name__ == '__main__':
    unittest.main()
