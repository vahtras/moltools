import unittest, os
import numpy as np

from molecules import Cluster, Property, Water
from template import Template

FILE = os.path.join( os.path.dirname(__file__), 'tip3p44_10qm.mol' )

class WaterTest( unittest.TestCase ):
    def setUp(self):
        self.c = Cluster.get_water_cluster(
                FILE,
                in_AA = False,
                out_AA = False,
                N_waters = 10)
        

    def test_size(self):
        c = Cluster.get_water_cluster(
                FILE,
                in_AA = False,
                out_AA = False,
                N_waters = 10)

# Read in distributed properties, transform to atomic sites from waters euler angles
        for wat in c:
            t1, t2, t3  = wat.get_euler()
            kwargs_dict = Template().get( *("TIP3P", "HF", "ANOPVDZ",
                True , "0.0" ))
            for at in wat:
                Property.add_prop_from_template( at, kwargs_dict )
                Property.transform_ut_properties( at.Property, t1, t2 ,t3)

# Read in the properties for the oxygen atom, the projected dipole vector should
# be the same
        for wat in c:
            d = np.zeros( (3,) )
            a = np.zeros( (3,3,) )
            b = np.zeros( (3,3,3,) )

            kwargs_dict = Template().get( *("TIP3P", "HF", "ANOPVDZ",
                False , "0.0" ))
            d_ref = kwargs_dict[('O1','dipole')]
            b_ref = Water.ut_3_square( kwargs_dict[('O1','beta')] )
            for at in wat:
                d += at.Property["dipole"] + at.Property["charge"] * at.r
                a += Water.ut_2_square( at.Property["alpha"] )
                b += Water.ut_3_square( at.Property["beta"])

            self.eq( np.dot( np.einsum('iij->j',b), d)/np.linalg.norm(d), 
                    np.dot( np.einsum('iij->j',b_ref),d_ref)/np.linalg.norm(d_ref) )




    def eq(self, a, b,):
        np.testing.assert_almost_equal( a, b, decimal = 2 )



if __name__ == '__main__':
    unittest.main()
