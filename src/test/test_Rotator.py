import unittest
import numpy as np

from molecules import Rotator

class RotatorTest( unittest.TestCase ):

    def setUp(self):
        pass

    def test_ut2s_alpha(self):
        v = np.full( (6,), 5 )
        assert Rotator.ut2s( v ).shape == (3 ,3 )
        np.testing.assert_allclose( Rotator.ut2s( v ), np.full( (3,3,), 5 ) )

    def test_ut2s_beta(self):
        v = np.full( (10,), -3 )
        assert Rotator.ut2s( v ).shape == (3 ,3, 3 )
        np.testing.assert_allclose( Rotator.ut2s( v ), np.full( (3,3,3,), -3 ) )

    def test_s2ut_alpha(self):
        v = np.full( (3,3), 83 )
        assert Rotator.s2ut( v ).shape == (6,)
        np.testing.assert_allclose( Rotator.s2ut( v ), np.full( (6,), 83 ) )

    def test_s2ut_beta(self):
        v = np.full( (3,3,3), -5 )
        assert Rotator.s2ut( v ).shape == (10,)
        np.testing.assert_allclose( Rotator.s2ut( v ), np.full( (10,), -5 ) )

if __name__ == '__main__':
    unittest.main()
