

__all__ = [ 'Rotator' ]

import numpy as np
import utilz 


class Rotator(object):
    """
**Container class for rotational operations on points, vectors, and tensors.**
"""

    def __init__(self):
        pass

    class RotatorError( Exception ):
        def __init__(self):
            pass

    @staticmethod
    def b_hrs( b):
        if b.shape == (10,):
            b = Rotator.ut_3_square(b)
        elif b.shape != (3,3,3,):
            print "supplied wrong beta"
            raise RotatorError

        zzz = Rotator.rot_avg( b )
        xzz = Rotator.rot_avg( b, car1=0 )
        return np.sqrt( zzz + xzz )

    @staticmethod
    def dr( b ):
        if b.shape == (10,):
            b = Rotator.ut_3_square(b)
        elif b.shape != (3,3,3,):
            print "supplied wrong beta"
            raise SystemExit

        zzz = Rotator.rot_avg( b )
        xzz = Rotator.rot_avg( b, car1=0 )
        return zzz / xzz

    @staticmethod
    def rot_avg( beta, car1 = 2, car2 = 2, car3 = 2):
        try:
            import h5py
        except ImportError:
            raise SystemExit('Install h5py in order to calculate euler rotational averages')
        """
        Requires euler.h5 binary file containing rotational angle products
        Define it as in current script directory + euler.h5
        """
        b_new = np.zeros( (3,3,3,) )
        """given beta in molecular frame, convert to exp. reference"""
        vec = h5py.File( os.path.join(os.path.dirname( os.path.realpath( __file__ )), 'euler.h5' ), 'r')['data'].value
        for X in range(3):
            if X != car1:
                continue
            for Y in range(3):
                if Y != car2:
                    continue
                for Z in range(3):
                    if Z != car3:
                        continue
                    for x1 in range(3):
                        for y1 in range(3):
                            for z1 in range(3):
                                for x2 in range(3):
                                    for y2 in range(3):
                                        for z2 in range(3):
                                            b_new[X,Y,Z] += vec[X,Y,Z,x1,y1,z1,x2,y2,z2] * beta[x1,y1,z1] * beta[x2,y2,z2]
        return b_new[ car1, car2, car3 ]

    @staticmethod
    def transform_1( qm_dipole, t1, t2, t3 ):
        """
Rotate vector around z-axis clockwise by :math:`\\rho_{1}`, around the y-axis counter-clockwise by :math:`\\rho_2`, and finally clockwise around the z-axis by :math:`\\rho_3`.

.. code:: python

    >>> import numpy as np
    >>> d = np.array( [ 1, 0, 0] )
    >>> print Rotator.transform_1( d, 0, numpy.pi/2, 0 )
    [ 0.0, 0.0, 1.0 ]
"""
        d_new1 = np.zeros([3]) #will be returned
        d_new2 = np.zeros([3]) #will be returned
        d_new3 = np.zeros([3]) #will be returned

        rz  = Rotator.get_Rz( t1 )
        ryi = Rotator.get_Ry_inv( t2 )
        rz2 = Rotator.get_Rz( t3 )

        for i in range(3):
            for x in range(3):
                d_new1[i] += rz[i][x] * qm_dipole[x]
        for i in range(3):
            for x in range(3):
                d_new2[i] += ryi[i][x] * d_new1[x]
        for i in range(3):
            for x in range(3):
                d_new3[i] += rz2[i][x] * d_new2[x]
        return d_new3
    @staticmethod
    def transform_2( qm_alpha, t1, t2 , t3 ):
        a_new1 = np.zeros([3,3]) #will be calculated
        a_new2 = np.zeros([3,3]) #will be calculated
        a_new3 = np.zeros([3,3]) #will be calculated

        rz  = Rotator.get_Rz( t1 )
        ryi = Rotator.get_Ry_inv( t2 )
        rz2 = Rotator.get_Rz( t3 )

        for i in range(3):
            for j in range(3):
                for x in range(3):
                    for y in range(3):
                        a_new1[i][j] += rz[i][x] * rz[j][y] * qm_alpha[x][y]

        for i in range(3):
            for j in range(3):
                for x in range(3):
                    for y in range(3):
                        a_new2[i][j] += ryi[i][x] * ryi[j][y] * a_new1[x][y]

        for i in range(3):
            for j in range(3):
                for x in range(3):
                    for y in range(3):
                        a_new3[i][j] += rz2[i][x] * rz2[j][y] * a_new2[x][y]

        return a_new3

    @staticmethod
    def inv_3( beta, t1, t2, t3):
        """Will inversely rotate tensor """
        assert beta.shape == (3,3,3)
        r1 = Rotator.get_Rz_inv( t1 )
        r2 = Rotator.get_Ry( t2 )
        r3 = Rotator.get_Rz_inv( t3 )
        return reduce(lambda a,x: np.einsum('ia,jb,kc,abc', x, x, x, a), [r1,r2,r3], beta )

    @staticmethod
    def transform_3( qm_beta, t1, t2, t3 ):
        b_new1 = np.zeros([3,3,3]) #will be calculated
        b_new2 = np.zeros([3,3,3]) #will be calculated
        b_new3 = np.zeros([3,3,3]) #will be calculated

        rz =  Rotator.get_Rz( t1 )
        ryi = Rotator.get_Ry_inv( t2 )
        rz2 = Rotator.get_Rz( t3 )

        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for x in range(3):
                        for y in range(3):
                            for z in range(3):
                                b_new1[i][j][k] += rz[i][x] * rz[j][y] * rz[k][z] * qm_beta[x][y][z]
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for x in range(3):
                        for y in range(3):
                            for z in range(3):
                                b_new2[i][j][k] += ryi[i][x] * ryi[j][y] * ryi[k][z] * b_new1[x][y][z]
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for x in range(3):
                        for y in range(3):
                            for z in range(3):
                                b_new3[i][j][k] += rz2[i][x] * rz2[j][y] * rz2[k][z] * b_new2[x][y][z]
        return b_new3

    @staticmethod
    def square_2_ut(alpha):
        assert alpha.ndim == 2
        tmp_a = np.zeros( 6 )
        for index, (i, j ) in enumerate( utilz.upper_triangular(2) ):
            tmp_a[ index ] = (alpha[i, j] + alpha[ j, i]) / 2
        return tmp_a

    @staticmethod
    def get_Rz( theta ):
        vec = np.array(    [[ np.cos(theta),-np.sin(theta), 0],
                            [ np.sin(theta), np.cos(theta), 0],
                            [ 0,    0,  1]])
        return vec
    @staticmethod
    def get_Rz_inv( theta ):
        vec = np.array(     [[ np.cos(theta), np.sin(theta), 0],
                            [ -np.sin(theta), np.cos(theta), 0],
                            [ 0,             0,            1]])
        return vec
    @staticmethod
    def get_Ry( theta ):
        vec = np.array(    [[ np.cos(theta),0, np.sin(theta)],
                            [ 0,    1,  0],
                            [ -np.sin(theta), 0, np.cos(theta)]])
        return vec
    @staticmethod
    def get_Ry_inv( theta ):
        vec = np.array(    [[ np.cos(theta),0, -np.sin(theta)],
                            [ 0,    1,  0],
                            [ np.sin(theta), 0, np.cos(theta)]])
        return vec

    @staticmethod
    def tensor_to_ut( beta ):
# naive solution, transforms matrix B[ (x,y,z) ][ (xx, xy, xz, yy, yz, zz) ] into array
# Symmtrized UT array    B[ (xxx, xxy, xxz, xyy, xyz, xzz, yyy, yyz, yzz, zzz) ]
        new = np.zeros( (10) )
        new[ 0 ] = beta[0,0]
        new[ 1 ] = (beta[0,1] + beta[1,0] ) /2
        new[ 2 ] = (beta[0,2] + beta[2,0] ) /2
        new[ 3 ] = (beta[0,3] + beta[1,1] ) /2
        new[ 4 ] = (beta[0,4] + beta[1,2] + beta[2,1] ) /3
        new[ 5 ] = (beta[0,5] + beta[2,2] ) /2
        new[ 6 ] = beta[1,3]
        new[ 7 ] = (beta[1,4] + beta[2,3] ) /2
        new[ 8 ] = (beta[1,5] + beta[2,4] ) /2
        new[ 9 ] = beta[2,5]
        return new
    @staticmethod
    def square_3_ut(beta):
        assert beta.ndim == 3
        tmp_b = np.zeros( 10 )
        for index, (i, j, k ) in enumerate( utilz.upper_triangular(3) ):
            tmp_b[ index ] = ( \
                    beta[i, j, k] + beta[i, k, j] + \
                    beta[j, i, k] + beta[j, k, i] + \
                    beta[k, i, j] + beta[k, j, i] )/ 6
        return tmp_b

    @staticmethod
    def ut2s( vec ):
        if len( vec ) == 6:
            return Rotator.ut_2_square( vec )
        elif len( vec ) == 10:
            return Rotator.ut_3_square( vec )

    @staticmethod
    def s2ut( vec ):
        if vec.shape == (3,3,):
            return Rotator.square_2_ut( vec )
        elif vec.shape == (3,3,3,):
            return Rotator.square_3_ut( vec )


    @staticmethod
    def ut_2_square( alpha):
        assert len(alpha) == 6
        tmp_a = np.zeros( (3,3, ))
        for index, val in enumerate( utilz.upper_triangular(2) ) :
            tmp_a[ val[0], val[1] ] = alpha[ index ]
            tmp_a[ val[1], val[0] ] = alpha[ index ]
        return tmp_a

    @staticmethod
    def ut_3_square( beta ):
        assert len(beta) == 10
        tmp_b = np.zeros( (3,3,3, ))
        for index, (i, j, k ) in enumerate( utilz.upper_triangular(3) ) :
            tmp_b[ i, j ,k] = beta[ index ]
            tmp_b[ i, k ,j] = beta[ index] 
            tmp_b[ j, i, k] = beta [ index ]
            tmp_b[ j, k, i] = beta [ index ]
            tmp_b[ k, i, j] = beta [ index ]
            tmp_b[ k, j, i] = beta [ index ]
        return tmp_b


