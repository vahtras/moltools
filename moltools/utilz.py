#!/usr/bin/env python

__all__ = [ 'a0', 'o_filter', 'unique', 'splitter',
        's2ut', 'ut2s', 'nm_to_au' ]

"""

This will replace read_dal.py as general functions module with no dependencies
on other modules or files in this repository
"""

import os, sys, re, argparse, ctypes, multiprocessing, functools, warnings, itertools

import numpy as np
import matplotlib as mpl

#from .molecules import Atom

a0 = 0.52917721092
lab = [ "X", "Y", "Z"]
charge_dic = {"H1": 1.0 ,"H2":1.0 , "C1":6.0, "C7":6.0, "H3":1.0,
            "H4":1.0, "H6": 1.0, "H8":1.0, 
            "H9":1.0, "H10": 1.0, "H12":1.0, 
            "O5":8.0, "O11": 8.0,
            "H": 1.0, "C": 6.0, "N": 7.0, "O": 8.0, "S": 16.0}
mass_dict = {"H": 1.008,  "C": 6.0, "N": 7.0, "O": 15.999, "S": 16.0}
freq_dict = {"0.0": "static","0.0238927": "1907_nm", "0.0428227" : "1064_nm",
        "0.0773571" : "589_nm" }
allowed_elements = ( 'H', 'O' )
def nm_to_au( val ):
    conv = 45.563352491687866
    if val == 'inf':
        return "%.7f" % 0
    return "%.7f" %( conv / float(val) )

def au_to_nm( val ):
    conv = 45.563352491687866
    if np.allclose( 0.0, float(val), atol=1e-7 ):
        return 'inf'
    return "%d" % int( "%.0f" % ( conv/float(val) ))

def angle( p1, p2 ,p3 ):
    """Return the angle between 3 points"""
    v1 = p1 - p2
    v2 = p3 - p2
    return np.arccos( np.dot(v1,v2) / (np.linalg.norm(v1)* np.linalg.norm(v2)))


def make_para( shape = ( 0,) ):                                                      
    """Returns a np array that can be used for writing between different processes"""
    arr = np.zeros( shape )                                                          
    sb = multiprocessing.Array( ctypes.c_double, reduce(lambda a,b:a*b, arr.shape) )              
    sb = np.ctypeslib.as_array( sb.get_obj())                                        
    sb[:] = None                                                                     
    return sb.reshape( shape )  


def convex_hull_volume( pts):
    """Calculates the minimum convex set from points, in 
    3D ths corresponds to the minimum volume occupied by the points"""
    from scipy.spatial import Delaunay, ConvexHull
    def tetrahedron_volume(a, b, c, d):
        return np.abs(np.einsum('ij,ij->i', a-d, np.cross(b-d, c-d))) / 6
    ch = ConvexHull(pts)
    dt = Delaunay(pts[ch.vertices])
    tets = dt.points[dt.simplices]
    return np.sum(tetrahedron_volume(tets[:, 0], tets[:, 1],
        tets[:, 2], tets[:, 3]))

def rotate_point_by_two_points(p, p1, p2, theta = 0.0 ):
    """Rotate the point p around the line with point at p1 or p2 with 
    direction vector p2-p1"""
    origin = p1.copy()
    p, p1, p2 = map(lambda x:x-origin, [p, p1, p2] )
    r3, r2, r1 = get_euler( p2-p1, p-p1)

    p = np.einsum('ab,bc,cd,d', Rz_inv(r3), Ry(r2), Rz_inv(r1), p )
    p = np.einsum('ab,b', Rz(theta), p )
    p = np.einsum('ab,bc,cd,d', Rz(r1), Ry_inv(r2), Rz(r3), p )
    p += origin
    return p

def rotate_point_by_two_points2(p, p1, p2, theta = 0.0 ):
    """Rotate the point p around the line with point at p1 or p2 with 
    direction vector p2-p1"""
    return np.dot( R( p2-p1, theta ), p )



def rotate_point_around_cross(p1, p2, p3, theta = 0.0 ):
    """Rotate the point p1 clockwise 
    around the vector formed by crossing p1-p2 and p3-p2"""
    trans, r3, r2, r1 = center_and_xz( p2, p3, p1 )

    p1, p2, p3 = map(lambda x: x + trans, [p1, p2, p3] )

    p1, p2, p3 = map( lambda x: np.einsum('ab,bc,cd,d', Rz_inv(r3), Ry(r2), Rz_inv(r1), x ), [p1, p2, p3] )
    p1 = np.einsum('ab,b', Ry_inv(theta), p1 )
    p1, p2, p3 = map( lambda x: np.einsum('ab,bc,cd,d', Rz(r1), Ry_inv(r2), Rz(r3), x ), [p1, p2, p3] )

    p1, p2, p3 = map(lambda x: x - trans, [p1, p2, p3] )
    return p1

def reflect_point_by_three_points( p, p1, p2, p3 ):
    """ will reflect point p by 
    plane formed by points p1, p2 and p3"""
    v1 , v2 = p2 - p1, p3 - p1
    n = np.cross( v1, v2 )
    n = n / np.linalg.norm( n )
    r = np.dot( p, n ) * n / np.linalg.norm( n )
    r_orig = p - r
    return -r + r_orig

def reflect( plane = 'xz'):
    """Reflect by plane in cartesian 3D space"""
    vec = np.identity( 3 )
    opt_dict = { 'xy' : (2,2), 'xz' : (1, 1 ), 'yz' :(0,0,),
                 'yx' : (2,2), 'zx' : (1, 1 ), 'zy' :(0,0,) }
    if plane not in opt_dict:
        print 'Wrong plane given'
        raise SystemExit
    vec[ opt_dict[plane] ] = -1
    return vec
    
def Rz( theta ):
    vec = np.array(    [[ np.cos(theta),-np.sin(theta), 0],
                        [ np.sin(theta), np.cos(theta), 0],
                        [ 0,    0,  1]])
    return vec
def Rz_inv( theta ):
    vec = np.array(     [[ np.cos(theta), np.sin(theta), 0],
                        [ -np.sin(theta), np.cos(theta), 0],
                        [ 0,             0,            1]])
    return vec

def Ry( theta ):
    vec = np.array(    [[ np.cos(theta),0, np.sin(theta)],
                        [ 0,    1,  0],
                        [ -np.sin(theta), 0, np.cos(theta)]])
    return vec

def Ry_inv( theta ):
    vec = np.array(    [[ np.cos(theta),0, -np.sin(theta)],
                        [ 0,    1,  0],
                        [ np.sin(theta), 0, np.cos(theta)]])
    return vec

def Rx( theta ):
    vec = np.array(    [[ 1, 0.0, 0.0 ],
                        [ 0.0, np.cos(theta), -np.sin(theta)],
                        [ 0.0, np.sin(theta), np.cos(theta)]] )
    return vec

def Rx( theta ):
    vec = np.array(    [[ 1, 0.0, 0.0 ],
                        [ 0.0, np.cos(theta), np.sin(theta)],
                        [ 0.0, -np.sin(theta), np.cos(theta)]] )
    return vec




def R( axis, theta ):
    ux, uy, uz = axis
    st = np.sin(theta)
    ct = np.cos(theta)
    omct = 1 - ct
    R = np.array([[ct + ux*ux*omct, ux*uy*omct - uz*st, ux*uz*omct + uy*st],
            [ux*uy*omct + uz*st, ct + uy*uy*omct, uy*uz*omct - ux*st],
            [ux*uz*omct - uy*st, uy*uz*omct + ux*st, ct + uz*uz*omct]])
    return R


def scale_vec_to_abs( vec, value = 1.0 ):
    """Given vector, scale it to final value """
    assert isinstance( vec, np.ndarray )
    return vec * value / np.linalg.norm( vec )

def center_of_nuclei_charge( p_n, q_n ):
    assert isinstance( p_n, np.ndarray)
    assert isinstance( q_n, np.ndarray)
    coc = np.einsum( 'ij,i', p_n , q_n ) / np.sum(q_n )
    return coc

def electric_monopole_moment( q_e ):
    """given arrays"""
    return np.sum( q_e )

def electric_dipole_moment( p_n, q_n, p_e, q_e ):
    """given arrays """
    assert isinstance( p_e, np.ndarray)
    assert isinstance( q_e, np.ndarray)
    coc = center_of_nuclei_charge( p_n, q_n )
    return np.einsum('ij,i', (p_e - coc), q_e )

def electric_quadrupole_moment( p_n, q_n, p_e, q_e ):
    pass

def dipole( r, r_n, r_e ):
    """Return the total dipole moment given vector r,
    the positive charges at each point, and negative charges at each point

    For 
    
    """

    assert isinstance( r, np.ndarray )
    assert isinstance( r_n, np.ndarray )
    assert isinstance( r_e, np.ndarray )

#Center of nuclei charge chosen as gauge origin to set nuclei dipole moment to zero
    coc = np.einsum('ij,i', r, r_n )/ reduce(lambda a,x:a+x,r_n)
    p = np.einsum('ij,i', r - coc , r_n )

    return p

def center_and_xz(p1 = np.array([0,0,0]),
        p2 = np.array([0,0,1]),
        p3 = np.array([1,0,0]), ):
    """Given three points, will return one translation vector
    and 3 rotation matrices Rz**-1, Ry, and Rz**-1
    Required to put p1 in origo, with p2 and p3 in the xz plane and
    p2-p1 as the z axis.
    """
    t_v = -p1.copy()
    v2 = p2 - p1
    v3 = p3 - p1
    r1, r2, r3 = get_euler( v2, v3 )
    return t_v, r1, r2, r3

def get_t_and_rho(p1 = np.array([0,0,0]),
        p2 = np.array([0,0,1]),
        p3 = np.array([1,0,0]),
        plane = 'xz' ):
    """Given three points, will return one translation vector
    and 3 rotation matrices Rz**-1, Ry, and Rz**-1
    Required to put p1 in origo, with p2 and p3 in the xz plane and
    p2-p1 as the z axis.
    """
    t_v = -p1.copy()
    v2 = p2 - p1
    v3 = p3 - p1
    r1, r2, r3 = get_euler( v2, v3, plane = plane )
    return t_v, r1, r2, r3




def get_rotation( p1 = [0, 0, 0], p2 = [0, 0, 1], p3 = [1, 0, 0] ):
    """Give the rotation vector necessary to place points p1, p2, p3
    in cannonical cartesian unit basis
    
    with line p2-p1 as z axis, and line p3-p1 in the positive x-axis in zx-plane"""
    p1, p2, p3 = map( np.array, [p1, p2, p3] )
    r1 = p2 - p1
    r2 = p3 - p1

    ez = r1/ np.linalg.norm( r1 )
    y = np.cross( r1, r2 )
    ey = y / np.linalg.norm( y )
    ex = np.cross( ey, ez )
    return np.array( [ ex, ey, ez ] )

def get_euler( r1, r2, plane = 'xz' ):
    """Given two vectors, return 3 euler angles defined as follows
    
    1) Rotation around z axis clockwise,
    2) Rotation around y axis counter-clockwise,
    3) Rotation around z axis clockwise,

    Thus, the euler angles obtained for r1 and r2 will result in 
    r1 being the new z-axis (0,0,1) and r2 being the ne x-axis (1,0,0)

    get_euler( [0,0,1], [1,0,0] ) will thus produce 
    t1 = 0
    t2 = 0
    t3 = 0

    outputs vectors in order as they are needed to produce the input from
    the [0,0,1 and [1,0,0] unit vectors

    """

    rots = { 'xz' : (Rz_inv, Ry, Rz_inv ),
             'xy' : (Ry_inv, Rz, Ry_inv )
        }


    r1, r2 = map(lambda x: x.copy(), [r1, r2] )
# First angle
    t1 = np.arctan2( r1[1], r1[0] )
    if t1 < 0:
        t1 += 2 * np.pi

    r1, r2 = map( lambda x: np.einsum( 'ij,j', rots[ plane ][0]( t1 ), x), [r1, r2 ] )

# Second angle
    t2 = np.arctan2( -r1[0], r1[2] )
    if t2 < 0:
        t2 += 2 * np.pi

    r1, r2 = map( lambda x: np.einsum( 'ij,j', rots[ plane ][1]( t2 ), x), [r1, r2 ] )

# Third angle
    t3 = np.arctan2( r2[1], r2[0] )
    if t3 < 0:
        t3 += 2 * np.pi

    r1, r2 = map( lambda x: np.einsum( 'ij,j', rots[ plane ][2]( t3 ), x), [r1, r2 ] )
    return t3, t2, t1


def dipole_iso( d ):
    assert d.shape == (3,)
    return np.sqrt( (d**2).sum() )

def alpha_iso( a ):
    a = ut2s( a )
    assert a.shape == (3,3,)
    return a.trace()/3.0

def alpha_aniso( a ):
    a = ut2s( a )
    assert a.shape == (3,3,)
    return np.sqrt( (np.einsum('ij,ij', 3*a, a) - np.einsum('ii,jj', a, a )) /2 )

def alpha_aniso2( a ):
    a = ut2s( a )
    assert a.shape == (3,3,)
    tot = np.sqrt ((a[ 0, 0 ] - a[ 1, 1 ])**2 + (a[ 1, 1 ] - a[ 2, 2 ])**2 + (a[ 2, 2 ] - a[ 0, 0 ])**2 + 6*a[0, 2]**2 + 6*a[0, 1]**2 +6*a[1, 2]**2 )
    return np.sqrt(1/2.0)*tot




def b_para( b, p = None ):
    """Given beta in either UT or full tensor form, and dipole moment
    will return the projected beta vector, the rotationally invariant of the beta
    tensor"""
    b = ut2s( b )
    assert b.shape == (3,3,3)

    if p is not None:
        assert p.shape == (3,)
        b_p = 3.0/5.0*np.einsum('ijj,i', b, p )/np.linalg.norm( p )
    else:
        b_p = 0.0
        for i in range(3):
            b_p += b[2, i, i]
            b_p += b[i, 2, i]
            b_p += b[i, i, 2]
        b_p *= 1.0/5.0
    return b_p


def beta_vec( b ):
    b = ut2s( b )
    assert b.shape == (3,3,3)
    b_x =   b[0, 0, 0] +\
            b[0, 1, 1] + b[1, 0, 1] + b[1, 1, 0] +\
            b[0, 2, 2] + b[2, 0, 2] + b[2, 2, 0]

    b_y =   b[1, 1, 1] +\
            b[1, 0, 0] + b[0, 1, 0] + b[0, 0, 1] +\
            b[1, 2, 2] + b[2, 1, 2] + b[2, 2, 1]

    b_z =   b[2, 2, 2] +\
            b[2, 0, 0] + b[0, 2, 0] + b[0, 0, 2] +\
            b[2, 1, 1] + b[1, 2, 1] + b[1, 1, 2]
    return np.array( [b_x, b_y, b_z] )

def b_at_sphere( b, x, y, z ):
    """Given beta and 3 numpy arrays for meshed 3dgrid, 
    returns the beta dipole at each grid point """
    b_vec = np.zeros( x.shape + (3,) )
    for i, page in enumerate(x):
        for j, row in enumerate( page ):
            for k, col in enumerate( row ):
                tmp = x[ i, j, k], y[ i, j, k], z[ i, j, k]
                b_vec[i,j,k,:] = np.einsum( 'ijk,j,k', b, tmp, tmp )
    return b_vec


def E_at_sphere(r_min = 2.0,
        r_max = 5.0,
        r_points = 2,
        tau_points = 10,
        theta_points = 10 ):
    """Return x, y, z grid coordinate vectors along with vectors"""
    r = np.linspace( r_min, r_max, r_points )
    tau = np.linspace( 0, np.pi * 2, tau_points )
    theta = np.linspace( 0, np.pi, theta_points )
    R, TAU, THETA = np.meshgrid( r, tau, theta, )
    x = R * np.sin( THETA ) * np.cos( TAU )
    y = R * np.sin( THETA ) * np.sin( TAU )
    z = R * np.cos( THETA ) 


    return x, y, z


def polar_to_cartesian( r, tau, theta):
    x, y, z = r* np.sin( theta )*np.cos( tau ) \
           , r* np.sin(  theta )*np.sin( tau )  \
           , r* np.cos(  theta ) 
    return x , y , z

def splitter(arr, key = lambda x: x):
    """Given a list of items, return lists of different items based on key

    example simple
    -------
    >>> a = ['A', 'B', 'B' ]
    >>> splitter( a )
    [['A'], ['B', 'B']]


    For more advanced use, use lambdas

    example advanced
    -------
    >>> a = [ {'A': 1 }, {'A': 2}, {'A': 5} ]
    >>> splitter( a, lambda x: x['A'] < 3 )
    [[{'A': 1}, {'A': 2} ], [{'A': 5}]]


    """
    tmp = []
    have = []
    for i in arr:
        inew = key( i )
        if inew in have:
            for each in tmp:
                if key(each[0]) == inew:
                    each.append( i )
        else:
            have.append( inew )
            tmp.append( [] )
            tmp[-1].append( i )
    return tmp

def unique(arr, key = lambda x: x, get_original = False ):
    tmp = []
    have = []
    for i in arr:
        inew = key( i )
        try:
            if inew not in have:
                tmp.append( i )
                have.append( inew )
        except ValueError:
            if any( (inew != x).all() for x in have ):
                tmp.append( i )
                have.append( inew )
    if get_original:
        return tmp
    return have

def read_alpha( file_, freq = '0.0', in_AA = False, freqs = 1 ):
# If freqs > 1, will return a tuple of all alphas for each frequency
#
# Reading in Alpha tensor
    fre = freq[0:7]
    pat_alpha = re.compile(r'([XYZ])DIPLEN.*([XYZ])DIPLEN.*= *(-?\d*\.{1}\d+D*-?\+?\d*)')
    pat_new_freq = re.compile(r'FREQUENCY.*SECOND ORDER')
    alpha = np.zeros( [3,3,] )
    lab = ['X', 'Y', 'Z', ]

# For every new frequency, will append this one and store alpha in the last
# element, otherwise, first column is first frequency by default
    freqlist = None

    for i in open( file_ ).readlines():

        if pat_new_freq.search( i ):
            if freqlist is None:
                freqlist = []
            freqlist.append( np.zeros( (3,3 )) )

        if pat_alpha.search( i ):
            matched = pat_alpha.search(i).groups()
            if "D" in matched[2]:
                frac = float( matched[2].replace("D","E") )
            else:
                frac = float( matched[2] )

            A = matched[0]
            B = matched[1]

            alpha[ lab.index( A ) , lab.index( B ) ]  = frac
            freqlist[-1][lab.index( A ), lab.index( B ) ] = frac

            if A == "X" and B == "Y":
                alpha[ lab.index( B ) , lab.index( A ) ]  = frac
                freqlist[-1][lab.index( B ), lab.index( A ) ] = frac

            if A == "X" and B == "Z":
                alpha[ lab.index( B ) , lab.index( A ) ]  = frac
                freqlist[-1][lab.index( B ), lab.index( A ) ] = frac

            if A == "Y" and B == "Z":
                alpha[ lab.index( B ) , lab.index( A ) ]  = frac
                freqlist[-1][lab.index( B ), lab.index( A ) ] = frac

    if freqs > 1:
        return freqlist

    return alpha 

def read_energy( fname, calctype = 'HF' ):
    """Return the energy from dalton .out file fname"""

    for line in open(fname).readlines():
        if re.compile(r'.*Final.*energy').match(line):
            return line.split()[-1]
def converged( _str, ):
    pat = re.compile( r'Total wall time used in DALTON:' )
    for line in _str.split('\n'):
        if pat.search( line ):
            return True
    return False
def find_dir(name, path):
    for root, dirs, files in os.walk(path):
        if name in files:
            return root

def find_file(name, path):
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)


def upper_triangular(n, start=0):
    """Recursive generator for triangular looping of Carteesian tensor

Usage, form 2D-matrix from upper-triangular matrix represented by an array::

    ref = np.arange( 6 ) # Non-zero elements in 2-dimensional UT-tensor
    arr = np.zeros( (3, 3) ) # Target 
    for ind, (i, ii) in enumerate( upper_triangular(2) ):
        arr[ i, ii ] = ref[ ind ]
"""
    if n > 2:
        for i in range(start, 3):
            for j in upper_triangular(n-1, start=i):
                yield (i,) + j
    else:
        for i in range(start, 3):
            for j in range(i, 3):
                yield i, j

def ut_2_square( alpha ):
    assert len(alpha) == 6
    tmp_a = np.zeros( (3,3, ))
    for index, val in enumerate( upper_triangular(2) ) :
        tmp_a[ val[0], val[1] ] = alpha[ index ]
        tmp_a[ val[1], val[0] ] = alpha[ index ]
    return tmp_a

def ut_3_square( beta ):
    assert len(beta) == 10
    tmp_b = np.zeros( (3,3,3, ))
    for index, (i, j, k ) in enumerate( upper_triangular(3) ) :
        tmp_b[ i, j ,k] = beta[ index ]
        tmp_b[ i, k ,j] = beta[ index] 
        tmp_b[ j, i, k] = beta [ index ]
        tmp_b[ j, k, i] = beta [ index ]
        tmp_b[ k, i, j] = beta [ index ]
        tmp_b[ k, j, i] = beta [ index ]
    return tmp_b

def ut2s( vec ):
    """Transform upper triangular alpha or beta to square form"""
    if len( vec ) == 6:
        return ut_2_square( vec )
    elif len( vec ) == 10:
        return ut_3_square( vec )
    elif vec.shape == (3,3,):
        return vec
    elif vec.shape == (3,3,3,):
        return vec


def square_2_ut(alpha):
    assert alpha.shape == ( 3,3,)
    tmp_a = np.zeros( 6 )
    for index, (i, j ) in enumerate( upper_triangular(2) ):
        tmp_a[ index ] = (alpha[i, j] + alpha[ j, i]) / 2
    return tmp_a

def square_3_ut(beta):
    assert beta.shape == (3, 3, 3,)
    tmp_b = np.zeros( 10 )
    for index, (i, j, k ) in enumerate( upper_triangular(3) ):
        tmp_b[ index ] = ( \
                beta[i, j, k] + beta[i, k, j] + \
                beta[j, i, k] + beta[j, k, i] + \
                beta[k, i, j] + beta[k, j, i] )/ 6
    return tmp_b

def s2ut( vec ):
    if vec.shape == (3,3,):
        return square_2_ut( vec )
    elif vec.shape == (3,3,3,):
        return square_3_ut( vec )

def b_hrs_intensity( b):
    if b.shape == (10,):
        b = ut2s(b)
    elif b.shape != (3,3,3,):
        print "supplied wrong beta"
        raise SystemExit

    zzz = rot_avg( b )
    xzz = rot_avg( b, car1=0 )
    return np.sqrt( zzz + xzz )

def b_hrs_dep_ratio( b):
    if b.shape == (10,):
        b = ut2s(b)
    elif b.shape != (3,3,3,):
        print "supplied wrong beta"
        raise SystemExit

    zzz = rot_avg( b )
    xzz = rot_avg( b, car1=0 )
    return zzz/xzz

def rot_avg( beta, car1 = 2, car2 = 2, car3 = 2):
    """
    Requires euler.h5 binary file containing rotational angle products
    Define it as in current script directory + euler.h5
    """
    b_new = np.zeros( (3,3,3,) )
    try:
        import h5py
    except ImportError:
        raise SystemExit('Install h5py in order to calculate euler rotational averages')
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

def largest_triangle( _l ):
    """ Naive O( n^3 ) implementation of finding largest area triangle
    list _l have points and are thus shapes( n, 3 )"""

    largest = 0
    points = (0, 0, 0,)
    for p1, p2, p3 in itertools.product( _l, _l, _l ):
        if (p1 == p2).all() or (p1 == p3).all() or (p2 == p3).all():
            continue
        tri = np.linalg.norm( np.cross((p2-p1),(p3-p1)) ) /2
        if tri > largest:
            points = ( p1, p2, p3 )
            largest = tri
    return points

def file_is_used( _f ):
    """ Will return True if the file _f is open by an other process.

    Does not work for files opened in vim"""
    _f = os.path.realpath( _f )
    procs = [p for p in os.listdir( '/proc' ) if p.isdigit() ]
    all_fds = ["/proc/%s/fd" %p for p in procs ]

    user_fds = filter( lambda x: os.stat(x).st_uid == 1000, all_fds )
    all_files = map( lambda x: [ os.path.realpath(os.path.join(x, p)) for p in os.listdir(x) ] , user_fds )
    open_fds = reduce( lambda a, x: a+x, all_files )
    
    return _f in open_fds


# Special decorator used to cast time exception for functions taking more than 10 ms
class TimeException( Exception ):
    def __init__(self, val):
        self.val = val
def takes_time( func, *args, **kwargs ):
    def wrapped( *args, **kwargs ):
        d1 = time.time()
        val = func( *args, **kwargs )
        delta = time.time() - d1
        if delta > 0.01:
            raise TimeException( delta )
        return val
    return wrapped

def o_filter( 
        out_files, 
        vary = "r", 
        r = 0.0,
        tau = 0.0,
        theta = 0.0,
        rho1 = 0.0,
        rho2 = 0.0,
        rho3 = 0.0,
        ):
    p = re.compile(r'_(.*)-(.*)-(.*)-(.*)-(.*)-(.*).out')
    out = []
    for f in out_files:
        r_c, tau_c, theta_c, rho1_c, rho2_c, rho3_c = map(float, p.search(f).groups())
        if vary == "r":
            if (tau, theta, rho1, rho2, rho3) != (tau_c, theta_c, rho1_c, rho2_c, rho3_c):
                continue
        if vary == "tau":
            if (r, theta, rho1, rho2, rho3) != (r_c, theta_c, rho1_c, rho2_c, rho3_c):
                continue
        if vary == "theta":
            if (r, tau, rho1, rho2, rho3) != (r_c, tau_c, rho1_c, rho2_c, rho3_c):
                continue
        if vary == "rho1":
            if (r, tau, theta, rho2, rho3) != (r_c, tau_c, theta_c, rho2_c, rho3_c):
                continue
        if vary == "rho2":
            if (r, tau, theta, rho1, rho3) != (r_c, tau_c, theta_c, rho1_c, rho3_c):
                continue
        if vary == "rho3":
            if (r, tau, theta, rho1, rho2 ) != (r_c, tau_c, theta_c, rho1_c, rho2_c, ):
                continue
        out.append( f )
    return out

def r_scatt_int_nat( alpha ):
    """Given alpha tensor, returns the rayleigh scattering intensity,
    as scattered for natural light"""
    if alpha.shape == (6,):
        alpha = ut2s( alpha )
    assert alpha.shape == (3,3,)
    iso = alpha.trace()/3.0
    aniso = np.sqrt( (np.einsum('ij,ij', 3*alpha, alpha) - np.einsum('ii,jj',alpha,alpha )) /2 )
    return 45*iso**2 + 13*aniso**2

def r_scatt_int_pol_orth( alpha ):
    """Given alpha tensor, returns the rayleigh scattering intensity,
    as scattered for natural light"""
    if alpha.shape == (6,):
        alpha = ut2s( alpha )
    assert alpha.shape == (3,3,)
    iso = alpha.trace()/3.0
    aniso = np.sqrt( (3*np.einsum('ij,ij',alpha,alpha) -np.einsum('ii,jj',alpha,alpha )) /2 )
    return 45*iso**2 + 7*aniso**2

def r_scatt_int_pol_para( alpha ):
    """Given alpha tensor, returns the rayleigh scattering intensity,
    as scattered for natural light"""
    if alpha.shape == (6,):
        alpha = ut2s( alpha )
    assert alpha.shape == (3,3,)
    iso = alpha.trace()/3.0
    aniso = np.sqrt( (3*np.einsum('ij,ij',alpha,alpha) -np.einsum('ii,jj',alpha,alpha )) /2 )
    return 6*aniso**2



