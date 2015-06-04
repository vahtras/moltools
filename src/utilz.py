#!/usr/bin/env python


"""

This will replace read_dal.py as general functions module with no dependencies
on other modules or files in this repository
"""

import os,sys, re, argparse, ctypes, multiprocessing, functools
import numpy as np
import math as m

#from particles import *
from pd.gaussian import *

#import molecules 
#from template import Template
import matplotlib as mpl
#from use_calculator import *

import h5py


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
    au_to_nm = 45.563352491687866
    return "%.7f" %( au_to_nm / float(val) )

def make_para( shape = ( 0,) ):                                                      
    """Returns a np array that can be used for writing between different processes"""
    arr = np.zeros( shape )                                                          
    sb = multiprocessing.Array( ctypes.c_double, reduce(lambda a,b:a*b, arr.shape) )              
    sb = np.ctypeslib.as_array( sb.get_obj())                                        
    sb[:] = None                                                                     
    return sb.reshape( shape )  


def convex_hull_volume( pts):
    from scipy.spatial import Delaunay, ConvexHull
    def tetrahedron_volume(a, b, c, d):
        return np.abs(np.einsum('ij,ij->i', a-d, np.cross(b-d, c-d))) / 6
    ch = ConvexHull(pts)
    dt = Delaunay(pts[ch.vertices])
    tets = dt.points[dt.simplices]
    return np.sum(tetrahedron_volume(tets[:, 0], tets[:, 1],
        tets[:, 2], tets[:, 3]))

def rotate_point_by_two_points(p, p1, p2, theta):
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
def reflect_point_by_three_points( p, p1, p2, p3 ):
    """ will reflect point p by 
    plane formed by points p1, p2 and p3"""
    v1 , v2 = p2 - p1, p3 - p1
    n = np.cross( v1, v2 )
    n = n / np.linalg.norm( n )
    r = np.dot( p, n ) * n / np.linalg.norm( n )
    r_orig = p - r
    return -r + r_orig

def S( plane = 'xz'):
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



def scale_vec_to_abs( vec, value = 1.0 ):
    """Given vector, scale it to final value """
    assert isinstance( vec, np.ndarray )
    return vec * value / np.linalg.norm( vec )

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
    r3, r2, r1 = get_euler( v2, v3 )
    return t_v, r1, r2, r3



def get_euler( r1, r2 ):
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
    the z and x-axis

    """
    r1, r2 = map(lambda x: x.copy(), [r1, r2] )
# First angle
    t1 = np.arctan2( r1[1], r1[0] )
    if t1 < 0:
        t1 += 2 * np.pi

    r1, r2 = map( lambda x: np.einsum( 'ij,j', Rz_inv( t1 ), x), [r1, r2 ] )

# Second angle
    t2 = np.arctan2( -r1[0], r1[2] )
    if t2 < 0:
        t2 += 2 * np.pi

    r1, r2 = map( lambda x: np.einsum( 'ij,j', Ry( t2 ), x), [r1, r2 ] )

# Third angle
    t3 = np.arctan2( r2[1], r2[0] )
    if t3 < 0:
        t3 += 2 * np.pi

    r1, r2 = map( lambda x: np.einsum( 'ij,j', Rz_inv( t3 ), x), [r1, r2 ] )
    return t3, t2, t1

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
        if inew not in have:
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
def read_beta_hf( file_, freq = "0.0",  in_AA = False, out_AA = False ):
    nuc_dip = np.zeros(3)
    el_dip = np.zeros(3)
    alpha = np.zeros([3,3])
    beta = np.zeros([3,3,3])
    tmp = []
    atoms = []
    missing = {}
    exists = {}
    lab = ["X", "Y", "Z"]

    pat_Q = re.compile(r'Total charge of the molecule')
    pat_xyz = re.compile(r'^\s*(\w+)\s+(-*\d*\.+\d+)\s+(-*\d*\.+\d+)\s+(-*\d*\.+\d+) *$')
    pat_pol = re.compile(r'([XYZ])DIPLEN.*total.*:')
#Special xyz hack for camb3lyp output from akka dalton to find atoms
    pat_akka_xyz = re.compile(r'^\s*(\w+)\s+:\s+\d\s+x\s+(-*\d*\.+\d+)\s+\d\s+y\s+(-*\d*\.+\d+)\s+\d\s+z\s+(-*\d*\.+\d+) *$')

    pat_labels_xyz = re.compile(r'^\s*(\S+)\s+(-*\d*\.+\d+)\s+(-*\d*\.+\d+)\s+(-*\d*\.+\d+) *$')
# Reading in dipole and charge
    for i in open( file_ ).readlines():
        if pat_Q.search( i ):
            Q = float(i.split()[-1])
        if pat_xyz.match(i):
            f = pat_xyz.match(i).groups()
            matched = pat_xyz.match(i).groups()
#Skip coordinates in out file that are for MM region from QMMM
            kwargs = { "AA": in_AA, "element" :  matched[0], "x" : matched[1],
                    "y" : matched[2], "z" : matched[3] }
            tmpAtom = molecules.Atom( **kwargs )
            atoms.append( tmpAtom )

        elif pat_akka_xyz.match(i):
            f = pat_akka_xyz.match(i).groups()
            matched = pat_akka_xyz.match(i).groups()
#Skip coordinates in out file that are for MM region from QMMM
            kwargs = { "AA": in_AA, "element" :  matched[0], "x" : matched[1],
                    "y" : matched[2], "z" : matched[3] }
            tmpAtom = molecules.Atom( **kwargs )
            atoms.append( tmpAtom )

        elif pat_labels_xyz.match(i):
            f = pat_labels_xyz.match(i).groups()
            matched = pat_labels_xyz.match(i).groups()
            lab = matched[0]
            if len(lab.split('-')) == 4:
                element = "H"
            else:
                element = lab.split('-')[2][0]
            kwargs = { "AA": in_AA, "element" :  element, "x" : matched[1],
                    "y" : matched[2], "z" : matched[3] }
            tmpAtom = molecules.Atom( **kwargs )
            atoms.append( tmpAtom )

        if pat_pol.search(i):
            if pat_pol.search(i).group(1) == "X":
                try:
                    if "D" in i.split()[3]:
                        frac = float(i.split()[3].replace("D","E"))
                    else:
                        frac = float(i.split()[3])
                except IndexError:
                    if "D" in i.split()[2]:
                        frac = float( i.split()[2].strip(":").replace("D","E"))
                    else:
                        frac = float( i.split()[2].strip(":"))
                el_dip[0] += frac
            if pat_pol.search(i).group(1) == "Y":
                try:
                    if "D" in i.split()[3]:
                        frac = float(i.split()[3].replace("D","E"))
                    else:
                        frac = float(i.split()[3])
                except IndexError:
                    if "D" in i.split()[2]:
                        frac = float( i.split()[2].strip(":").replace("D","E"))
                    else:
                        frac = float( i.split()[2].strip(":"))
                el_dip[1] += frac
            if pat_pol.search(i).group(1) == "Z":
                try:
                    if "D" in i.split()[3]:
                        frac = float(i.split()[3].replace("D","E"))
                    else:
                        frac = float(i.split()[3])
                except IndexError:
                    if "D" in i.split()[2]:
                        frac = float( i.split()[2].strip(":").replace("D","E"))
                    else:
                        frac = float( i.split()[2].strip(":"))
                el_dip[2] += frac


#Set center of nuceli charge to 0
    coc = sum([ x.r * charge_dic[x.element] for x in atoms ]) /\
            sum([ charge_dic[x.element] for x in atoms ])
    for i in atoms:
        nuc_dip += charge_dic[ i.element ] * (i.r - coc )

    if in_AA:
# Make sure center of charge is in Atomic units to give correct electronic dipole
        coc /= a0

# Reading in Alfa and Beta tensor

    fre = str("%.5f" % float(freq))
    pat_alpha = re.compile(r'@.*QRLRVE.*([XYZ])DIPLEN.*([XYZ])DIPLEN.*%s' %fre)
    alpha = np.zeros( [3,3,] )
    lab = ['X', 'Y', 'Z', ]

    for i in open( file_ ).readlines():
        if pat_alpha.match( i ):
            try:
                if "D" in i.split()[-1]:
                    frac = float( i.split()[-1].replace("D","E") )
                else:
                    frac = float( i.split()[-1] )
            except IndexError:
                if "D" in i.split()[-1]:
                    frac = float( i.split()[-1].strip("=").replace("D","E") )
                else:
                    frac = float( i.split()[-1].strip("=") )
            A = pat_alpha.match(i).groups(1)[0]
            B = pat_alpha.match(i).groups(1)[1]
            alpha[ lab.index( A ) , lab.index( B ) ]  = frac
            if A == "X" and B == "Y":
                alpha[ lab.index( B ) , lab.index( A ) ]  = frac
            if A == "X" and B == "Z":
                alpha[ lab.index( B ) , lab.index( A ) ]  = frac
            if A == "Y" and B == "Z":
                alpha[ lab.index( B ) , lab.index( A ) ]  = frac

    pat_beta = re.compile(r'@ B-freq')
    for i in open( file_ ).readlines():
        if pat_beta.match(i):
            try:
                if i.split()[7].lstrip("beta") in exists:
                    continue
                exists[ i.split()[7].lstrip("beta") ] = float(i.split()[9] )
            except ValueError:
                a, b, c = i.split()[9].lstrip("beta").strip("()").split(",")
                if i.split()[7].lstrip("beta") in missing:
                    continue
                missing[ i.split()[7].lstrip("beta") ] =  "(%s;%s,%s)"%(a,b,c)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                try:
                    beta[i][j][k] = exists[ "(%s;%s,%s)" %(lab[i],lab[j],lab[k])]
                except KeyError:
                    beta[i][j][k] = exists[ missing["(%s;%s,%s)"%(lab[i],lab[j],lab[k]) ] ]
    N_el = sum([charge_dic[at.element] for at in atoms]) - Q
    tot_dip = el_dip - coc * N_el

    return atoms, tot_dip, alpha , beta

def read_props_qmmm( file_, freq = "0.0",  in_AA = False ):
    """ Same as read_beta_hf but skips coordinates not in allowd_elements
    """
    nuc_dip = np.zeros(3)
    el_dip = np.zeros(3)
    alpha = np.zeros([3,3])
    beta = np.zeros([3,3,3])
    tmp = []
    atoms = []
    missing = {}
    exists = {}
    lab = ["X", "Y", "Z"]

    pat_xyz = re.compile(r'^\s*(\w+)\s+(-*\d*\.+\d+)\s+(-*\d*\.+\d+)\s+(-*\d*\.+\d+) *$')
    pat_pol = re.compile(r'([XYZ])DIPLEN.*total.*:')

# Reading in dipole
    for i in open( file_ ).readlines():
        if pat_xyz.match(i):
            f = pat_xyz.match(i).groups()
            
            matched = pat_xyz.match(i).groups()
#Skip coordinates in out file that are for MM region from QMMM
            if matched[0] not in allowed_elements:
                continue

            kwargs = { "element" :  matched[0], "x" : matched[1],
                    "y" : matched[2], "z" : matched[3] }
            tmpAtom = molecules.Atom( **kwargs )
            atoms.append( tmpAtom )

        if pat_pol.search(i):
            if pat_pol.search(i).group(1) == "X":
                try:
                    if "D" in i.split()[3]:
                        frac = float(i.split()[3].replace("D","E"))
                    else:
                        frac = float(i.split()[3])
                except IndexError:
                    if "D" in i.split()[2]:
                        frac = float( i.split()[2].strip(":").replace("D","E"))
                    else:
                        frac = float( i.split()[2].strip(":"))
                el_dip[0] += frac
            if pat_pol.search(i).group(1) == "Y":
                try:
                    if "D" in i.split()[3]:
                        frac = float(i.split()[3].replace("D","E"))
                    else:
                        frac = float(i.split()[3])
                except IndexError:
                    if "D" in i.split()[2]:
                        frac = float( i.split()[2].strip(":").replace("D","E"))
                    else:
                        frac = float( i.split()[2].strip(":"))
                el_dip[1] += frac
            if pat_pol.search(i).group(1) == "Z":
                try:
                    if "D" in i.split()[3]:
                        frac = float(i.split()[3].replace("D","E"))
                    else:
                        frac = float(i.split()[3])
                except IndexError:
                    if "D" in i.split()[2]:
                        frac = float( i.split()[2].strip(":").replace("D","E"))
                    else:
                        frac = float( i.split()[2].strip(":"))
                el_dip[2] += frac

    for i in atoms:
        nuc_dip[0] += charge_dic[ i.element ] * i.x
        nuc_dip[1] += charge_dic[ i.element ] * i.y
        nuc_dip[2] += charge_dic[ i.element ] * i.z

# Reading in Alfa and Beta tensor

    fre = str("%.5f" % float(freq))
    pat_alpha = re.compile(r'@.*QRLRVE.*([XYZ])DIPLEN.*([XYZ])DIPLEN.*%s' %fre)
    alpha = np.zeros( [3,3,] )
    lab = ['X', 'Y', 'Z', ]

    for i in open( file_ ).readlines():
        if pat_alpha.match( i ):
            try:
                if "D" in i.split()[-1]:
                    frac = float( i.split()[-1].replace("D","E") )
                else:
                    frac = float( i.split()[-1] )
            except IndexError:
                if "D" in i.split()[-1]:
                    frac = float( i.split()[-1].strip("=").replace("D","E") )
                else:
                    frac = float( i.split()[-1].strip("=") )
            A = pat_alpha.match(i).groups(1)[0]
            B = pat_alpha.match(i).groups(1)[1]
            alpha[ lab.index( A ) , lab.index( B ) ]  = frac
            if A == "X" and B == "Y":
                alpha[ lab.index( B ) , lab.index( A ) ]  = frac
            if A == "X" and B == "Z":
                alpha[ lab.index( B ) , lab.index( A ) ]  = frac
            if A == "Y" and B == "Z":
                alpha[ lab.index( B ) , lab.index( A ) ]  = frac

    pat_beta = re.compile(r'@ B-freq')
    for i in open( file_ ).readlines():
        if pat_beta.match(i):
            try:
                if i.split()[7].lstrip("beta") in exists:
                    continue
                exists[ i.split()[7].lstrip("beta") ] = float(i.split()[9] )
            except ValueError:
                a, b, c = i.split()[9].lstrip("beta").strip("()").split(",")
                if i.split()[7].lstrip("beta") in missing:
                    continue
                missing[ i.split()[7].lstrip("beta") ] =  "(%s;%s,%s)"%(a,b,c)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                try:
                    beta[i][j][k] = exists[ "(%s;%s,%s)" %(lab[i],lab[j],lab[k])]
                except KeyError:
                    beta[i][j][k] = exists[ missing["(%s;%s,%s)"%(lab[i],lab[j],lab[k]) ] ]
    if in_AA:
        nuc_dip /= a0
    tot_dip = nuc_dip - el_dip

    return atoms, nuc_dip - el_dip, alpha , beta

def find_dir(name, path):
    for root, dirs, files in os.walk(path):
        if name in files:
            return root

def find_file(name, path):
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)

# Functions related to transforming between upper triangular and square form
def upper_triangular(n, start=0):
    """Recursive generator for triangular looping of Carteesian tensor"""
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
