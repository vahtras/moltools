#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""
The molecules modules serves as an interface to write water molecule input files using predefined geometries, to be used with the DALTON qm package.
"""

import re, os, itertools, warnings, subprocess, shutil, logging, string
import h5py

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt
import numpy as np
import cPickle as pickle
import copy as copymod


import read_dal
import utilz

from pd import gaussian
from template import Template
from generator import Generator

from loprop.loprop import *

a0 = 0.52917721092
au_nm_conv = 45.563352491
elem_array = ['X', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne']

charge_dict = {"H": 1.0, "C": 6.0, "N": 7.0, "O": 8.0, "S": 16.0,
        "P" : 15, "X" : 0.0 }
# from TIP3P charge defs.
el_charge_dict = {"H": .417, "O": -0.834 , "X" : 0.417 , 'S': -0.25}
mass_dict = {"H": 1.008,  "C": 12.0, "N": 14.01, "O": 15.999, "S": 32.066,
        "X" : 1.008, 'P' : 30.974 }

color_dict = { "X": 'black' ,"H":'brown', "N":'blue',"C":'green',"P":'black', "O":'red', 'S' : 'yellow'}

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

class Property( dict ):
    """
**An object representing properties as numpy.ndarray types mapped to by python dictionaries.**

**Supports up to quadrupoles, upper triangular polarizability and upper trianguler hyperpolarizability**

.. code:: python
    
    >>> p = Property()
    >>> print p["charge"]
    [0.0]

    >>> print p["dipole"]
    [0.0, 0.0, 0.0]

"""
    def __init__(self):

        self["charge"] = 0.0
        self["dipole"] = np.zeros( 3 )
        self["quadrupole"] = np.zeros( 6 )
        self["alpha"] =  np.zeros( 6 ) 
        self["beta"] =  np.zeros( 10 ) 

    def copy_property(self):
        p = Property()
        p["charge"] =      self["charge"]
        p["dipole"] =      self["dipole"].copy()
        p["quadrupole"] =  self["quadrupole"].copy()
        p["alpha"] =       self["alpha"].copy()
        p["beta"] =        self["beta"].copy()
        return p

    def __add__(self, other):
        tmp = Property()
        for i, prop in enumerate(self):
            tmp[prop] = np.array( self[prop] ) + np.array(other[prop] )
        return tmp
    def __sub__(self, other):
        assert isinstance( other, Property)
        tmp = Property()
        for i, prop in enumerate(self):
            tmp[prop] = np.array( self[prop] ) - np.array(other[prop] )
        return tmp

    @property
    def q(self):
        return self['charge']
    @property
    def d(self):
        return self['dipole']
    @property
    def Q(self):
        return self['quadrupole']
    @property
    def a(self):
        return self['alpha']
    @property
    def b(self):
        return self['beta']
    @q.setter
    def q(self, val):
        self['charge'] = val
    @d.setter
    def d(self, val):
        assert val.shape == (3,)
        self['dipole'] = val
    @Q.setter
    def Q(self, val):
        assert val.shape == (6,)
        self['quadrupole'] = val
    @a.setter
    def a(self, val):
        assert val.shape == (6,)
        self['alpha'] = val
    @b.setter
    def b(self, val):
        assert val.shape == (10,)
        self['beta'] = val

    @property
    def b_proj(self):
        """
        Rotationally invariant property

        Beta projected on the dipole moment vector for a whole molecule / segment.
        Should not be used if only for an atom
        """
        return np.einsum('ijj,i', utilz.ut2s(self.b), self.d)/np.linalg.norm(self.d) #* #self.d / np.linalg.norm( self.d )

    def potline(self, max_l =2 , pol= 22, hyper=1, fmt = "%.5f "):
        string = ""
        if 0  <= max_l :
            string += fmt % self["charge"]
        if max_l >= 1 :
            string += fmt*3 %( self["dipole"][0], self["dipole"][1], self["dipole"][2] )
        if max_l >= 2 :
            string += fmt*6  %( 
                    self["quadrupole"][0], self["quadrupole"][1], self["quadrupole"][2] ,
                    self["quadrupole"][3], self["quadrupole"][4], self["quadrupole"][5] )
        if pol == 1:
            string += fmt %( 
                    float(self["alpha"][0] + self["alpha"][3] + self["alpha"][5])/3,
                    )
        elif pol %10 == 2 :
            string += fmt * 6 %( 
                    self["alpha"][0], self["alpha"][1], self["alpha"][2] ,
                    self["alpha"][3], self["alpha"][4], self["alpha"][5] )
        if hyper == 1:
            string += fmt*10 %( 
                    self["beta"][0], self["beta"][1], self["beta"][2] ,
                    self["beta"][3], self["beta"][4], self["beta"][5] ,
                    self["beta"][6], self["beta"][7], self["beta"][8] ,
                    self["beta"][9])
        return string


    @staticmethod
    def from_propline( st, maxl = 2, pol = 22, hyper = 2 ):
        """
Given dalton POT, returns class Property that can be attached to Atom.

Convinience function for generating properties for the class Molecule directly by 
invoking dalton on a supercomputer.

    >>> p = Property.from_propline( "1 0.0 0.0 0.0 -0.25", maxl = 0 )
    >>> at.Property = p
        """
        st = map( float, st.split()[4:] )
        p = Property()

        p['charge'] = st.pop(0)
        if maxl > 0:
            for i in range(3):
                p['dipole'][i] = st.pop(0)
        if maxl > 1:
            for i in range(6):
                p['quadrupole'][i] = st.pop(0)

        if pol == 1:
            iso = st.pop(0)
            p['alpha'][0] = iso
            p['alpha'][4] = iso
            p['alpha'][6] = iso
        elif pol%10 == 2:
            for i in range(6):
                p['alpha'][i] = st.pop(0)
        if hyper == 2:
            for i in range(10):
                p['beta'][i] = st.pop(0)
        return p

    @staticmethod
    def from_template( at_string, template ):
        """Given string for atomic label, and dictionary template for the
        molecule,
        will return all properties found in template.py for this template
        """
        all_props = [ 'charge', 'dipole', 'quadrupole', 'alpha', 'beta' ]


        for p in all_props:
            if (at_string, p ) not in template:
                raise RuntimeWarning("'( %s, %s )' not found in provided template" %(at_string,p))



        p = Property() 
        for key in template:
            if key[0] == at_string:
                for each in all_props:
                    p[ each ] = template[ (at_string, each ) ]
        return p


    @staticmethod
    def add_prop_from_template( at, wat_templ ):

        """
Puts properties read from the :ref:`template` module into the :ref:`atom` at.

    
    >>> #Dist = True is the default, properties obtained using LoProp
    >>> temp = template.Template().get( dist = False ) 
    >>> w = Water.get_standard() 
    >>> Property.add_prop_from_template( w.o, temp )
    >>> print w.o.Property["dipole"]
    [0.0, 0.0, 0.78719]

"""
        p = Property()
        for i, keys in enumerate( wat_templ ):
            if keys[0] == ( at.element + str(at.order) ):
                p[keys[1]] = np.array( wat_templ[ keys ] )
        at.Property = p
#backwards compatible fix since charge was len 1 ndarray and now scalar
        at.Molecule.Property = True

    def inv_rotate( self, t1, t2, t3 ):
        """Rotate all properties by t1, t2, t3
        t1 negative rotation around Z-axis
        t2 positiv rotation around Y-axis
        t3 negative rotation around Z-axis
        """
        p = Property()
        r1 = utilz.Rz_inv(t1)
        r2 = utilz.Ry(t2)
        r3 = utilz.Rz_inv(t3)
        p.q = self.q
        p.d = np.einsum('ab,bc,cd,d', r3, r2, r1, self.d )
        p.a = utilz.s2ut( np.einsum('ec,fd,ca,db,ai,bj,ij', r3, r3, r2, r2, r1, r1, utilz.ut2s(self.a) ) )
        p.Q = utilz.s2ut( np.einsum('ec,fd,ca,db,ai,bj,ij', r3, r3, r2, r2, r1, r1, utilz.ut2s(self.Q) ) )
        p.b = utilz.s2ut( np.einsum('Id,Je,Kf,da,eb,fc,ai,bj,ck,ijk', r3, r3, r3, r2, r2, r2, r1, r1, r1, utilz.ut2s(self.b) ) )
        return p

    def rotate( self, t1, t2, t3 ):
        """Rotate all properties by t1, t2, t3
        t1 positive rotation around Z-axis
        t2 negative rotation around Y-axis
        t3 positive rotation around Z-axis
        """
        p = Property()
        r1 = utilz.Rz(t1)
        r2 = utilz.Ry_inv(t2)
        r3 = utilz.Rz(t3)
        p.q = self.q
        p.d = np.einsum('ab,bc,cd,d', r3, r2, r1, self.d )
        p.a = utilz.s2ut( np.einsum('ec,fd,ca,db,ai,bj,ij', r3, r3, r2, r2, r1, r1, utilz.ut2s(self.a) ) )
        p.Q = utilz.s2ut( np.einsum('ec,fd,ca,db,ai,bj,ij', r3, r3, r2, r2, r1, r1, utilz.ut2s(self.Q) ) )
        p.b = utilz.s2ut( np.einsum('Id,Je,Kf,da,eb,fc,ai,bj,ck,ijk', r3, r3, r3, r2, r2, r2, r1, r1, r1, utilz.ut2s(self.b) ) )
        return p

    def transform_by_matrix(self, matrix):
        """docstring for by_matrix"""
        assert matrix.shape == (3,3,)
        p = Property()
        p.q = self.q
        p.d = np.einsum( 'ij,j', matrix, self.d )
        p.a = utilz.s2ut(np.einsum( 'ai,bj,ij', matrix, matrix, utilz.ut2s(self.a) ))
        p.Q = utilz.s2ut(np.einsum( 'ai,bj,ij', matrix, matrix, utilz.ut2s(self.Q) ))
        p.b = utilz.s2ut(np.einsum( 'ai,bj,ck,ijk', matrix, matrix, matrix, utilz.ut2s(self.b) ))
        return p


    def transform_ut_properties( self, t1, t2, t3):
        """
Rotate all the properties of each atom by 3 euler angles.

    >>> w = Water.get_standard()
    >>> w.rotate( 0, np.pi/2, 0 )  #Rotate counter-clockwise by 90 degrees around y-axis
    >>> temp = template.Template().get() #Default template
    >>> Property.add_prop_from_template( w.o, temp )
    >>> print w.o.Property[ "dipole" ]
    array([ 0.   , -0.   ,  0.298])

#Dipole moment of oxygen atom pointing in positive z-direction

    >>> r1, r2, r3 = w.get_euler()
    >>> w.o.Property.transform_ut_properties( r1, r2, r3 )
    >>> print w.o.Property[ "dipole" ]
    [ -2.98000000e-01   3.64944746e-17   1.82472373e-17]

#Dipole moment of oxygen atom now pointing in negative x-direction

"""


        if self.has_key( "dipole" ):
            self["dipole"] = Rotator.transform_1( self["dipole"] , t1, t2, t3 )
        if self.has_key( "quadrupole" ):
            self["quadrupole"] = self.transform_ut_2( self["quadrupole"], t1, t2, t3 )
        if self.has_key( "alpha" ):
            self["alpha"] = self.transform_ut_2( self["alpha"],t1, t2, t3 )
        if self.has_key( "beta" ):
            self["beta"] = self.transform_ut_3( self["beta"], t1, t2, t3 )

    def transform_ut_2( self, prop, t1, t2 ,t3 ):
        tmp = Rotator.ut_2_square( prop )
        tmp = Rotator.transform_2( tmp , t1 ,t2 ,t3 )
        tmp = Rotator.square_2_ut( tmp )
        return tmp

    def transform_ut_3( self, prop, t1, t2 ,t3 ):
        tmp = Rotator.ut_3_square( prop )
        tmp = Rotator.transform_3( tmp, t1 ,t2 ,t3 )
        tmp = Rotator.square_3_ut( tmp )
        return  tmp 


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
        for index, (i, j ) in enumerate( upper_triangular(2) ):
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
        for index, (i, j, k ) in enumerate( upper_triangular(3) ):
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
        for index, val in enumerate( upper_triangular(2) ) :
            tmp_a[ val[0], val[1] ] = alpha[ index ]
            tmp_a[ val[1], val[0] ] = alpha[ index ]
        return tmp_a

    @staticmethod
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


class Atom(object):

    """
    **Object representation of atoms.**
    """
    def __init__(self, *args, **kwargs ):
        """
Initialize either directly:

.. code:: python

    >>> a = Atom( x = 0, y = 0, z = 0, 'element' = 'H', AA = True )

... or with pre-defined key-word arguments:

.. code:: python

    >>> kwargs = { 'x' : 0, 'y' : 0, 'z' : 0, 'element' : 'H', "AA" : True }
    >>> a = Atom( **kwargs )

List of key-word arguments:

======== ======== ========
Keyword  Default  Type
======== ======== ========
x        0.0      float
y        0.0      float
z        0.0      float
element  X        string
name     1-XXX-X1 string
pdb_name  X1       string
number   0        int
AA       True     bool
======== ======== ========

        """
#Element one-key char
        self._element = "X"

#Order in xyz files
        self._order = None

#Name is custom name, for water use O1, H2 (positive x ax), H3
        self._name = None

#Label is custom name, for water use O1, H2 (positive x ax), H3
        self._label = None

        self.x = 0.0
        self.y = 0.0
        self.z = 0.0


# Use populate_bonds in class Molecule to attach all atoms to their neighbours
        self.bonds = {}
        self.angles = {}
        self.dihedral = {}

        self._q = None


        self._number = None
        self._res_id = 0
        self._atom_id = None

        self.in_water = False

#Connectivity
        self.Molecule = False
        self.Cluster = None

#QM border settings
        self._in_qm = False
        self._in_mm = False
        self._in_qmmm = False

#Property set to true if atoms have properties
        self.Property = Property()
        self.AA = True

        if kwargs != {}:
            self.AA = bool( kwargs.get( "AA", False ) )
            self.x = float( kwargs.get( "x", 0.0 ))
            self.y = float( kwargs.get( "y", 0.0 ))
            self.z = float( kwargs.get( "z", 0.0 ))
            self.element = kwargs.get( "element", "X" )
            self.number = kwargs.get( "number", 0 )
            self.pdb_name = kwargs.get( "pdb_name", 'X1' )
            self.order = kwargs.get( "order", 0 )
            self.in_qm = kwargs.get( "in_qm", False )
            self.in_mm = kwargs.get( "in_mm", False )
            self.in_qmmm = kwargs.get( "in_qmmm", False )
            self._res_id = kwargs.get( "res_id", 0 )
        self._mass = None

# property setters and getters
    @property
    def element(self):
        if self._element:
            return self._element
    @element.setter
    def element(self, val):
        self._element = val
    @property
    def name(self):
        if self._name:
            return self._name
        if self._element and self._order:
            return self._element + str(self._order)
        return 'X'
    @name.setter
    def name(self,val):
        self._name = val
    @property
    def label(self):
        if self._label is not None:
            return self._label
        self._label = self.element
        return self._label
    @label.setter
    def label(self, val):
        self._label = val

    @property
    def order(self):
        if self._order:
            return self._order
        return 0
    @order.setter
    def order(self, val):
        self._order = val

    @property
    def com(self):
        return self.r
    def get_mol_line(self):
        return "{0:15s}{1:10.5f}{2:10.5f}{3:10.5f}\n".format( self.label, self.x, self.y, self.z ) 

    def in_mm(self):
        return self.Molecule.in_mm
    def in_qm(self):
        return self.Molecule.in_qm

    def pdb_string(self):
        x, y, z = map( lambda x: string.rjust( "%.3f"%x, 8)[:8], self.r )
        """Return pdb line as specified by the PDB standard"""
        st = "{0:6s}{1:5s}{2:1s}{3:4s}{4:1s}{5:3s}{6:1s}{7:1s}{8:4s}{9:1s}{10:3s}{11:8s}{12:8s}{13:8s}{14:22s}{15:2s}\n".format( "ATOM", 
                str(self.atom_id),
                "",
                self.pdb_name,
                "",
                self.Molecule.res_name,
                "",
                self.Molecule.Cluster._chain_id,
                string.rjust( str(self.Molecule.res_id), 4),
                "",
                "",
                x,
                y,
                z,
                "",
                string.rjust(self.element,2)
                )        
        return st

    @property
    def atom_id(self):
        if self._atom_id is not None:
            return self._atom_id
        self._atom_id = self.Molecule.index( self )
        return self._atom_id

    @property
    def p(self):
        """Wrapper to access class Property object attached to the atom"""
        return self.Property
    @p.setter
    def p(self, val):
        assert isinstance( val, Property )
        self.Property = val

    def move_closer_to_atom(self, atom, final_value):
        """Given other atom, will move this one so that it is at final value"""
        vec = self.r - atom.r
        origin = self.r
        new = utilz.scale_vec_to_abs( vec, value = final_value )
        self.x, self.y, self.z = atom.r + new

    def plot(self ):
        """
Plot Atom in a 3D frame

.. code:: python

    >>> a = Atom( element = 'H' )
    >>> a.plot()
    
"""

#Plot water molecule in green and  nice xyz axis
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d' )
        ax.plot( [0, 1, 0, 0, 0, 0], [0,0 ,0,1,0,0], [0,0,0,0,0,1] )
        ax.text( 1.1, 0, 0, "X", color = 'red' )
        ax.text( 0, 1.1, 0, "Y", color = 'red' )
        ax.text( 0, 0, 1.1, "Z", color = 'red' )

        ax.plot( [self.x], [self.y], [self.z], self.Molecule.style[self.element], linewidth= self.Molecule.linewidth[self.element] )
        ax.set_zlim3d( -5,5)
        plt.xlim(-5,5)
        plt.ylim(-5,5)
        plt.show()



    def __len__( self ):
        return 1
    def __iter__(self):
        yield self

    def scale_angle(self, scale = 1.0):
        """scales only angle
        
        defaults to 1.0"""

        if len(self.angles) > 1:
            warnings.warn("Did not scale %s since it had %d angles" %(self,len(self.angles)), Warning)
            return
        Rz, Rzi = Rotator.get_Rz, Rotator.get_Rz_inv
        Ry, Ryi = Rotator.get_Ry, Rotator.get_Ry_inv

        for at2, at3 in self.angles:
            r3 = self.bonds[at2].bonds[at3].r - self.bonds[at2].r 
            t = self.bonds[at2].r.copy()
            for i in self.Molecule:
                i.x, i.y, i.z = i.r - t

            rho1 = np.math.atan2( r3[1], r3[0] )
            for i in self.Molecule:
                r3 = np.dot( Rzi(rho1), r3 )
                i.x, i.y, i.z = np.dot( Rzi(rho1), i.r )

            rho2 = np.math.atan2( -r3[0], r3[2] )
            for i in self.Molecule:
                r3 = np.dot( Ry(rho2), r3 )
                i.x, i.y, i.z = np.dot( Ry(rho2), i.r )

            rho3 = np.math.atan2( r3[1], r3[0] )
            for i in self.Molecule:
                i.x, i.y, i.z = np.dot( Rzi(rho3), i.r )
            theta = scale*self.get_angle( self.bonds[at2], self.angles[(at2,at3)] )

            bond = self.get_bond( self.bonds[at2] )
            self.x = bond * np.sin(theta)
            np.testing.assert_almost_equal( self.y, 0 )
            self.z = bond* np.cos( theta)

            for i in self.Molecule:
                i.x, i.y, i.z = np.dot( Rz(rho3), i.r )
                i.x, i.y, i.z = np.dot( Ryi(rho2), i.r )
                i.x, i.y, i.z = np.dot( Rz(rho1), i.r )
                i.x, i.y, i.z = i.r + t

    def get_bond(self, other):
        return np.linalg.norm(other.r - self.r)

    def get_angle(self, at1, at2 ):
        r1 = self.r - at1.r
        r2 = at2.r - at1.r
        n1 = np.linalg.norm(r1)
        n2 = np.linalg.norm(r2)
        deg = np.arccos(np.dot(r1,r2)/(n1*n2)) 
        return deg

    def scale_bond(self, scale = 1.0):
        """scales only bond by a scalefactor 
        
        scale defaults to 1.0"""

        if len(self.bonds) > 1:
            warnings.warn("Did not scale %s since it had %d bonds" %(self,len(self.bonds)), Warning)
            return
        for at in self.bonds:
            self.translate( self.bonds[at].r + (self.r - self.bonds[at].r)*scale )

    def copy(self):
        return self.copy_atom()
    def copy_self(self):
        return self.copy_atom()
    def copy_atom(self):
        a = Atom( **{'x':self.x, 'y':self.y, 'z':self.z,'AA':self.AA,
            'element':self.element,'name':self.name,'number':self.number,
            'pdb_name':self.pdb_name} )
        a._res_id = self.res_id
        a._atom_id = self.atom_id
        a.Property = self.Property.copy_property()
        return a

    @property
    def r(self):
        return np.array( [ self.x, self.y, self.z ] )

    @property
    def q(self):
        if self.Property:
            return self.Property["charge"]
        if self._q is not None:
            return self._q
        self._q = el_charge_dict[ self.element ]
        return self._q

    @property
    def mass(self):
        if self._mass is not None:
            return self._mass
        self._mass = mass_dict[ self.element ]
        return self._mass

    @property
    def res_id(self):
        if self.Molecule:
            return self.Molecule.res_id
        if self._res_id:
            return self._res_id
        return 0

    def potline(self, max_l=2, pol=22, hyper=1):
        return  "{0:4} {1:10f} {2:10f} {3:10f} ".format( \
                str(self.Molecule.cluster_order), self.x, self.y, self.z ) + self.Property.potline( max_l, pol, hyper ) + "\n"

#Atom string method
    def __str__(self):
        return "%s %f %f %f" %(self.label, self.x, self.y, self.z)

    def __add__(self, other):
        return Molecule( self, other )

    def dist_to_atom(self, other):

        """
Return the distance between two atoms

.. code:: python

   >>> H1 = Atom( z = 0 )
   >>> H2 = Atom( z = 1 )
   >>> print H1.dist_to_atom( H2 )
   1.0

"""
        return np.sqrt( (self.x - other.x)**2 + (self.y -other.y)**2 + (self.z -other.z)**2 )



    def dist_to_point(self, other):
        """
Return the distance to a point

.. code:: python

   >>> a = Atom( z = 0 )
   >>> print H1.dist_to_point( [0, 3, 4] )
   5.0

"""
        return np.sqrt( (self.x - other[0])**2 + (self.y -other[1])**2 + (self.z -other[2])**2 )

    def to_AU(self):
        if self.AA:
            self.x /= a0
            self.y /= a0
            self.z /= a0
            self.AA = False

    def to_AA(self):
        if not self.AA:
            self.x *= a0
            self.y *= a0
            self.z *= a0
            self.AA = True

class Molecule( list ):
    """
**Inherits list methods, specific molecules will inherit from this class.**
"""

    def __init__(self , *args, **kwargs):
        super( Molecule, self).__init__()
#Bond dict defined in angstromg, if molecule is in AU will be different later
        self.bonding_cutoff = { 
                ('X','X') : 1.0,
                ('H','H') : 0.8,
                ('H','C') : 1.101,
                ('H','N') : 1.1,
                ('H','O') : 1.1,
                ('H','P') : 1.1,
                ('H','S') : 1.3,
                ('C','C') : 1.66,
                ('C','N') : 1.60,
                ('C','O') : 1.5,
                ('C','P') : 2.0,
                ('C','S') : 2.0,
                ('N','N') : 1.5,
                ('N','O') : 1.5,
                ('N','P') : 1.5,
                ('N','S') : 1.5,
                ('O','O') : 1.5,
                ('O','P') : 2.0,
                ('O','S') : 1.75,
                ('P','P') : 1.5,
                ('P','S') : 2.0,
                ('S','S') : 2.1,
            }
        for key1, key2 in self.bonding_cutoff.keys():
            self.bonding_cutoff[ (key2, key1)] = self.bonding_cutoff[ (key1, key2) ]

# Dictionary with bonds
        self.bond_dict = {}

#center will be defined for all molecules after all atoms are added
#depends on which molecule
        self._res_name = None
        self._res_id = 0
        self._r = None
        self._com = None
        self.Cluster = None
        self.no_hydrogens = True

# This will be set True if attaching LoProp properties
        self.LoProp = False

# Used for printing template properties, defines which atoms are centered around zx plane
        self._origo_z_x = None

# For plotting different elements:
        self.style = { "X": 'ko' ,"H":'wo', "N":'bo',"C":'go',"P":'ko', "O":'ro',
                'S' : 'yo'}
        self.linewidth = {"X":25,"H":25, "N": 30, "C": 30, "O":40, "P" : 40,
                'S' : 45 }

# Make emptpy, beware that this causes molecules to give zero dipole momnet
# before template is loaded
        self.Property = None

#By default, in no region
        self._in_qm = False
        self._in_mm = False
        self._in_qmmm = False

#By default, AU 
        self.AA = False

#if supplied a dictionary with options, gather these in self.info
        self.info = {}

        if type(args) == tuple:
            if len(args) == 1:
                if type(args[0]) == list:
                    for i in args[0]:
                        self.add( i )
                else:
                    self.add( args[0] )
            else:
                for at in args:
                    self.add( at )

        if kwargs != {} :
            for i in kwargs:
                self.info[ i ] = kwargs[ i ]
            self.AA = kwargs.get( "AA" , False )

    def connect_everything(self):
        for a in self:
            a.Molecule = self

    def reflect(self, key = lambda x:(x[0].r, x[1].r, x[2].r),
            plane = 'xz'):
        """docstring for reflect"""
        p1, p2, p3 = key( self )
        origin = p1.copy()
        t, r1, r2, r3 = utilz.center_and_xz( p1, p2, p3 )
        self.t( -t )

        R1_inv = utilz.Rz_inv( r1 )
        R2     = utilz.Ry( r2 )
        R3_inv = utilz.Rz_inv( r3 )
        S = utilz.S( plane )
        R1 = np.einsum( 'ij->ji', R3_inv )
        R2_inv = np.einsum( 'ij->ji', R2 )
        R3 = np.einsum( 'ij->ji', R1_inv )

# Rotate each atom to fit key orientation, followed by plane reflection and rotation back
        for at in self:
            at.x, at.y, at.z = np.einsum( 'ab, bc, cd, de, ef, fg, gi, i', R3, R2_inv, R1, S, R3_inv, R2, R1_inv, at.r )

            at.p = at.p.transform_by_matrix( R1_inv )
            at.p = at.p.transform_by_matrix( R2 )
            at.p = at.p.transform_by_matrix( R3_inv )
            at.p = at.p.transform_by_matrix( S )
            at.p.transform_ut_properties( r3, r2, r1)

    def t(self, *args):
        """Wrapper function for self.translate_by_r"""
        return self.translate_by_r( *args )

# Slicing the molecule givees back a molecule, but only accessing one index gives atom
    def __add__(self, other):
        return Molecule( list.__add__( self, other))

    def __getslice__(self, i, j):
        return self.__getitem__(slice(i, j))
                
    def __getitem__(self, item):
        if isinstance( item, slice ):
            result = list.__getitem__(self, item)
            try:
                return Molecule(result)
            except TypeError:
                return result
        else:
            return super(Molecule,self).__getitem__( item )
                
    #@property
    #def origo_z_x(self):
    #    if self._origo_z_x is None:
    #        logging.error('Did not set orientation for molecule in xz plane!')
    #        raise SystemExit
    #    return self._origo_z_x

    #def center_zx(self, p1, p2, p3):
    #    """Given 3 atoms will set p1 as origo, p2 as z axis and p3 in zx-plane"""
    #    assert isinstance( p1, Atom )
    #    assert isinstance( p2, Atom )
    #    assert isinstance( p3, Atom )
    #    self._origo_z_x = (p1.r.copy,(),())

    def inv_rotate(self, t1, t2, t3):
        """rotate all atoms and properties as
        1) inverse Z rotation by t1
        2) positive Y rotation by t2
        3) inverse Z rotation by t3
        """
#Put back in original point
        r1 = utilz.Rz_inv(t1)
        r2 = utilz.Ry(t2)
        r3 = utilz.Rz_inv(t3)
        for at in self:
            at.x, at.y, at.z = np.einsum('ab,bc,cd,d', r3, r2, r1, at.r )
            at.p = at.p.inv_rotate( t1, t2, t3 )
    def rotate(self, t1, t2, t3):
        """Rotate atomss and properties by t1, t2, t3
        t1 positive rotation around Z-axis
        t2 negative rotation around Y-axis
        t3 positive rotation around Z-axis
        """
        r1 = utilz.Rz(t1)
        r2 = utilz.Ry_inv(t2)
        r3 = utilz.Rz(t3)
        for at in self:
            at.x, at.y, at.z = np.einsum('ab,bc,cd,d', r3, r2, r1, at.r )
            at.p = at.p.rotate( t1, t2, t3 )

    def template(self, max_l = 0, pol = 1, hyp = 0,
            label_func = lambda x: x.pdb_name ):
        """Write out the template with properties for molecule"""
        if len(label_func.func_code.co_varnames) != 1:
            print "Unsupported multi key function"
            raise SystemExit
        st_label = "_".join( label_func.func_code.co_names )
        st = "{\n"
        st += "'meta' : { 'label' : '%s', },\n" % st_label
        #st += "'origo' : '%s',\n" % self.origo_z_x[0]
        #st += "'z' : '%s',\n" % self.origo_z_x[1]
        #st += "'x' : '%s',\n}" % self.origo_z_x[2]

        for at in self:
            st += "( {0:5s}, {1:8s}) : {2:2.5f},\n".format( "'" +label_func(at)  +"'" , "'charge'", at.p.q )
            tmp = "( {0:5s}, {1:8s}) : [%s],\n"%(reduce(lambda a,x:a+x,map(lambda x: " {%d:1.5f}, " %x, range(2,5) )))

            st += tmp.format( "'" + label_func(at) + "'", "'dipole'", *(at.p.d.tolist()) )
            tmp = "( {0:5s}, {1:8s}) : [%s],\n"%(reduce(lambda a,x:a+x,map(lambda x: " {%d:1.5f}, " %x, range(2,8) )))

            st += tmp.format(  "'" + label_func(at) + "'", "'quadrupole'",*at.p.Q.tolist() )
            tmp = "( {0:5s}, {1:8s}) : [%s],\n"%(reduce(lambda a,x:a+x,map(lambda x: " {%d:1.5f}, " %x, range(2,8) )))

            st += tmp.format( "'" + label_func(at) + "'", "'alpha'",*at.p.a )
            tmp = "( {0:5s}, {1:8s}) : [%s],\n"%(reduce(lambda a,x:a+x,map(lambda x: " {%d:1.5f}, " %x, range(2,12) )))

            st += tmp.format( "'" + label_func(at) +"'", "'beta'", *at.p.b )

        st += '\n}'
        return st

    def reorder(self):
        for at in self:
            at.order = self.index(at)

    @property
    def res_id(self):
        if self._res_id is not None:
            return self._res_id
        else:
            return 1

    @property
    def res_name(self):
        if self._res_name is not None:
            return self._res_name
        self._res_name = "MOL"
        return self._res_name
    
    def exclists(self, max_len = None):
        """Gives the exclusion string for each atom/bond in molecule for dalton PE
        file

        max_len is needed if other molecules are involved which have more atoms than this.
        Set max_len to the length of largest MM-region residue so that zeros
        is appended to the smaller segments

        for a water molecule:
            water.exclists()
            [(2, 3), (1, 3), (1, 2)]

            water.exclists(max_len = 4)
            [(2, 3, 0), (1, 3, 0), (1, 2, 0)]
        
        """
        if max_len is None:
            max_len = len(self)
        tmp = [tuple([j]+[i.cluster_order for i in self if i.cluster_order != j ]) for j in map(lambda x:x.cluster_order,self)]
        diff = max_len - len(self)
        for app in range(diff):
            for ind, each in enumerate(tmp):
                tmp[ind] = each + (0,)
        return tmp

    def get_xyz_string(self, ):
        st = "%d\n\n" % len(self)
        for i in self:
            st += "{0:10s} {1:10f} {2:10f} {3:10f}\n".format(\
                    i.element, i.x,  i.y , i.z )
        return st
    def append(self, atom ):
        super( Molecule, self).append( atom )

    def add(self, item):
        if isinstance( item, Atom ):
            self.append( item )
        else:
            logging.warning('Passed %s instance to Molecule' %type(item) )

    def plot_2d(self, key = None, ):
        """Plots the 2D projetion of the molecule onto the y plane,
        PLAN: TODO:
        Later implement internal plane projection on arbitray two-vector plane
        """
        fig, ax = plt.subplots()
        #norm  = self.get_internal_plane()
        self.populate_bonds()

        norm = np.array( [0, 0, 1] )

        if key is None:
            key = lambda x: x.element

        for at in self:
            x, y = (at.r - np.einsum('i,i', norm, at.r)*norm)[:-1]
            ax.scatter( x, y, color = color_dict[at.element], linewidth=4 )
            ax.annotate(  key( at ) , (x, y ), (x, y+0.1) )
        
        x = map(lambda x: (x.r - np.einsum('i,i', norm, x.r)*norm)[0], self)
        y = map(lambda x: (x.r - np.einsum('i,i', norm, x.r)*norm)[1], self)

        #plot the bonds
#Plot bonds
        for each in self.bond_dict:
            for key in self.bond_dict[ each ]:
                ax.plot( [key.x, each.x],
                         [key.y, each.y],
                         color = 'black', linewidth = 0.25 )

        for at in self:
            pass
        ax.set_title( 'Projection of molecule on yx-plane' )
        ax.set_xlabel( 'x-axis' )
        ax.set_ylabel( 'y-axis' )

        ax.set_xlim( min(x) - 1.0, max(x) + 1.0 )
        ax.set_ylim( min(y) - 1.0, max(y) + 1.0 )

        plt.show()

    def get_internal_plane(self):
        """Returns the normal to the internal plane defined as the plane
        spanned by the two vectors which are between the two largest inter-atomic
        distance from each other in the x, y, or z direction"""
        minx, miny, minz = np.inf, np.inf, np.inf
        maxx, maxy, maxz = -np.inf, -np.inf, -np.inf
        prj = [np.array([i, j, k]) for i, j, k in zip([1, 0, 0],[0, 1, 0], [0, 0, 1] ) ]
        for a1 in self:
            dx, dy, dz = a1.r
            if dx < minx: minx = dx
            if dy < miny: miny = dy
            if dz < minz: minz = dz
            if dx > maxx: maxx = dx
            if dy > maxy: maxy = dy
            if dz > maxz: maxz = dz
        i1, i2 = map(lambda x:x[0], sorted( zip( range(3), [maxx-minx, maxy-miny, maxz-minz] ), key = lambda x: x[1])[1:] )

        return np.cross(prj[i1], prj[i2])

    def __str__(self):
        return "MOL-%d" %self.res_id

    @property
    def in_qm(self):
        return self._in_qm
    @property
    def in_mm(self):
        return self._in_mm
    @in_qm.setter
    def in_qm(self, val):
        self._in_qm = val
    @in_mm.setter
    def in_mm(self, val):
        self._in_mm = val

    @property
    def label(self):
        if self._label is not None:
            return self._label
        return self.element + self.order
 
    def save(self, fname = "molecule.p"):
        pickle.dump( self, open( fname, 'wb' ), protocol = 2 )

    @staticmethod
    def load(fname = 'molecule.p'):
        if not os.path.isfile( fname):
            raise IOError
        return pickle.load( open(fname, 'rb' ) )
    
    @property
    def b_proj(self):
        b, p = Rotator.ut_3_square(self.sum_property['beta']), self.sum_property['dipole']
        return np.einsum( 'ijj,i', b, p )/np.linalg.norm( p )

    def attach_properties(self, 
            model = "TIP3P",
            method = "HF",
            basis = "ANOPVDZ",
            loprop = True,
            freq = "0.0",
            euler_key = lambda x: (x[0].r, x[1].r, x[2].r),
            template_key = lambda x: x.pdb_name,
            force_template = False):
        """
Attach property for Molecule method, by default TIP3P/HF/ANOPVDZ, static
        """
        if isinstance(self, Water) and not force_template:
            template_key = lambda x: x.element + str(x.order)

        if isinstance(self, Water):
            euler_key = lambda x: (x.o.r, (x.h1.r-x.h2.r)/2 + x.h2.r, x.h1.r)

        templ = Template().get( *(model, method, basis, loprop, freq) )
        for at in self:
            at.p = Property.from_template( template_key(at), templ )
        t1, t2, t3 = self.get_euler( key = euler_key )
        for at in self:
            at.p.transform_ut_properties( t1, t2, t3 )
        if loprop:
            self.LoProp = True
        else:
            self.LoProp = False
        self.Property = True

    def dist_to_point( self , point ):
        return np.sqrt(np.sum((self.com - np.array(point))**2))

    def get_euler(self, key = lambda x: (x[0].r, x[1].r, x[2].r),
            rot_type = None):
        """Takes a function as an argument.

        By default, the defauling function picks the point of the first 3 atoms in the Molecule.

        If the euler angles are sought after for specific atomic types, supply:

        get_euler( key = lambda x: (x.N.r, x.C.r, x.CA.r)

        where .N, .C, and .CA are wrapper functions for retrieving atoms with pdb_name of respective atom

        For more verbosity:
        get_euler( key = lambda x: (x.get_atom_by_pdbname('N').r,
             x.get_atom_by_pdbname('C').r,
             x.get_atom_by_pdbname('CA').r)
        
        """
        if rot_type == 'water':
            key = lambda x: (x[0].r, x[2].r + (x[1].r-x[2].r)/2, x[1].r)
        try:
            p1, p2, p3 = key( self )
        except IndexError:
            logging.error('Tried to get euler on Molecule with too few atoms')
        t, r1, r2, r3 = utilz.center_and_xz( p1, p2, p3 )
        return r1, r2, r3

    def props_from_targz(self,
            f_ = None,
            tmpdir = None,
            maxl = 2,
            pol = 22,
            hyper = 2,
            decimal = 5,
            ):
        if tmpdir is None:
            tmpdir = "/tmp"
        if f_ is None:
            print "Supply .tar.gz file from dalton quadratic .QLOP calculation"
            return
        import tarfile
#Using Olavs external scripts
        tarfile.open( f_, 'r:gz' ).extractall( tmpdir )
        try:
            outpot = MolFrag( tmpdir = tmpdir,
                    max_l = maxl,
                    pol = pol,
                    pf = penalty_function( 2.0 ),
                    freqs = None,
                    ).output_potential_file(
                            maxl = maxl,
                            pol = pol,
                            hyper = hyper,
                            decimal = decimal,
                            )
        except:
            print tmpdir
        lines = [ " ".join(l.split()) for l in outpot.split('\n') if len(l.split()) > 4 ]
        if not len(lines) == len(self):
            print "Something went wrong in MolFrag output, check length of molecule and the molfile it produces"
            raise SystemExit
        for at, prop in zip(self, lines):
            at.Property = Property.from_propline( prop ,
                    maxl = maxl,
                    pol = pol,
                    hyper = hyper )

    def props_from_qm(self,
            tmpdir = None,
            dalpath = None,
            procs = 4,
            decimal = 5,
            maxl = 2,
            pol = 22,
            hyper = 2,
            method = 'hflin',
            basis = ['ano-1 2', 'ano-1 4 3 1', 'ano-2 5 4 1' ],
            dalexe = None,
            basdir = '/home/x_ignha/repos/dalton/basis',
            log = None,
            keep_outfile = False,
            freq = None
            ):
        """
        Will generate a .mol file of itself, run a DALTON calculation as a
        childed process, get the properties back and put them on all atoms.

        Might take long time for large residues.
        """
        if freq == None:
            freq = 0.0
        if log:
            logging.basicConfig( filename=log, level=logging.DEBUG )

        if tmpdir is None:
            print "Warning, must set tmpdir explicitly, since this function removes all files from the dalton generated real tmpdir."
            print "For security reasons will not run this if not set tmpdir"
            print "Explicitly by user"
            return

        remove_files = ['RSP_VEC', 'AOPROPER', 'SIRIFC' ]
        if 'lin' in method:
            hyper = 0
#Specific for triolith host, will remove in slurm environment leftover RSP
#files if they exist in tmp dir
        if os.environ.has_key( 'SLURM_JOB_NAME' ):
#Set allocated temporary directory
            tmpdir = os.environ['SNIC_TMP']
            for f_ in [f for f in os.listdir(tmpdir) if "RSP" in f]:
                if os.path.isfile( os.path.join( tmpdir, f_ ) ):
                    os.remove( os.path.join( tmpdir, f_) )
        else:
            tmpdir = os.path.join( tmpdir, str(os.getpid()) )
            if not os.path.isdir( tmpdir ):
                os.mkdir( tmpdir )

        dal = 'dalton.dal'
        mol = 'molecule.mol'
        dal_full, mol_full = map( lambda x: os.path.join( tmpdir, x ), [dal,mol])
        if method == 'hfqua':
            open( dal, 'w').write( Generator.get_hfqua_dal( ) )
        elif method == 'b3lypqua':
            open( dal, 'w').write( Generator.get_b3lypqua_dal( ) )
        elif method == 'hflin':
            if freq:
                open( dal, 'w').write( Generator.get_hflin_freq_dal( freq = freq ) )
            else:
                open( dal, 'w').write( Generator.get_hflin_dal( ) )
        elif method == 'b3lyplin':
            if freq:
                open( dal, 'w').write( Generator.get_b3lyplin_freq_dal( freq = freq ) )
            else:
                open( dal, 'w').write( Generator.get_b3lyplin_dal( ) )
        else:
            print "wrong calculation type specified"
            return
        open( mol, 'w').write( self.get_mol_string( basis = basis) )

#Make sure that the external dalton script copies the .out and .tar.gz
#files from calculation to current directory once child process finishes

        if dalexe is not None:
#On triolith modern dalton can only be run this custom way
            p = subprocess.Popen(['sbcast', dal,
                os.path.join( tmpdir , 'DALTON.INP')],
                stdout = subprocess.PIPE )
            out, err = p.communicate()
            p = subprocess.Popen(['sbcast', mol,
                os.path.join( tmpdir, 'MOLECULE.INP'), ],
                stdout = subprocess.PIPE )
            out, err = p.communicate()
        elif os.path.isfile( dalpath ):
            dalton = dalpath
        else:
            print "set env variable DALTON to dalton script, \
             or supply the script to props_from_qm directly as  \
             dalpath = <path-to-dalscript> "
            raise SystemExit

        if dalexe:
#Run as dalton executable directly in the dir with *.INP files
            os.chdir( tmpdir )
            p = subprocess.Popen([ 
                "WORK_MEM_MB=1024",
                "WRKMEM=$(($WORK_MEM_MB*131072))"
                "DALTON_TMPDIR=%s"%tmpdir,
                "BASDIR=%s" %basdir,
                "mpprun",
                "-np",
                "%d" %procs,
                dalexe], stdout = subprocess.PIPE )
            out, err = p.communicate()

            tar = "final.tar.gz"
            of = "DALTON.OUT"
            p = subprocess.Popen(['tar',
                'cvfz',
                'AOONEINT','AOPROPER','DALTON.BAS',
                'SIRIFC','RSPVEC','SIRIUS.RST',
                tar
                ],
                stdout = subprocess.PIPE)
        else:
#Run as dalton script
            p = subprocess.Popen([dalton, 
                '-N', str(procs), '-noarch', '-D', '-t', tmpdir,
                dal, mol
                ], stdout = subprocess.PIPE,
                )
            p.stdout
            lines_iterator = iter(p.stdout.readline, b"")
            for line in lines_iterator:
                logging.debug( line ) 
            out, err = p.communicate()


            of = "DALTON.OUT"
            tar = "dalton_molecule.tar.gz"
#Need to do this since late dalton scripts appends the tmp with seperate PID
            real_tmp = utilz.find_dir( of, tmpdir )
            of, tar = map( lambda x: os.path.join( real_tmp, x ), [of, tar ] )
            #print out, err, real_tmp
            #raise SystemExit

        if 'lin' in method:
            at, p, a = read_dal.read_alpha( of )
        else:
            at, p, a, b = read_dal.read_beta_hf( of )

#Using Olavs external scripts
        try:
            outpot = MolFrag( tmpdir = real_tmp,
                    max_l = maxl,
                    pol = pol,
                    pf = penalty_function( 2.0 ),
                    freqs = [ freq ],
                    ).output_potential_file(
                            maxl = maxl,
                            pol = pol,
                            hyper = hyper,
                            decimal = decimal,
                            )
        except UnboundLocalError:
            logging.error("Some error in LoProp, check so that the latest source is in PYTHONPATH")
            print real_tmp

        try:
            lines = [ " ".join(l.split()) for l in outpot.split('\n') if len(l.split()) > 4 ]
        except:
            logging.error("Might be that DALTON was not properly executed. If MPI version, make sure all proper mpi shared objects are in your path.")

        if not len(lines) == len(self):
            print "Something went wrong in MolFrag output, check length of molecule and the molfile it produces"
            raise SystemExit
        f_at = lambda x: map(float,x.get_mol_line().split()[1:])
        f_prop = lambda x: map(float,x.split()[1:4])

        assert len( self ) == len( lines )


        for at, prop in zip(sorted(self, key = f_at), sorted( lines, key = f_prop )):
            at.Property = Property.from_propline( prop ,
                    maxl = maxl,
                    pol = pol,
                    hyper = hyper )
        self.LoProp = True


#So that we do not pollute current directory with dalton outputs
#Also remove confliction inter-dalton calculation files
# For triolith specific calculations, remove all files in tmp
        for f_ in [os.path.join(real_tmp, f) for f in os.listdir(real_tmp) if os.path.isfile( os.path.join( real_tmp, f)) ]:
            os.remove( f_ )

        try:
            for f_ in [f for f in os.listdir(real_tmp) if os.path.isfile(f) ]:
                print os.path.join(real_tmp, f_) 
            for f_ in [mol, dal]:
                print f_
                os.remove( f_ )
        except OSError:
            logging.error('Something wrint in trying to remove files in real_tmp')
        if not keep_outfile:
            for fil in [f for f in os.listdir(os.getcwd()) if "dalton_molecule" in f]:
                os.remove( os.path.join( os.getcwd(), fil ) )
        #try:
        #    os.remove( tar )
        #except OSError:
        #    pass
        #os.remove( of )
        #shutil.rmtree( tmpdir )
        #for f_ in [mol, dal]:
        #    try:
        #        os.remove( f_ )
        #    except OSError:
        #        pass

    @classmethod
    def from_string(cls, fil):
        """Given .xyz file return a Molecule with these atoms"""
        rel = open(fil).readlines()[2:]
        m = cls()
        for i in range(len(rel)):
            m.append( Atom(**{'element':rel[i].split()[0],
                'x':rel[i].split()[1],
                'y':rel[i].split()[2],
                'z':rel[i].split()[3].strip(),
                'number' : i+1,
                }) )
        return m

    def custom_names(self):
        for i in self:
            i.name = i.element + str(i.number)

    def populate_bonds(self):
#Implement later that it can only be called once
        bond_dict = {}
        for i, at in enumerate( self ):
            bond_dict[ at ] = []

        if self.AA:
            conv = 1.0
        else:
            conv = 1/a0
        for i in range(len(self)):
            for j in range( i + 1, len(self) ):
                if self[i].dist_to_atom( self[j] ) < conv*self.bonding_cutoff[ (self[i].element, self[j].element) ]:

                    self[i].bonds[ self[j].name ] = self[j]
                    self[j].bonds[ self[i].name ] = self[i]
                    bond_dict[ self[i] ].append( self[j] )
                    bond_dict[ self[j] ].append( self[i] )
        self.bond_dict = bond_dict

    def populate_angles(self):
# Must be run after populate_bonds
        for at1 in self:
            for at2 in [at2 for at2 in at1.bonds.values()]:
                for at3 in [at3 for at3 in at2.bonds.values() if at3 is not at1 ]:
                    at1.angles[(at2.name,at3.name)] = at3

    @staticmethod
    def from_charmm_file( f):
        """Return molecule just by bonds from a charmm force field file
        
        Good for generation of dihedrals
        """
        m = Molecule()
        reg_at = re.compile(r'ATOM\s\w+\s+\w+\s+-*\d{1}.\d+\s')
        reg_bond = re.compile(r'BOND\s+')
        reg_el = re.compile(r'(^\w).*')

        ats = []
        for i in open(f).readlines():
            if reg_at.match( i ):
                m.append( Atom( **{ 'name':i.split()[1], 'element': reg_el.match(i.split()[1]).group(1) } ))

        for i in open(f).readlines():
            if reg_bond.match( i ):
                el1 = reg_el.match( i.split()[1]).group(1)
                el2 = reg_el.match( i.split()[2]).group(1)

                m.bond_dict[ (i.split()[1], i.split()[2] ) ] = m.bonding_cutoff[ \
                        (el1, el2) ]
                m.bond_dict[ (i.split()[2], i.split()[1] ) ] = m.bonding_cutoff[ \
                        (el2, el1) ]


        for i in range(len(m)):
            for j in range(len(m)):
                if m.bond_dict.has_key( (m[i].name, m[j].name) ) or \
                    m.bond_dict.has_key( (m[j].name, m[i].name) ) :
                    m[i].bonds.append( m[j] )
                    m[j].bonds.append( m[i] )

        dih = m.find_dihedrals()

        full_charm = ""
        skip = []

        aname_to_atype, atype_to_aname, atype_to_anumber, atype_dihed, anumber_to_atype = m.atom_map_from_string(f)

        for at in m:
            at.number = atype_to_anumber[ at.name ]


        for at in m:
            for targ in at.dihedral:
                l1 =  tuple( at.dihedral[targ] ) 
                l2 =  tuple( reversed(at.dihedral[targ] ) )
                t1 =  tuple( map( lambda x: aname_to_atype[x], at.dihedral[targ] ) )
                t2 =  tuple( reversed(map( lambda x: aname_to_atype[x], at.dihedral[targ] ) ))
                if t1 in atype_dihed:
                    if (l1 in skip) or (l2 in skip):
                        continue
                    skip.append( l1 )
                    pre_str = "\t".join( map( lambda x: str(atype_to_anumber[x]),
                        at.dihedral[targ] ))
                    full_charm += pre_str + "\t" + atype_dihed[ t1 ] + '\n'
                    continue
                if t2 in atype_dihed:
                    if (l2 in skip) or (l2 in skip):
                        continue
                    skip.append( l2 )
                    pre_str = "\t".join( map( lambda x: str(atype_to_anumber[x]),
                        at.dihedral[targ] ))
                    full_charm += pre_str + "\t" + atype_dihed[ t2 ] + '\n'
        return full_charm

    def find_dihedrals(self):

        dihed = []
        for at1 in self:
            if at1.bonds == []:
                continue
            for at2 in at1.bonds:
                if at2.bonds == []:
                    continue
                for at3 in [a for a in at2.bonds if a != at1]:
                    if at3.bonds == []:
                        continue
                    for at4 in [a for a in at3.bonds if a != at2]:
                        dihed.append( [at1.name, at2.name, at3.name, at4.name] )
                        at1.dihedral[ (at3.name, at4.name) ] = (at1.name,at2.name, at3.name, at4.name)
        return dihed

    @property
    def p(self):
        return self.sum_property

    @property
    def sum_property(self):
        """
Return the sum properties of all properties in molecules

.. code:: python
    >>> wat
        """
        el_dip = np.array([ (at.r-self.coc)*at.Property['charge'] for mol in self for at in mol])
        nuc_dip = np.array([ (at.r-self.coc)*charge_dict[at.element] for mol in self for at in mol])
        dip_lop = np.array([at.Property['dipole'] for mol in self for at in mol])
        dip = el_dip + nuc_dip
        d = (dip + dip_lop).sum(axis=0)
        p = Property()
        for at in self:
            p += at.Property
        p['dipole'] = d
        return p

#Vector pointing to center of atom position
    @property
    def r(self):
        """
Center of coordinate

.. code:: python

   >>> m = Molecule()
   >>> m.append( Atom(element = 'H', z = 1) )
   >>> m.append( Atom(element = 'O', z = 0) )
   >>> print m.r
   0.5

"""
        return  np.array([at.r for at in self]).sum(axis = 0) / len(self)

    def translate_by_r( self, *args ):
        """Will translate all atoms by vector r"""
        if type(args[0]) == int:
            r = np.array( args )
            assert len(r) == 3 
        elif type( args[0] ) == list:
            r = np.array( args[0] )
            assert len(r) == 3 
        else:
            r = args[0]
        for at in self:
            at.x += r[0]
            at.y += r[1]
            at.z += r[2]
        return self
    def translate(self, r):
        """
Translate molecules center-of-mass to position r

.. code:: python

    >>> m = Molecule()
    >>> m.append( Atom(element = 'H', z = 1) )
    >>> m.append( Atom(element = 'H', z = 0) )
    >>> print m.com
    [0, 0, 0.5 ]
    >>> m.translate( [0, 3, 5] )
    >>> print m.com
    [0, 3, 5 ]
    
"""
        vec = r - self.com
        for at in self:
            at.x = vec[0] + at.x 
            at.y = vec[1] + at.y 
            at.z = vec[2] + at.z 
        return self


    def translate_coc(self, r):
        vec = r - self.coc
        for at in self:
            at.x = vec[0] + at.x 
            at.y = vec[1] + at.y 
            at.z = vec[2] + at.z 
        return self

    @staticmethod
    def atom_map_from_string( fil ):
        aname_to_atype = {}
        atype_to_aname = {}
        aname_to_anumber = {}
        atype_dihedral_dict = {}
        anumber_to_atype = {}

        reg = re.compile(r'ATOM\s\w+\s+\w+\s+-*\d{1}.\d+\s')

        reg_dihed = re.compile (r'\w+\s+\w+\s+\w+\s+\w+\s+-*\d.\d+\s+\d{1}')

        cnt = 1
        for i in open(fil).readlines():
            if reg.match(i):
                aname_to_atype[ i.split()[1] ] = i.split()[2]
                atype_to_aname[ i.split()[2] ] = i.split()[1]
                aname_to_anumber[ i.split()[1] ] = cnt
                anumber_to_atype[ cnt ] = i.split()[1]
                cnt += 1
            if reg_dihed.match( i ):
                atype_dihedral_dict[(i.split()[0], i.split()[1],
                    i.split()[2], i.split()[3])] = " ".join( i.split()[4:] )
        return aname_to_atype, atype_to_aname, aname_to_anumber, atype_dihedral_dict, anumber_to_atype


#Center of nuclei charge
    @property
    def coc(self):
        """
Return center of charge

.. code:: python

    >>> m = Molecule()
    >>> m.add_atom( Atom( z : 0.11, element : 'H' ) )
    >>> m.coc
    [0., 0., 0.11]

        """
        return sum( [at.r * charge_dict[at.element] for at in self])\
                /sum( map(float,[charge_dict[at.element] for at in self]) )

    @property
    def com(self):
        return np.array([at.mass*at.r for at in self]).sum(axis=0) / np.array([at.mass for at in self]).sum()

    def dist_to_mol(self, other):
        """
Distance to other molecule, measured by center-of-mass

.. code:: python

    >>> m1 = Molecule( )
    >>> m2 = Molecule( )
    >>> m1.append( Atom()) ; m1.append( Atom( z = 1) )
    >>> m2.append( Atom(x = 1)) ; m2.append( Atom( x = 1, z = 1) )
    >>> print m1.dist_to_mol( m2 )
    1.0
    
"""
        return np.sqrt( ((self.com - other.com)**2 ).sum(axis=0) )


    def plot(self, copy = True, center = False, d = False ):
        """
Plot Molecule in a 3D frame

.. code:: python

    >>> m = Molecule()
    >>> m.append( Atom(element = 'H', x = 1, z = 1) )
    >>> m.append( Atom(element = 'H', x =-1, z = 1) )
    >>> m.append( Atom(element = 'O', z = 0) )
    >>> m.plot()
    
"""

#Make a copy in order to not change original, and perform plot on it
        if copy:
            copy = copymod.deepcopy( self )
        else:
            copy = self

        if center:
            copy.center()

#Plot water molecule in green and  nice xyz axis
        copy.populate_bonds()
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d' )
#Plot bonds
        for each in copy.bond_dict:
            for key in copy.bond_dict[ each ]:
                ax.plot( [key.x, each.x],
                         [key.y, each.y],
                         [key.z, each.z], color = 'black' )

        ax.plot( [0, 1, 0, 0, 0, 0], [0,0 ,0,1,0,0], [0,0,0,0,0,1] )
        ax.text( 1.1, 0, 0, "X", color = 'red' )
        ax.text( 0, 1.1, 0, "Y", color = 'red' )
        ax.text( 0, 0, 1.1, "Z", color = 'red' )
        if d:
            x = copy.coc[0]
            y = copy.coc[1]
            z = copy.coc[2]
            p = copy.p
            ax.plot( [x,x+p[0]], [y,y+p[1]], [z,z+p[2]], 'k-', linewidth = 3 )
            ax.plot( [p[0]],[p[1]],[p[2]],'ko', markersize = 5, linewidth = 5 )
        for i in copy:
            ax.plot( [i.x], [i.y], [i.z], copy.style[i.element], linewidth= copy.linewidth[i.element] )

        ax.set_zlim3d( -5,5)
        plt.xlim(-5,5)
        plt.ylim(-5,5)
        plt.show()

    def get_mol_string(self, basis = ("ano-1 2", "ano-1 4 3 1",
        "ano-2 5 4 1" ) ):
        if len( basis ) > 1:
            el_to_rowind = {"H" : 0, "C" : 1, "O" : 1, "N" : 1,
                    "S" : 2, "P" : 2}
        else:
            el_to_rowind = {"H" : 0, "C" : 0, "O" : 0, "N" : 0, "S" : 0 }
        st = ""
        s_ = ""
        if self.AA: s_ += " Angstrom"
        uni = utilz.unique([ at.element for at in self])
        st += "ATOMBASIS\n\n\nAtomtypes=%d Charge=0 Nosymm%s\n" %(len(uni), s_)
        for el in uni:
            st += "Charge=%s Atoms=%d Basis=%s\n" %( str(charge_dict[el]),
                    len( [all_el for all_el in self if (all_el.element == el)] ),
                    basis[ el_to_rowind[el] ])
            for i in [all_el for all_el in self if (all_el.element == el) ]:
                st += i.get_mol_line()
        return st

    def get_pdb_string(self):
        st = """"""
        for at in self:
            st += at.pdb_string()
        return st

    def copy_self(self):
        return self.copy()
    def copy(self):
        return copy.deepcopy(self)

    def get_inp_string(self, method ='B3LYP', basis = "6-31+g*", procs= 8):

        """Write gaussian .inp file for geometry optimization"""
        st = r"%" + "Nprocshared=%d\n" %procs
        st += r"%" + "Mem=20MW\n"
        st += "#p %s/%s opt " %(method,basis)
        if not self.AA:
            st += "units=au " 
        st += '\n\ncomment\n\n'
        st += "%d %d\n" %( self.q, 1 )
        for i in self:
            st += "%s %.5f %.5f %.5f\n" %(i.element, i.x, i.y, i.z ) 
        st+= '\n\n\n'
        return st

    def center(self):

        """
Center molecule with center-of-mass in origo

.. code:: python

    >>> m.com
    [0., 0., 1.,]

    >>> m.center()
    >>> m.com
    [0., 0., 0.,]

"""
        tmp = np.array( [0,0,0] )
        self.translate( tmp )

    @staticmethod
    def from_mol_file( molfile, in_AA = False, out_AA = False):
        """
Read in molecule given .mol file and unit specification.

.. code:: python

    >>> m = ( "water.mol", AA = True )
    >>> for at in m:
            print at.element
    H
    H
    O
    
"""
        pat_labels_xyz = re.compile(r'^\s*(\S+-+\S+)\s+(-*\d*\.+\d+)\s+(-*\d*\.+\d+)\s+(-*\d*\.+\d+) *$')
        pat_xyz = re.compile(r'^\s*(\w+)\s+(-*\d*\.+\d+)\s+(-*\d*\.+\d+)\s+(-*\d*\.+\d+) *$')
        tmp_molecule = Molecule( AA = in_AA )
        for i in open( molfile ).readlines():
            if pat_xyz.search(i):
                matched = pat_xyz.match(i).groups()
                pd = matched[0].split('-')[-1]
                kwargs = { "AA": in_AA, 
                        "element" : matched[0][0],
                        "name" :  matched[0], "x" : matched[1],
                        "y" : matched[2], "z" : matched[3], 
                        "pdb_name" : pd }
                tmpAtom = Atom( **kwargs )
                tmp_molecule.append( tmpAtom )
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
                tmpAtom = Atom( **kwargs )
                tmp_molecule.append( tmpAtom )
        if in_AA:
            if not out_AA:
                tmp_molecule.to_AU()
        return tmp_molecule


    @staticmethod
    def from_xyz_string( _string, in_AA = True, out_AA = True ):
        m = Molecule( AA = in_AA )
        for ind, i in enumerate( _string.split('\n') ):
            if ind in [0, 1]: 
                continue

            label = i.split()[0]
            if len(label) == 1:
                elem = label
            else:
                elem = label[0]
            x = i.split()[1]
            y = i.split()[2]
            z = i.split()[3]
            at = Atom( **{"element":elem,
                "x" : x,
                "y" : y,
                "z" : z,
                "AA" : in_AA,
#Order later used to read in templates
                "order" : ind - 1,
                'name' : elem + str(ind-1)
                })
            m.append( at )
            at.Molecule = m

        if in_AA:
            if not out_AA:
                m.to_AU()
        return m



    @staticmethod
    def from_xyz( f, in_AA = True, out_AA = True ):
        """
Read in molecule from .xyz file given unit specifications.
Resulting molecule will be in either atomic units [ out_AA = False ], or in 
Angstrom [ out_AA = True ]

.. code:: python

    >>> m = ( "water.mol", in_AA = True, out_AA = False )
    >>> for at in m:
            print at.z
    H
    H
    O
    
"""
        if not os.path.isfile( f ):
            raise IOError

        fil = open(f).readlines()
        m = Molecule( AA = in_AA )
        for ind, i in enumerate( fil ):
            if ind in [0, 1]: 
                continue

            elem = i.split()[0]
            x = i.split()[1]
            y = i.split()[2]
            z = i.split()[3]
            at = Atom( **{"element":elem,
                "x" : x,
                "y" : y,
                "z" : z,
                "AA" : in_AA,
#Order later used to read in templates
                "order" : ind - 1,
                'name' : elem + str(ind-1)
                })
            m.append( at )
            at.Molecule = m

        if in_AA:
            if not out_AA:
                m.to_AU()
        return m

    def to_AU(self):
        if self.AA:
            for at in self:
                at.to_AU()
            self.AA = False

    def to_AA(self):
        if not self.AA:
            for at in self:
                at.to_AA()
            self.AA = True

    def get_atom_by_pdbname(self, label, dup = False):
        at = []
        for i in self:
            if i.pdb_name == label:
                at.append(i)
        if len(at) > 1 and not dup:
            print "Warning: Duplicate with pdb name %s, returning %s" %(label, at[0].label )
            return at[0]
        elif len(at) > 1 and dup:
            return at

        elif len(at) == 0:
            print "No %s in %s" %(label, self)
            return
        return at[0]

class Water( Molecule ):
    """
**Derives all general methods from Molecule.**

**Specific for water is the get_euler method, which defines which water orientation is the reference position.**
"""

    def __init__(self , *args, **kwargs):
        super(Water, self).__init__( *args, **kwargs )
        self.atoms = 0

        self.no_hydrogens = True
        self.h1 = False
        self.h2 = False
        self.o  = False

        self.AA = False

        self._coc = None

        self.in_qm = False
        self.in_mm = False
        self.in_qmmm = False

        if kwargs is not {}:
            self.AA = kwargs.get( "AA", False )

    def copy(self):
        return self.copy_water()
    def copy_self(self):
        return self.copy_water()
    def copy_water(self):
        w = Water()
        [w.append(i.copy_atom()) for i in self]
        return w
                
#Must override the parent Molecule to return proper class type water
    def __getitem__(self, item):
        if isinstance( item, slice ):
            result = list.__getitem__(self, item)
            try:
                return Molecule(result)
            except TypeError:
                return result
        else:
            return super(Molecule,self).__getitem__( item )




    @staticmethod
    def get_standard( AA = False,
            model = 'tip3p',
            worst = False):
        """
Return water molecule from specified template with :math:`r=0.9572` Angstrom and 
:math:`\\theta=104.52` degrees.

.. code:: python

    >>> m = Water.get_standard()

"""
#Geometrical parameters
        center = [0, 0, 0]
        model = model.lower()
        if model == 'tip3p':
            r_oh = 0.95720
            a_hoh =  104.52
        elif model == 'spc':
            r_oh = 1.00
            a_hoh =  109.47
        r_oh = r_oh / a0
        d = (90 - a_hoh/2 ) * np.pi / 180
        origin = np.array( [ 0, 0, 0] )

        h1 = Atom( AA = AA, element = "H" )
        h2 = Atom( AA = AA, element = "H" )
        o =  Atom( AA = AA, element = "O" )

        o.x = center[0]
        o.y = center[1]
        o.z = center[2] 

        h1.x = (center[0] + r_oh * np.cos(d))
        h1.y = center[1] 
        h1.z = (center[2] + r_oh* np.sin(d))

        h2.x = (center[0] - r_oh * np.cos(d)) 
        h2.y = center[1] 
        h2.z = (center[2] + r_oh* np.sin(d))
        o.order = 1
        h1.order = 2
        h2.order = 3
        w = Water( AA = AA)
        w.append( o )
        w.append( h1 )
        w.append( h2 )
        w.o.pdb_name = 'OW'
        w.h1.pdb_name = 'HW1'
        w.h2.pdb_name = 'HW2'
        if worst:
            w.populate_bonds()
            w.populate_angles()
            w.h1.scale_angle( 0.988 )
            w.h1.scale_bond( 0.985 )
            w.h2.scale_bond( 1.015 )
            w.inv_rotate()
        return w

    def is_worst(self):
        self.populate_bonds()
        self.populate_angles()
        r1 = self.h1.dist_to_atom( self.o )
        r2 = self.h2.dist_to_atom( self.o )
        a = self.h1.get_angle( self.o, self.h2 ) /np.pi * 180
        if np.allclose( [r1, r2, a], [1.7817, 1.8360, 103.2658 ] ,atol = 0.0001) or np.allclose( [r2, r1, a], [1.7817, 1.8360, 103.2658 ] ,atol = 0.0001) :
            return True
        return False

    def translate_o(self, r):
        vec = r - self.o.r
        for at in self:
            at.x = vec[0] + at.x 
            at.y = vec[1] + at.y 
            at.z = vec[2] + at.z 
        return self
#Center of oxygen
    @property
    def coo(self):
        return self.o.r
    def append(self, atom):
        """
Override list append method, will add up to 3 atoms,
1 must be oxygen, 2 must be hydrogens.

.. code:: python

    >>> m = Water()
    >>> m.append( Atom( z = 0.11, element = 'H' ) )
    >>> m.coc

"""
        if len(self) > 3:
            print "tried to add additional atoms to water, exiting"
            raise SystemExit

        if not isinstance( atom, Atom ):
            print "wrong class passed to water append"
            raise SystemExit

        if atom.element == "H":
            if self.no_hydrogens:
                self.h1 = atom
                atom.Molecule = self
                self.no_hydrogens = False
            else:
                self.h2 = atom
                atom.Molecule = self
        if atom.element == "O":
            self.o = atom
            atom.Molecule = self
            atom.name = "O1"
#Add the atom
        super( Water, self).append(atom)

#Define water center, by default set it to center of nuclei
        if (self.h1 and self.h2 and self.o):
            pass

        if self.res_id:
            if self.res_id != atom.res_id:
                print "Tried to add %s to %s, exiting" %(atom, self)
                raise SystemExit
        else:
#Initialize water res_id from atomic res_id
            self._res_id = atom.res_id

        if len(self) == 3:
            hyd1, hyd2 = [i for i in self if i.element == "H" ]
            d1 = hyd1.dist_to_point( [1,1,1] )
            d2 = hyd2.dist_to_point( [1,1,1] )
            if d1 < d2:
                self.h1 = hyd1
                self.h2 = hyd2
                hyd1.name = "H2"
                hyd2.name = "H3"
            else:
                self.h1 = hyd2
                self.h2 = hyd1
                hyd1.name = "H2"
                hyd2.name = "H3"

#Water string method
    def __str__(self):
        return "WAT" + str(self.res_id) 
    def dist_to_point( self , point ):
        return np.sqrt(np.sum((self.coo - np.array(point))**2))

    def dist_to_water(self, other):
        return np.sqrt(np.sum((self.coo - other.coo)**2) )

    def get_xyz_string(self, ):
        st = "%d\n\n" % len(self)
        for i in self:
            st += "{0:10s} {1:10f} {2:10f} {3:10f}\n".format(\
                    i.element, i.x,  i.y , i.z )
        return st

    @staticmethod
    def read_waters( fname , in_AA = True, out_AA = True , N_waters = 1):
        """From file with name fname, return a list of all waters encountered"""
        atoms = []
        if fname.endswith( ".xyz" ) or fname.endswith(".mol"):
            pat_xyz = re.compile(r'^\s*(\w+)\s+(-*\d*.+\d+)\s+(-*\d*.+\d+)\s+(-*\d*.+\d+) *$')
            for i in open( fname ).readlines():
                if pat_xyz.match(i):
                    f = pat_xyz.match(i).groups()
                    matched = pat_xyz.match(i).groups()
                    kwargs = { "element" :  matched[0], "x" : matched[1],
                            "y" : matched[2], "z" : matched[3] }
                    tmpAtom = Atom( **kwargs )
                    atoms.append( tmpAtom )
        elif fname.endswith( ".pdb" ):
            pat1 = re.compile(r'^(ATOM|HETATM)')
            for i in open( fname ).readlines():
                if pat1.search(i):
                    if ( i[11:16].strip() == "SW") or (i[11:16] == "DW"):
                        continue
                    kwargs = {
                            "AA" : in_AA,
                            "x" : float(i[30:38].strip()),
                            "y" : float(i[38:46].strip()),
                            "z" : float(i[46:54].strip()),
                            "element": i[11:16].strip()[0] }
                    tmpAtom = Atom( **kwargs )
                    atoms.append( tmpAtom )
        elif fname.endswith( ".out" ):
            pat_xyz = re.compile(r'^(\w+)\s+(-*\d*.+\d+)\s+(-*\d*.+\d+)\s+(-*\d*.+\d+) *$')
            for i in open( fname ).readlines():
                if pat_xyz.match(i):
                    f = pat_xyz.match(i).groups()
                    tmpAtom = Atom(f[0][0], float(f[1]), float(f[2]), float(f[3]), 0)
                    atoms.append( tmpAtom )
#loop over oxygen and hydrogen and if they are closer than 1 A add them to a water
        waters = []
        cnt = 1
        if fname.endswith( ".xyz" ) or fname.endswith(".mol"):
            for i in atoms:
                if i.element == "H":
                    continue
                if i.in_water:
                    continue
                tmp = Water()
                i.in_water = True
                tmp.append( i )
                for j in atoms:
                    if j.element == "O":
                        continue
                    if j.in_water:
                        continue
#If in angstrom
                    if in_AA:
                        if i.dist_to_atom(j) < 1.1:
                            tmp.append ( j )
                            j.in_water = True
                    else:
                        if i.dist_to_atom(j) < 1.1/a0:
                            tmp.append ( j )
                            j.in_water = True
                tmp.res_id = cnt
                cnt += 1
                waters.append( tmp )
        elif fname.endswith( ".pdb" ):
#Find out the size of the box encompassing all atoms
            xmin = 10000.0; ymin = 10000.0; zmin = 10000.0; 
            xmax = -10000.0; ymax = -10000.0; zmax = -10000.0; 
            for i in atoms:
                if i.x < xmin:
                    xmin = i.x
                if i.y < ymin:
                    ymin = i.y
                if i.z < zmin:
                    zmin = i.z
                if i.x > xmax:
                    xmax = i.x
                if i.y > ymax:
                    ymax = i.y
                if i.z > zmax:
                    zmax = i.z
            center = np.array([ xmax - xmin, ymax -ymin, zmax- zmin]) /2.0
            wlist = []
            for i in atoms:
                if i.element == "H":
                    continue
                if i.in_water:
                    continue
                tmp = Water()
#__Water__.append() method will update the waters residue number and center coordinate
#When all atoms are there
#Right now NOT center-of-mass
                i.in_water= True
                tmp.append(i)
                for j in atoms:
                    if j.element == "O":
                        continue
                    if j.in_water:
                        continue
#1.05 because sometimes spc water lengths can be over 1.01
                    if in_AA:
                        if i.dist_to_atom(j) <= 1.05:
                            j.in_water = True
                            tmp.append( j )
                    else:
                        if i.dist_to_atom(j) <= 1.05/a0:
                            j.in_water = True
                            tmp.append( j )
                tmp.res_id = cnt
                cnt += 1
                wlist.append( tmp )
            wlist.sort( key = lambda x: x.dist_to_point( center ))
            center_water = wlist[0]
            cent_wlist = wlist[1:]
            cent_wlist.sort( key= lambda x: x.dist_to_water( center_water) )

            if N_waters< 1:
                print "Please choose at least one water molecule"
                raise SystemExit
            waters = [center_water] + cent_wlist[ 0 : N_waters - 1 ]

        elif fname.endswith( ".out" ):
            for i in atoms:
                if i.element == "H":
                    continue
                if i.in_water:
                    continue
                tmp = Water()
                i.in_water = True
                tmp.append( i )
                for j in atoms:
                    if j.element == "O":
                        continue
                    if j.in_water:
                        continue
#If in cartesian:
                    if i.AA:
                        if i.dist(j) < 1.0:
                            tmp.append ( j )
                            j.in_water = True
                    else:
                        if i.dist(j) < 1.0/a0:
                            tmp.append ( j )
                            j.in_water = True
                tmp.res_id = cnt
                cnt += 1
                waters.append( tmp )
        for wat in waters:
            for atom in wat:
                atom._res_id = wat.res_id
        if in_AA:
            if not out_AA:
                for wat in waters:
                    wat.to_AU()
        return waters
     
    @staticmethod
    def get_string_from_waters( waters, max_l = 1, pol = 2 , hyper = 0, dist = False, AA = False ):
        """ Converts list of waters into Olav string for hyperpolarizable .pot"""
# If the properties are in distributed form, I. E. starts from Oxygen, then H in +x and H -x

        if AA:
            str_ = "AA"
        else:
            str_ = "AU"
        string = "%s\n%d %d %d %d\n" % ( str_, len(waters)*3,
                max_l, pol, hyper )
        for i in waters:
            for at in i:
                string +=  " ".join([str(at.res_id)] + map(str,at.r)) + " "
                string += at.Property.potline( max_l=max_l, pol=pol, hyper= hyper)
                string += '\n'
        return string




class Methanol(Molecule):
    """
Not yet implemented, only needs get_euler and z-matrix to be specific.
    """

    def __init__(self, *args, **kwargs):
        super( Methanol, self).__init__(**kwargs)

    def append(self, atom):
        """Typical append for each seperate molecule class"""
        if len(self) == 6:
            print "tried to add additional atoms to methanol, exiting"
            raise SystemExit

        if not isinstance( atom, Atom ):
            print "wront class passed to methanol append"
            raise SystemExit

        if atom.element == "H":
            if self.no_hydrogens:
                self.h1 = atom
                self.no_hydrogens = False
            else:
                self.h2 = atom

        if atom.element == "O":
            self.o = atom
#Add the atom
        super( Molecule, self).append(atom)

#Define methanol center, by default set it to center of C=O bond

        if len(self) == 6:
            self.center = sum([ i.r for i in self ] ) / len(self)
            #hc = charge_dict[ self.h1.element ]
            #oc = charge_dict[ self.h1.element ]
            #self.coc = np.array([ self.h1.x * hc  + self.h2.x *hc + self.o.x *oc,  \
            #    self.h1.y *hc + self.h2.y *hc + self.o.y *oc , \
            #    self.h1.z *hc + self.h2.z *hc + self.o.z *oc ]) /( 2*hc +oc)
        if self.res_id:
            if self.res_id != atom.res_id:
                print "Tried to add %s to %s, exiting" %(atom, self)
                raise SystemExit
        else:
#Initialize res_id from atomic res_id
            self._res_id = atom.res_id
#Also calculate center now
        if len(self) == 6:
            h1, h2, h3, h4 = [i for i in self if i.element == "H" ]
            #print "All hyds added for methanol %s" %str(self)
            #raise SystemExit
            #d1 = hyd1.dist_to_point( [1,1,1] )
            #d2 = hyd2.dist_to_point( [1,1,1] )
            #if d1 < d2:
            #    self.h1 = hyd1
            #    self.h2 = hyd2
            #else:
            #    self.h1 = hyd2
            #    self.h2 = hyd1


class Cluster(list):
    """
**Molecule container which groups molecules into quantum mechanics, molecular mechanics, and qm/mm regions for easy in generating input files for QM/MM.**
"""
    def __init__(self, *args, **kwargs):
        """ Can be a list of molecules or a list of clusters"""

        super(Cluster, self).__init__()
        self._chain_id = None
        self.Property = None
        self.atom_list = []
        if type(args) == tuple:
            if len(args) == 1:
                if type(args[0]) == list:
                    for i in args[0]:
                        self.add( i )
                else:
                    self.add( args[0] )
            else:
                for item in args:
                    self.add( item )

    def connect_everything(self):
        self.connect_atoms_to_cluster()
        self.connect_atoms_to_molecules()
        self.connect_molecules_to_cluster()

    def connect_atoms_to_molecules(self):
        for res in [mol for mol in self if isinstance(mol,Molecule) ]:
            for at in res:
                at.Molecule = res
    def connect_atoms_to_cluster(self):
        for cluster in self:
            for res in [res for res in cluster if isinstance(res, Molecule)]:
                for at in res:
                    at.Cluster = cluster

    def connect_molecules_to_cluster(self):
        for res in [mol for mol in self if isinstance(mol,Molecule) ]:
            res.Cluster = self

    def get_pdb_string(self):
        st = ""
        for at in [at for mol in self for at in mol]:
            st += at.pdb_string()
        return st

    @property
    def density(self):
        """Return the density in SI units kg/m^3"""
        N_A = 6.02214129e+23

#g/mol
        M = sum( [at.mass for mol in self for at in mol] )
#AU**3 or AA**3
        v = utilz.convex_hull_volume( np.array([at.r for mol in self for at in mol]))
        if not self.AA:
            v /= a0**3
#g mol**-1 -> g
        m = M / N_A
#AA**-3 -> m**-3
        v *= 1e-30
#g -> kg
        m *= 1e-3
#kg m**-3 = g cm**-3
        return m/v

    def g_list_from_damped(self, 
            max_l = 1,
            pol = 22,
            hyp = 2,
            rq = 1e-9,
            rp = 1e-9,
            AA_cutoff = 1.5,
            nullify = False,
            model = 'Thole'
            ):
        """Given cutoff in Angstromgs, will return a GassuanQuadrupoleList
        where atomic within AA_cutoff between different interacting segments
        
        has a damped gaussian """
        from pd.gaussian import GaussianQuadrupoleList
        from pd.thole import TholeList

        opts = { 'Thole' : TholeList, 'Gaussian' :GaussianQuadrupoleList,
                't' : TholeList, 'g' :GaussianQuadrupoleList,
                }

        aa = self.AA
        self.to_AU()
        g = opts[model].from_string( self.get_qmmm_pot_string() )
        for atom, res in map( lambda x: [x, x.residue], self.min_dist_atoms_seperate_res(AA_cutoff) ):
            ind = reduce( lambda a, x: a + len(x), res.Chain[:res.order_nr],0)+atom.order_nr
            g[ ind ]._R_q = rq
            g[ ind ]._R_p = rp
            if nullify:
                g[ ind ]._q = 0.0
                g[ ind ]._p0 = np.zeros( (3,) )
                g[ ind ]._a0 = np.zeros( (3,3,) )
                g[ ind ]._Q0 = np.zeros( (3,3,) )
                g[ ind ]._b0 = np.zeros( (3,3,3,) )
        if aa:
            self.to_AA()
        return g

# Slicing the cluster givees back a cluster, but only accessing one index gives molecule
    def __add__(self, other):
        return Cluster(list.__add__(self, other))
    def __getslice__(self, i, j):
        return self.__getitem__(slice(i, j))
    def __getitem__(self, item):
        if isinstance( item, slice ):
            result = list.__getitem__(self, item)
            try:
                return Cluster(result)
            except TypeError:
                return result
        else:
            return super(Cluster,self).__getitem__( item )

    @property
    def AA(self):
        AA = [at for res in self for at in res ][0].AA
        for each in [at for res in self for at in res ]:
            try:
                assert each.AA == AA
            except AssertionError:
                logging.error("All objects in cluster are not of same unit")
        return AA

#Cluster string method
    def __str__(self):
        return " ".join( [ str(i) for i in self ] )
    
    def save(self, fname = "cluster.p"):
        pickle.dump( self, open(fname, 'wb' ), protocol = 2 )

    @staticmethod
    def load(fname = 'cluster.p'):
        if not os.path.isfile( fname):
            raise IOError
        return pickle.load( open(fname, 'rb' ) )
    
    def get_dalton_qmmm(self, max_l = 2, pol =2, hyp = 0, qmmm_type = 'peqm',
        basis = ("ano-1 2 1", "ano-1 3 2 1", "ano-2 5 4 1" ) ):
        """Generate DALTON.INP, MOLECULE.INP and POTENTIAL.INP for cluster"""
        dal, mol, pot = ['' for i in range(3)]
        if qmmm_type == 'peqm':
            dal = Generator.get_pe_b3lyp_dal(max_l = max_l, AA = self.AA)
            pot = self.get_pe_pot_string(max_l = max_l,
                    pol = pol,
                    hyp = hyp)
        elif qmmm_type == 'qmmm':
            dal = Generator.get_qmmm_b3lyp_dal()
            pot = self.get_qmmm_pot_string(max_l = max_l,
                    pol = pol,
                    hyp = hyp,
                    ignore_qmmm = False)
        mol = self.get_qm_mol_string( basis = basis )
        return dal, mol, pot


    
    @staticmethod
    def get_all_molecules_from_file(fil,
            in_AA = False,
            out_AA = False,
            ):
#dont add atoms to molecule if they are within 1.5 angstrom
        max_dist = 1.5
        if in_AA :
            max_dist / a0

        """Given pdb/mol/xyz  file return a Cluster with all seperate molecules"""
        if fil.endswith('.xyz'):
            with open(fil,'r') as f_:
                pass
    def min_dist_atoms_seperate_res(self, AA_cutoff = 1.5 ):
        """Return list of atoms which have an other atom closer than 1.5 AA to them
        and are not in the same residue
        
        """
        tmp = []
        ats = self.min_dist_atoms( AA_cutoff = AA_cutoff )
        for i, at in enumerate( ats ):
            for j in range( i, len( ats ) ):
                if ats[i].res_id == ats[j].res_id:
                    continue
                if ats[i].dist_to_atom( ats[j] ) < AA_cutoff:
                    tmp.append( ats[i] )
                    tmp.append( ats[j] )
        return read_dal.unique( tmp )



     
    def min_dist_atoms(self, AA_cutoff = 1.5):
        """Return list of atoms which have an other atom closer than 1.5 AA to them
        
.. code:: python
    
    >>> c = Cluster()
    >>> c.add( Water.get_standard())
    >>> for at in c.min_dist_atoms():
            print at
        O1 0.000000 0.000000 0.000000
        H2 1.430429 0.000000 1.107157
        H3 -1.430429 0.000000 1.107157
        """
        if not self.AA:
            AA_cutoff /= a0
        N_ats = reduce( lambda a,x: a + len(x) , [res for res in self], 0 )
        d_mat = np.full( (N_ats, N_ats ), np.inf )

        ats = [at for res in self for at in res]
        for i1, at1 in enumerate( ats ):
            for i2, at2 in enumerate( ats ):
                if at1 == at2:
                    continue
                d_mat [i1, i2] = at1.dist_to_atom( at2 )
        x, y = np.where( d_mat < AA_cutoff )[0], np.where( d_mat < AA_cutoff) [1]
        min_ats = []

        for xi, zi in zip( x, y ):
            min_ats.append( ats[xi] )
            min_ats.append( ats[zi] )

        return read_dal.unique(min_ats)


    def min_dist_coo(self):
        dist = np.zeros( (len(self),len(self)) )
        new = np.zeros( len(self) - 1 )

        for i in range(len(self)):
            for j in range( i, len(self)):
                if i == j:
                    continue
                dist[i,j] = np.linalg.norm(self[i].coo - self[j].coo)
        for i in range( len(dist) - 1 ):
            c = dist[i,i+1:].copy()
            c.sort()
            new[i] = c[0]
        new.sort()
        return new

    def min_dist_com(self):
        dist = np.zeros( len(self) )
        for i in range(len(self)):
            for j in range(i ,len(self)):
                if i == j:
                    continue
                dist[i] = ( np.linalg.norm(self[i].com - self[j].com) )
        dist.sort()
        return dist

    def __eq__(self, other):
        """docstring for __eq__ """
        if not len(self) == len(other):
            return False
        for i, (m1, m2) in enumerate( zip(self, other) ):
            if m1 != m2:
                return False
        return True
        
    def get_qm_xyz_string(self, AA = False):
# If basis len is more than one, treat it like molecular ano type
# Where first element is for first row of elements

        st = "%d\n\n" % sum( [ len(m) for m in self if m.in_qm ] )
        for i in [all_el for mol in self for all_el in mol if mol.in_qm]:
            st += "{0:5s}{1:10.5f}{2:10.5f}{3:10.5f}\n".format( i.element, i.x, i.y, i.z )
        return st
# Molecules iterator
    @property
    def molecules(self):
        for mol in [mol for mol in self if isinstance(mol,Molecule)]:
            yield mol
# Atoms iterator
    @property
    def atoms(self):
        for atom in [at for mol in self for at in mol]:
            yield atom
   
# Specifi
    @property
    def coc(self):
    #obj should be atom
        return sum( [at.r * charge_dict[at.element] for mol in self for at in mol])\
                /sum( map(float,[charge_dict[at.element] for mol in self for at in mol]) )


    def plot(self, copy = True, center = False ):
        """
Plot Cluster a 3D frame in the cluster

.. code:: python

    >>> m = Cluster()
    >>> m.add_atom( Atom(element = 'H', x = 1, z = 1) )
    >>> m.add_atom( Atom(element = 'H', x =-1, z = 1) )
    >>> m.add_atom( Atom(element = 'O', z = 0) )
    >>> m.plot()
    
"""

#Make a copy in order to not change original, and perform plot on it
        if copy:
            copy = copymod.deepcopy( self )
            copy.connect_everything()
        else:
            copy = self
            copy.connect_everything()
        if center:
            copy.translate( -self.com )

        for mol in [mol for mol in copy if isinstance( mol, Molecule) ]:
            mol.populate_bonds()

#Plot in nice xyz axis
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d' )

#Plot bonds
        for mol in [mol for mol in copy if isinstance( mol, Molecule) ]:
            for atom in mol:
                for key in mol.bond_dict[ atom ]:
                    ax.plot( [key.x, atom.x],
                             [key.y, atom.y],
                             [key.z, atom.z], color = 'black' )



        ax.plot( [0, 1, 0, 0, 0, 0], [0,0 ,0,1,0,0], [0,0,0,0,0,1] )
        ax.text( 1.1, 0, 0, "X", color = 'red' )
        ax.text( 0, 1.1, 0, "Y", color = 'red' )
        ax.text( 0, 0, 1.1, "Z", color = 'red' )

        for i in copy:
            for j in i:
                ax.plot( [j.x], [j.y], [j.z], j.Molecule.style[j.element], linewidth= j.Molecule.linewidth[j.element] )
        ax.set_zlim3d( -5,5)
        plt.xlim(-5,5)
        plt.ylim(-5,5)
        plt.show()



    def get_qm_mol_string(self, basis = ("ano-1 2 1", "ano-1 3 2 1", "ano-2 5 4 1" ) , AA = False):
# If basis len is more than one, treat it like molecular ano type
# Where first element is for first row of elements

        if len( basis ) > 1:
            # Set row index number to periodic table one
            el_to_rowind = {"H" : 0, "C" : 1, "O" : 1, "N" : 1 , "S" : 2  }
        else:
            # Keep all 0, since basis is only len 1
            el_to_rowind = {"H" : 0, "C" : 0, "O" : 0, "N" : 0, "S" : 0 }

        st = ""
        comm1 = "QM: " + " ".join( [ str(m) for m in self if m.in_qm] )[:72]
        comm2 = "MM: " + " ".join( [ str(m) for m in self if m.in_mm] )[:73]
        uni = utilz.unique([ at.element for mol in self for at in mol if mol.in_qm])
        s_ = ""
        if AA: s_ += "Angstrom"

        st += "ATOMBASIS\n%s\n%s\nAtomtypes=%d Charge=0 Nosymm %s\n" %( \
                comm1,
                comm2,
                len(uni),
                s_)
        for el in uni:
            st += "Charge=%s Atoms=%d Basis=%s\n" %( str(charge_dict[el]),
                    len( [all_el for mol in self for all_el in mol if ((all_el.element == el) and mol.in_qm )] ),
                     basis[ el_to_rowind[el] ] )
            for i in [all_el for mol in self for all_el in mol if ((all_el.element == el) and mol.in_qm) ]:
                st += "{0:5s}{1:10.5f}{2:10.5f}{3:10.5f}\n".format( i.element, i.x, i.y, i.z )
        return st


# Specific output for PEQM calculation in dalton, all molecules exclude itself
    def get_pe_pot_string( self, max_l = 0, pol = 1, hyp = 0, out_AA = False ):
        max_len = max([len(mol) for mol in self if mol.in_mm])
        self.order_mm_atoms()
        st = r'!%s' % (self ) + '\n'
        st += r'@COORDINATES' + '\n'
        st += '%d\n' % sum([len(i) for i in self if i.in_mm ])
        if self.AA:
            st += "AA\n"
        else:
            st += "AU\n"
        #st += '%d\n' %len(mol)
        for mol in [m for m in self if m.in_mm]:
            for at in mol:
                st += "%s %.5f %.5f %.5f\n" % (at.element, \
                        at.x, at.y, at.z )

        st += r'@MULTIPOLES'  + '\n'
        st += 'ORDER 0\n'
        st += '%d\n' % sum([len(i) for i in self if i.in_mm ])
        for mol in [m for m in self if m.in_mm]:
            for at in mol:
                st += "%d %.5f\n" % (at.cluster_order, at.p.q )

        if max_l >= 1:
            st += 'ORDER 1\n'
            st += '%d\n' % sum([len(i) for i in self if i.in_mm ])
            for mol in [m for m in self if m.in_mm]:
                for at in mol:
                    st += ("{0:1d} %s\n" %("".join(["{%d:.5f} "%i for i in range(1,4)]))).format( *tuple([at.cluster_order] + at.Property["dipole"].tolist() ))
        if max_l >= 2:
            st += 'ORDER 2\n'
            st += '%d\n' % sum([len(i) for i in self if i.in_mm ])
            for mol in [m for m in self if m.in_mm]:
                for at in mol:
                    st += ("{0:1d} %s\n" %("".join(["{%d:.5f} "%i for i in range(1,7)]))).format( *tuple([at.cluster_order] + at.Property["quadrupole"].tolist() ))

        if pol > 1:
            st += r'@POLARIZABILITIES' + '\n'
            st += 'ORDER 1 1\n'
            st += '%d\n' % sum([len(i) for i in self if i.in_mm ])
            if pol % 2 == 0:
                for mol in [m for m in self if m.in_mm]:
                    for at in mol:
                        st += ("{0:1d} %s\n" %("".join(["{%d:.5f} "%i for i in range(1,7)]))).format( *tuple([at.cluster_order] + at.Property["alpha"].tolist() ))

        st += 'EXCLISTS\n%d %d\n' %( sum([len(i) for i in self if i.in_mm ])
     , len(mol))
        for mol in [m for m in self if m.in_mm]:
            for each in mol.exclists( max_len = max_len ):
                st += ("%s\n"%("".join(["{%d:1d} "%i for i in range(max_len)]))).format( *each )
        return st
# This is the old *QMMM input style in dalton, also valid for PointDipoleList
    def get_qmmm_pot_string( self, max_l = 1,
            pol = 22,
            hyp = 1,
# If complicated molecule, set dummy_pd to a coordinate to place the net property
            dummy_pd = False,
#Set ignore_qmmm to false to only write qmmm .pot file for molecues in mm region
            ignore_qmmm = True ):

#make sure that each unique residue is in seperate residue in POT output
        self.order_mm_mols()

# We need to check that it is not in LoProp mode
        if dummy_pd:
            assert self.LoProp == False

        if self.AA:
            st = "AA\n"
        else:
            st = "AU\n"
# Old qmmm format requires integer at end to properly read charges
        if hyp == 0:
            hyp_int = 1
        if ignore_qmmm:
            st += "%d %d %d %d\n" % (sum([len(i) for i in self ]), 
                    max_l, pol, 1 )
            st += "".join( [at.potline(max_l, pol, hyp) for mol in self for at in mol ] )
        else:
            st += "%d %d %d %d\n" % (sum([len(i) for i in self if i.in_mm ]), 
                    max_l, pol, 1 )
            st += "".join( [at.potline(max_l, pol, hyp) for mol in self for at in mol if mol.in_mm] )
        return st

    def get_xyz_string_qmmm(self, both= False, qm_region = False, mm_region = False ):
        ats = []
        if qm_region:
            st = "%d\n\n" % sum([len(i) for i in self if i.in_qm ])
            ats = [at for mol in self for at in mol if mol.in_qm]
        if mm_region:
            st = "%d\n\n" % sum([len(i) for i in self if i.in_mm ])
            ats = [at for mol in self for at in mol if mol.in_mm]
        if qm_region and mm_region:
            st = "%d\n\n" % sum([len(i) for i in self if i.in_qmmm ])
            ats = [at for mol in self for at in mol if mol.in_qmmm]
        if both:
            ats = [at for mol in self for at in mol ]
        st = "%d\n\n" %len(ats)
        for i in ats:
            st += "{0:10s} {1:10f} {2:10f} {3:10f}\n".format(\
                    i.element, i.x,  i.y , i.z )
        return st

    def get_xyz_string(self, ):
        st = "%d\n\n" % sum([len(i) for i in self ])
        for mol in self:
            for i in mol:
                st += "{0:10s} {1:10f} {2:10f} {3:10f}\n".format(\
                        i.element, i.x,  i.y , i.z )
        return st
    def order_mm_atoms(self):
        cnt = 1
        for a in [at for m in self for at in m if m.in_mm]:
            a.cluster_order = cnt
            cnt += 1
    def order_mm_mols(self):
        cnt = 1
        for mol in self:
            mol.cluster_order = cnt
            cnt += 1

    def update_water_props(self, model = "TIP3P",
            method = "HF", basis = "ANOPVDZ", dist = True,
            freq = "0.0"):
        from template import Template

        kwargs_dict = Template().get( *(model, method, basis,
            dist , freq ))
        for wat in self:
            t1, t2, t3 = wat.get_euler()
            for at in wat:
                Property.add_prop_from_template( at, kwargs_dict )
                at.Property.transform_ut_properties( t1, t2, t3)


    @staticmethod
    def get_water_cluster( fname , in_AA = False, out_AA = False , N_waters = 1000 ):
        """
Return a cluster of water molecules given file.

.. code:: python

    >>> c = Cluster.get_water_cluster( 'somefile.mol' , in_AA = False, out_AA = False, N_waters = 10 )
    >>> print len( c )
    10

"""
        atoms = []
        c = Cluster()
        if fname.endswith( ".xyz" ) or fname.endswith(".mol"):
            pat_xyz = re.compile(r'^\s*(\w+)\s+(-*\d*.+\d+)\s+(-*\d*.+\d+)\s+(-*\d*.+\d+) *$')
            for i in open( fname ).readlines():
                if pat_xyz.match(i):
                    f = pat_xyz.match(i).groups()
                    matched = pat_xyz.match(i).groups()
                    kwargs = { "element" :  matched[0], "x" : matched[1],
                            "y" : matched[2], "z" : matched[3] }
                    tmpAtom = Atom( **kwargs )
                    atoms.append( tmpAtom )

        elif fname.endswith( ".pdb" ):
            pat1 = re.compile(r'^(ATOM|HETATM)')
#Temporary atom numbering so that it is compatible with PEQM reader in dalton
            for i in open( fname ).readlines():
                if pat1.search(i):
                    if ( i[11:16].strip() == "SW") or (i[11:16] == "DW") \
                            or (i[11:16].strip() == "MW"):
                        continue
                    kwargs = {
                            "AA" : in_AA,
                            "x" : float(i[30:38].strip()),
                            "y" : float(i[38:46].strip()),
                            "z" : float(i[46:54].strip()),
                            "element": i[11:16].strip()[0],
                            "number" : i[6:11].strip()  }
                    tmpAtom = Atom( **kwargs )
                    atoms.append( tmpAtom )
        elif fname.endswith( ".out" ):
            pat_xyz = re.compile(r'^(\w+)\s+(-*\d*.+\d+)\s+(-*\d*.+\d+)\s+(-*\d*.+\d+) *$')
            for i in open( fname ).readlines():
                if pat_xyz.match(i):
                    f = pat_xyz.match(i).groups()
                    tmpAtom = Atom(f[0][0], float(f[1]), float(f[2]), float(f[3]), 0)
                    atoms.append( tmpAtom )
        elif fname.endswith( '.log' ):
            pat_atoms = re.compile ( r'NAtoms=\s+(\d+)' )
            pat_xyz = re.compile ( r'Standard ori' )
            lines = open(fname).readlines()
            conf = []
            for line in lines:
                if pat_atoms.search( line ):
                    N = int(pat_atoms.search( line ).group(1))
                    break
            for ind, line in enumerate(lines):
                if pat_xyz.search( line ):
                    conf.append( "\n".join( map( lambda x:x.strip('\n'), lines[ ind +5 : ind + 5 + N ]) ) )
            for each in conf[-1].split('\n'):
                spl = each.split()
                tmpAtom = Atom( element = elem_array[ int(spl[1]) ],
                        AA = True,
                        x = float(spl[3]),
                        y = float(spl[4]),
                        z = float(spl[5]),)
                atoms.append( tmpAtom )

#loop over oxygen and hydrogen and if they are closer than 1 A add them to a water
        waters = []
        cnt = 1
        if fname.endswith(".log") or fname.endswith( ".xyz" ) or fname.endswith(".mol"):
            xmin = 10000.0; ymin = 10000.0; zmin = 10000.0; 
            xmax = -10000.0; ymax = -10000.0; zmax = -10000.0; 
            for i in atoms:
                if i.x < xmin:
                    xmin = i.x
                if i.y < ymin:
                    ymin = i.y
                if i.z < zmin:
                    zmin = i.z
                if i.x > xmax:
                    xmax = i.x
                if i.y > ymax:
                    ymax = i.y
                if i.z > zmax:
                    zmax = i.z
            center = np.array([ xmax - xmin, ymax -ymin, zmax- zmin]) /2.0
            for i in atoms:
                if i.element == "H":
                    continue
                if i.in_water:
                    continue
                tmp = Water( AA = in_AA)
#Gaussian output seems to have Angstrom always
                if fname.endswith( '.log' ):
                    tmp.AA = True
                i.in_water = True
                tmp.append( i )
                for j in atoms:
                    if j.element == "O":
                        continue
                    if j.in_water:
                        continue
#If in angstrom
                    if in_AA:
                        if i.dist_to_atom(j) < 1.1:
                            tmp.append ( j )
                            j.in_water = True
                    else:
                        if i.dist_to_atom(j) < 1.1/a0:
                            tmp.append ( j )
                            j.in_water = True
                tmp._res_id = cnt
                cnt += 1
                waters.append( tmp )
            waters.sort( key = lambda x: x.dist_to_point( center ))
            center_water = waters[0]
            cent_wlist = waters[1:]
            cent_wlist.sort( key= lambda x: x.dist_to_water( center_water) )

            if N_waters < 1:
                print "WARNING ; chose too few waters in Cluster.get_water_cluster"
                raise SystemExit
# Ensure that cluster has ordered water structure from first index
            waters = [center_water] + cent_wlist[ 0 : N_waters - 1 ]
            for i in waters:
                c.append(i)
        elif fname.endswith( ".pdb" ):
#Find out the size of the box encompassing all atoms
            xmin = 10000.0; ymin = 10000.0; zmin = 10000.0; 
            xmax = -10000.0; ymax = -10000.0; zmax = -10000.0; 
            for i in atoms:
                if i.x < xmin:
                    xmin = i.x
                if i.y < ymin:
                    ymin = i.y
                if i.z < zmin:
                    zmin = i.z
                if i.x > xmax:
                    xmax = i.x
                if i.y > ymax:
                    ymax = i.y
                if i.z > zmax:
                    zmax = i.z
            center = np.array([ xmax - xmin, ymax -ymin, zmax- zmin]) /2.0
            wlist = []
            for i in atoms:
                if i.element == "H":
                    continue
                if i.in_water:
                    continue
                tmp = Water( AA = in_AA )
#__Water__.append() method will update the waters residue number and center coordinate
#When all atoms are there
#Right now NOT center-of-mass
                i.in_water= True
                tmp.append(i)
                for j in atoms:
                    if j.element == "O":
                        continue
                    if j.in_water:
                        continue
#1.05 because sometimes spc water lengths can be over 1.01
                    if in_AA:
                        if i.dist_to_atom(j) <= 1.05:
                            j.in_water = True
                            tmp.append( j )
                    else:
                        if i.dist_to_atom(j) <= 1.05/a0:
                            j.in_water = True
                            tmp.append( j )
                tmp._res_id = cnt
                cnt += 1
                wlist.append( tmp )

            wlist.sort( key = lambda x: x.dist_to_point( center ))
            center_water = wlist[0]
            cent_wlist = wlist[1:]
            cent_wlist.sort( key= lambda x: x.dist_to_water( center_water) )


            if N_waters < 1:
                print "WARNING ; chose too few waters in Cluster.get_water_cluster"
                raise SystemExit

# Ensure that cluster has ordered water structure from first index
            waters = [center_water] + cent_wlist[ 0 : N_waters - 1 ]
            for i in waters:
                c.append(i)
        elif fname.endswith( ".out" ):
            for i in atoms:
                if i.element == "H":
                    continue
                if i.in_water:
                    continue
                tmp = Water()
                i.in_water = True
                tmp.append( i )
                for j in atoms:
                    if j.element == "O":
                        continue
                    if j.in_water:
                        continue
#If in cartesian:
                    if i.AA:
                        if i.dist_to_atom(j) < 1.0:
                            tmp.append ( j )
                            j.in_water = True
                    else:
                        if i.dist_to_atom(j) < 1.0/a0:
                            tmp.append ( j )
                            j.in_water = True
                tmp.res_id = cnt
                cnt += 1
                waters.append( tmp )
        for wat in c:
            for atom in wat:
                atom._res_id = wat.res_id

        if in_AA or fname.endswith( '.log' ):
            if not out_AA:
                for wat in c:
                    wat.to_AU()
        if not in_AA:
            if out_AA:
                for wat in c:
                    wat.to_AA()
        for wat in c:
            wat.o.order = 1
            wat.h1.order = 2
            wat.h2.order = 3
        c.set_qm_mm(100)
        return c


    def mol_too_close(self, mol, dist = 2.5):
        for mols in self:
            for ats in mols:
                for at in mol:
                    if at.dist_to_atom( ats ) < dist:
                        return True
        return False

    def attach_properties(self, 
            model = "TIP3P",
            method = "HF",
            basis = "ANOPVDZ",
            loprop = True,
            freq = "0.0"):
        """
Attach property to all atoms and oxygens, by default TIP3P/HF/ANOPVDZ, static
        """
        templ = Template().get( *(model, method, basis, loprop, freq) )
        for mol in self:
            for at in mol:
                Property.add_prop_from_template( at, templ )
            t1, t2, t3 = mol.get_euler()
            for at in mol:
                at.Property.transform_ut_properties( t3, t2, t1 )
            if loprop:
                mol.LoProp = True
            else:
                mol.LoProp = False
        self.Property = True

    def add(self, item ):
        if isinstance( item , Molecule ):
            self.add_mol( item )
        elif isinstance( item, Atom) :
            self.add_atom( item )
        else:
            logging.warning( 'Tried to pass other instance than Atom or Molecule to Cluster' )

    def add_mol(self, mol, ):
        if isinstance( mol , Molecule ):
            super( Cluster, self ).append( mol )
            mol.Cluster = self
        elif type( mol ) == list:
            for each in mol:
                each.in_mm = in_mm
                each.in_qm = in_qm
                each.in_qmmm = in_qmmm
                super( Cluster, self ).append( each )
                each.Cluster = each

    def add_atom(self, *at):
        for i, iat in enumerate(at):
            self.append( iat )
            iat.Cluster = self

    def reset_qm_mm(self):
        """Set all atoms to neither qm not mm"""
        for at in self:
            at.in_mm = False
            at.in_qm = False

    def set_qm_mm(self, N_qm = 1, N_mm = 0):
        """First set all waters to False for security """
        for i in self:
            i.in_qm = False
            i.in_mm = False

        """Set the first N_qm in qm region and rest N_mm in mm region"""
        for i in self[ 0 : N_qm  ]:
            i.in_qm = True
        for i in self[ N_qm  : N_qm + N_mm ]:
            i.in_mm = True
    def copy_cluster(self):
        tmp_c = Cluster()
        for res in self:
            tmp_c.add( res.copy_self() )
        return tmp_c

    @classmethod
    def from_pot_string( cls, pot, in_AA = False, out_AA = False ):
        c = cls()
        lines = pot.split('\n' )

        pat_N = re.compile( '@COORDINATES\n\s*(\d+)' )
        pat_xyz = re.compile(r'^\s*([A-Z])+\s+(-*\d*.+\d+)\s+(-*\d*.+\d+)\s+(-*\d*.+\d+) *$')
        for lin in lines:
            if pat_xyz.match(lin):
                x, y, z  = map(float, pat_xyz.match(lin).groups()[1:])
                c.add_atom( Atom( x = x, y = y, z = z, AA = in_AA ) )
        if in_AA and not out_AA:
            c.to_AU()
        elif out_AA and not in_AA:
            c.to_AA()
        return c


    def props_from_qm(self,
            tmpdir = None,
            dalpath = None,
            procs = 4,
            decimal = 5,
            maxl = 2,
            pol = 22,
            hyper = 2,
            method = 'hflin',
            env = os.environ,
            basis = ['ano-1 2', 'ano-1 4 3 1', 'ano-2 5 4 1' ],
            dalexe = None,
            basdir = '/home/x_ignha/repos/dalton/basis',
            log = None,
            keep_outfile = False,
            freq = None,
            ):
        """Put properties on all class/sublass of Molecules in Cluster"""
        for mol in [m for m in self if isinstance(m, Molecule)]:
            mol.props_from_qm( tmpdir = tmpdir,
                    dalpath = dalpath,
                    procs = procs,
                    decimal = decimal,
                    maxl = maxl, 
                    pol = pol,
                    hyper = hyper,
                    method = method,
                    basis = basis,
                    dalexe = dalexe,
                    basdir = basdir,
                    log = log,
                    keep_outfile = False,
                    freq = freq,
                    )


    def get_inp_string(self, method ='B3LYP', basis = "6-31+g*", procs= 8):
        """Write gaussian .inp file for geometry optimization"""
        st = r"%" + "Nprocshared=%d\n" %procs
        st += r"%" + "Mem=20MW\n"
        st += "#p %s/%s opt " %(method,basis)
        if not self.AA:
            st += "units=au " 
        st += '\n\ncomment\n\n'
        st += "%d %d\n" %( self.sum_property['charge'][0], 1 )
        for i in [at for mol in self for at in mol]:
            st += "%s %.5f %.5f %.5f\n" %(i.element, i.x, i.y, i.z ) 
        st+= '\n\n\n'
        return st


    @property
    def com(self):
        if len(self) == 0:return np.zeros(3)
        return sum([at.r*at.mass for mol in self for at in mol]) / sum([at.mass for mol in self for at in mol] )

    @property
    def p(self):
        return self.sum_property

    @property
    def sum_property(self):
        """
Return the sum properties of all molecules in cluster
        """
        el_dip = np.array([ (at.r-self.coc)*at.Property['charge'] for mol in self for at in mol])
        nuc_dip = np.array([ (at.r-self.coc)*charge_dict[at.element] for mol in self for at in mol])
        dip_lop = np.array([at.Property['dipole'] for mol in self for at in mol])
        dip = el_dip + nuc_dip
        dip_tot = (dip + dip_lop).sum(axis=0)
        p = Property()
        for mol in self:
            for at in mol:
                p += at.Property
        p['dipole'] = dip_tot
        return p

    def to_AA(self):
        if not self.AA:
            for i in self:
                i.to_AA()

    def to_AU(self):
        if self.AA:
            for i in self:
                i.to_AU()

    def translate(self, r):
        """
Translate everythin in cluster by r

.. code:: python

    >>> m = Cluster()
    >>> m.append( Atom(element = 'H', z = 1) )
    >>> m.append( Atom(element = 'H', z = 0) )
    >>> print m.com
    [0, 0, 0.5 ]
    >>> m.translate( [0, 3, 5] )
    >>> print m.com
    [0, 3, 5.5 ]
    
"""
        for at in [at for mol in self for at in mol]:
            at.x = r[0] + at.x 
            at.y = r[1] + at.y 
            at.z = r[2] + at.z 


if __name__ == '__main__':
    print "no main method for module implemented"
