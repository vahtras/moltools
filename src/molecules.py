#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""
The molecules modules serves as an interface to write water molecule input files using predefined geometries, to be used with the DALTON qm package.
"""

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt

import numpy as np
import re, os, itertools, h5py, warnings

from template import Template

a0 = 0.52917721092

charge_dict = {"H": 1.0, "C": 6.0, "N": 7.0, "O": 8.0, "S": 16.0}
# from TIP3P charge defs.
el_charge_dict = {"H": .417, "O": -0.834 , "X" : 0.417}
mass_dict = {"H": 1.008,  "C": 12.0, "N": 14.01, "O": 15.999, "S": 32.066,
    "X" : 1.008 }

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

        self["charge"] = np.zeros( 1 )
        self["dipole"] = np.zeros( 3 )
        self["quadrupole"] = np.zeros( 6 )
        self["alpha"] =  np.zeros( 6 ) 
        self["beta"] =  np.zeros( 10 ) 

    def copy_property(self):
        p = Property()
        p["charge"] =      self["charge"].copy()
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

    def potline(self, max_l , pol, hyper, fmt = "%.5f "):
        string = ""
        if 0  <= max_l :
            string += fmt % tuple(self["charge"] )
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
        if pol == 2 :
            string += fmt * 6 %( 
                    self["alpha"][0], self["alpha"][1], self["alpha"][2] ,
                    self["alpha"][3], self["alpha"][4], self["alpha"][5] )
            return string
        if pol == 22 :
            string += fmt * 6%( 
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
    def add_prop_from_template( at, wat_templ ):

        """
Puts properties read from the :ref:`template` module into the :ref:`atom` at.

    
    >>> temp = template.Template().get() #Default template
    >>> w = Water.get_standard() #Default water
    >>> Property.add_prop_from_template( w.o, temp )
    >>> print w.o["dipole"]
    [0.0, 0.0, 0.78719]

"""

        p = Property()
        for i, keys in enumerate( wat_templ ):
            if keys[0] == at.name:
                p[keys[1]] = np.array( wat_templ[ keys ] )
        at.Property = p
        at.Molecule.Property = True

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

    @staticmethod
    def rot_avg( beta, car1 = 2, car2 = 2, car3 = 2):
        """
        Requires euler.h5 binary file containing rotational angle products
        """
        b_new = np.zeros( (3,3,3,) )
        """given beta in molecular frame, convert to exp. reference"""
        vec = h5py.File('euler.h5','r')['data'].value
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

    >>> d = numpy.array( [ 1, 0, 0] )
    >>> print Rotator.transform_1( d, numpy.pi/2, 0, 0 )
    [ 0., -1., 0. ]
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
pdbname  X1       string
number   0        int
AA       True     bool
======== ======== ========

        """
#Element one-key char
        self.element = "X"

#Name is custom name, for water use O1, H2 (positive x ax), H3
        self.name = None
#Label is custom name, for water use O1, H2 (positive x ax), H3
        self.label = ""

        self.x = 0.0
        self.y = 0.0
        self.z = 0.0


# Use populate_bonds in class Molecule to attach all atoms to their neighbours
        self.bonds = {}
        self.angles = {}
        self.dihedral = {}

        self._q = None

        self.cluster = None

        self._res_id = 0
        self.atom_id = None

        self.in_water = False
        self.Molecule = Molecule()

#Property set to true if atoms have properties
        self.Property = Property()
        self.AA = True

        if kwargs != {}:
            self.AA = bool( kwargs.get( "AA", True ) )
            self.x = float( kwargs.get( "x", 0.0 ))
            self.y = float( kwargs.get( "y", 0.0 ))
            self.z = float( kwargs.get( "z", 0.0 ))
            self.element = kwargs.get( "element", "X" )
            self.name = kwargs.get( "name", "1-XXX-X1" )
            self.number = kwargs.get( "number", 0 )
            self.pdbname = kwargs.get( "pdbname", 'X1' )
        self._mass = None

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

    def translate(self, r ):
        self.x, self.y, self.z = r

    def copy_atom(self):
        a = Atom( **{'x':self.x, 'y':self.y, 'z':self.z,'AA':self.AA,
            'element':self.element,'name':self.name,'number':self.number,
            'pdbname':self.pdbname} )
        a._res_id = self.res_id
        a.atom_id = self.atom_id
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
        return None

    def potline(self, max_l, pol, hyper):
        return  "{0:4} {1:10f} {2:10f} {3:10f} ".format( \
                str(self.res_id), self.x, self.y, self.z ) + self.Property.potline( max_l, pol, hyper ) + "\n"

    def __str__(self):
        return "%s %f %f %f" %(self.name, self.x, self.y, self.z)

    def __sub__(self, other ):
        return self.r - other.r
    def __add__(self, other ):
        return self.r + other.r

    def get_array(self):
        return np.array( self.r ).copy()

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

    def to_au(self):
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
#Bond dict defined in angstromg, if molecule is in AU will be different later
        self.bonding_cutoff = { ('H','H') : 1.1,
                ('H','C') : 1.1,
                ('H','N') : 1.1,
                ('H','O') : 1.1,
                ('H','P') : 1.1,
                ('C','H') : 1.1,
                ('C','C') : 1.5,
                ('C','N') : 1.5,
                ('C','O') : 1.5,
                ('C','P') : 2.0,
                ('N','H') : 1.1,
                ('N','C') : 1.5,
                ('N','N') : 1.5,
                ('N','O') : 1.5,
                ('N','P') : 1.5,
                ('O','H') : 1.1,
                ('O','N') : 1.5,
                ('O','C') : 1.5,
                ('O','O') : 1.5,
                ('O','P') : 2.0,
                ('P','H') : 1.1,
                ('P','N') : 1.5,
                ('P','O') : 1.5,
                ('P','O') : 2.0,
                ('P','P') : 1.5,
            }

# Dictionary with bonds
        self.bond_dict = {}

#center will be defined for all molecules after all atoms are added
#depends on which molecule
        self.res_id = 0
        self._r = None
        self._com = None
        self.cluster = None
        self.no_hydrogens = True

# For plotting different elements:
        self.style = {"H":'wo', "O":'ro'}
        self.linewidth = {"H":25, "O":40}

# Make emptpy, beware that this causes molecules to give zero dipole momnet
# before template is loaded
        self.Property = None

#By default, AU 
        self.AA = False

#if supplied a dictionary with options, gather these in self.info
        self.info = {}
        if kwargs != {} :
            for i in kwargs:
                self.info[ i ] = kwargs[ i ]

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
        if self.AA:
            conv = 1.0
        else:
            conv = 1/a0
        for i in range(len(self)):
            for j in range( i + 1, len(self)):
                if self[i].dist_to_atom( self[j] ) < conv*self.bonding_cutoff[ (self[i].element, self[j].element) ]:

                    self[i].bonds[ self[j].name ] = self[j]
                    self[j].bonds[ self[i].name ] = self[i]

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

#Dipole moment
    @property
    def p(self):
        """
Return the dipole moment

.. code:: python

   >>> m = Molecule()
   >>> m.append( Atom(element = 'H', z = 1) )
   >>> m.append( Atom(element = 'O', z = 0) )
   >>> print m.p
   -0.834

"""
        #if self.Property:
        #    el_dip = np.array([ (at.r-self.coc)*at.Property['charge'] for at in self ])
        #    nuc_dip = np.array([ (at.r-self.coc)*charge_dict[at.element] for at in self ])
        #    dip_lop = np.array([at.Property['dipole'] for at in self])
        #    dip = el_dip + nuc_dip
        #    return dip.sum(axis=0)+ dip_lop.sum(axis=0)

        return np.array([at.r*at.q for at in self]).sum(axis=0)

    @property
    def sum_property(self):
        """
Return the sum properties of all properties in molecules

.. code:: python
    >>> wat
        """
        p = Property()
        for at in self:
            p += at.Property
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

    def translate_o(self, r):
        vec = r - self.o.r
        for at in self:
            at.x = vec[0] + at.x 
            at.y = vec[1] + at.y 
            at.z = vec[2] + at.z 

    def translate_coc(self, r):
        vec = r - self.coc
        for at in self:
            at.x = vec[0] + at.x 
            at.y = vec[1] + at.y 
            at.z = vec[2] + at.z 

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

        if self.Property:
            pass

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


    def plot(self ):
        """
Plot the molecule in a 3D frame

.. code:: python

    >>> m = Molecule()
    >>> m.append( Atom(element = 'H', x = 1, z = 1) )
    >>> m.append( Atom(element = 'H', x =-1, z = 1) )
    >>> m.append( Atom(element = 'O', z = 0) )
    >>> m.plot()
    
"""

#Plot water molecule in green and  nice xyz axis
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d' )
        ax.plot( [0, 1, 0, 0, 0, 0], [0,0 ,0,1,0,0], [0,0,0,0,0,1] )
        ax.text( 1.1, 0, 0, "X", color = 'red' )
        ax.text( 0, 1.1, 0, "Y", color = 'red' )
        ax.text( 0, 0, 1.1, "Z", color = 'red' )
        x = self.coc[0]
        y = self.coc[1]
        z = self.coc[2]
        p = self.p

        ax.plot( [x,x+p[0]], [y,y+p[1]], [z,z+p[2]], '-', linewidth = 3 )
        for i in self:
            ax.plot( [i.x], [i.y], [i.z], self.style[i.element], linewidth= self.linewidth[i.element] )
        ax.set_zlim3d( -5,5)
        plt.xlim(-5,5)
        plt.ylim(-5,5)
        plt.show()

    def get_mol_string(self, basis = ("ano-1 2 1", "ano-1 3 2 1" ) ):
        if len( basis ) > 1:
            el_to_rowind = {"H" : 0, "C" : 1, "O" : 1, "N" : 1  }
        else:
            el_to_rowind = {"H" : 0, "C" : 0, "O" : 0, "N" : 0 }
        st = ""
        s_ = ""
        if self.AA: s_ += " Angstrom"
        uni = Molecule.unique([ at.element for at in self])
        st += "ATOMBASIS\n\n\nAtomtypes=%d Charge=0 Nosymm%s\n" %(len(uni), s_)
        for el in uni:
            st += "Charge=%s Atoms=%d Basis=%s\n" %( str(charge_dict[el]),
                    len( [all_el for all_el in self if (all_el.element == el)] ),
                    basis[ el_to_rowind[el] ])
            for i in [all_el for all_el in self if (all_el.element == el) ]:
                st += "%s %.5f %.5f %.5f\n" %(i.element, i.x, i.y, i.z ) 
        return st

    @staticmethod
    def unique(arr):
        tmp = []
        for i in arr:
            if i not in tmp:
                tmp.append(i)
        return tmp

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
    def from_mol_file( molfile, AA = False):
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
        pat_xyz = re.compile(r'^\s*(\S+)\s+(-*\d*\.{1}\d+)\s+(-*\d*\.{1}\d+)\s+(-*\d*\.{1}\d+) *$')
        tmp_molecule = Molecule()
        for i in open( molfile ).readlines():
            if pat_xyz.search(i):
                matched = pat_xyz.match(i).groups()
                pd = matched[0].split('-')[-1]
                kwargs = { "name" :  matched[0], "x" : matched[1],
                        "y" : matched[2], "z" : matched[3], "AA" : AA,
                        "pdbname" : pd }
                tmpAtom = Atom( **kwargs )
                tmp_molecule.append( tmpAtom )
        return tmp_molecule



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
            print "Error: Molecule.from_xyz recieved non-xyz file: %s" %f
            raise SystemExit

        fil = open(f).readlines()
        m = Molecule()
        for ind, i in enumerate( fil ):
            if ind in [0, 1]: continue

            elem = i.split()[0]
            x = i.split()[1]
            y = i.split()[2]
            z = i.split()[3]
            m.append( Atom( **{"element":elem,
                "x" : x,
                "y" : y,
                "z" : z,
                "AA" : in_AA,
                }))

        if in_AA:
            if not out_AA:
                m.to_AU()
        return m

    def to_AU(self):
        assert self.AA == True
        for at in self:
            at.x = at.x / a0
            at.y = at.y / a0
            at.z = at.z / a0
        self.AA = False

    def to_AA(self):
        assert self.AA == False
        for at in self:
            at.x *= a0
            at.y *= a0
            at.z *= a0
        self.AA = True

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

    def copy_water(self):
        w = Water()
        [w.append(i.copy_atom()) for i in self]
        return w
    @staticmethod
    def get_standard( AA = False):
        """
Return water molecule from specified template with :math:`r=0.972` Angstrom and 
:math:`\\theta=104.5` degrees.

.. code:: python

    >>> m = Water.get_standard()

"""
#Geometrical parameters
        center = [0, 0, 0]
        r_oh =  104.5
        a_hoh = 0.9720
        r_oh = r_oh / a0
        d = (90 - a_hoh/2 ) * np.pi / 180
        origin = np.array( [ 0, 0, 0] )

        h1 = Atom( **{ "AA" : AA, "element" : "H"} )
        h2 = Atom( **{ "AA" : AA, "element" : "H"} )
        o =  Atom( **{ "AA" : AA, "element" : "O"} )

        o.x = center[0]
        o.y = center[1]
        o.z = center[2] 

        h1.x = (center[0] + r_oh * np.cos(d))
        h1.y = center[1] 
        h1.z = (center[2] + r_oh* np.sin(d))

        h2.x = (center[0] - r_oh * np.cos(d)) 
        h2.y = center[1] 
        h2.z = (center[2] + r_oh* np.sin(d))

        w = Water()
        w.append( o )
        w.append( h1 )
        w.append( h2 )
        
        return w

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
    >>> m.add_atom( Atom( z : 0.11, element : 'H' ) )
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
        super( Water , self).append(atom)

#Define water center, by default set it to center of nuclei
        if (self.h1 and self.h2 and self.o):
            pass

        if self.res_id:
            if self.res_id != atom.res_id:
                print "Tried to add %s to %s, exiting" %(atom, self)
                raise SystemExit
        else:
#Initialize water res_id from atomic res_id
            self.res_id = atom.res_id

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

    def __str__(self):
        return "WAT" + str(self.res_id) 
    
    def exclists(self):
        tmp = []
        uniq = []
        for i in itertools.permutations( [at.number for at in self], len(self) ):
            if i[0] not in uniq:
                tmp.append(i)
                uniq.append( i[0] )
        return tmp

    def dist_to_point( self , point ):
        return np.sqrt(np.sum((self.coo - np.array(point))**2))

    def dist_to_water(self, other):
        return np.sqrt(np.sum((self.coo - other.coo)**2) )

    def get_euler(self):
        """
Returns the 3 euler angles required to rotate the water to given coordinate system.
The return values are ordered in :math:`\\rho_1`, :math:`\\rho_2` and :math:`\\rho_3`.

.. code:: python

    >>> w = Water()
    >>> w.append( Atom( x = 1, z = 1, element = 'H' ) )
    >>> w.append( Atom( x =-1, z = 1, element = 'H' ) )
    >>> w.append( Atom( x = 0, z = 1, element = 'O' ) )
    >>> r1, r2, r3 = w.get_euler()
    >>> print r1
    0.0


        """

        H1 = self.h1.r.copy()
        H2 = self.h2.r.copy()
        O1 = self.o.r.copy()

        dip = (-0.5*O1 + 0.25*H1 + 0.25 *H2).copy()

        origin = O1.copy()
        H1, H2, O1 = H1 - origin, H2 - origin, O1 - origin

        theta1 = np.arctan2( dip[1], dip[0])

        H1 =  np.dot( Rotator.get_Rz_inv( theta1 ) , H1 )
        H2 =  np.dot( Rotator.get_Rz_inv( theta1 ) , H2 )
        O1 =  np.dot( Rotator.get_Rz_inv( theta1 ) , O1 )

        dip = np.dot( Rotator.get_Rz_inv( theta1 ) , dip )

#Rotate by theta around y axis so that the dipole is in the z axis 
        theta2 = np.arctan2( -dip[0], dip[2] )

        H1 =  np.dot( Rotator.get_Ry( theta2 ) , H1 )
        H2 =  np.dot( Rotator.get_Ry( theta2 ) , H2 )
        O1 =  np.dot( Rotator.get_Ry( theta2 ) , O1 )

        dip = np.dot( Rotator.get_Ry( theta2 ) , dip )

#Rotate around Z axis so that hydrogens are in xz plane.
        if H2[1] >0:
            xc = H2[0]
            yc = H2[1]
        else:
            xc = H1[0]
            yc = H1[1]
        theta3 = np.arctan2( yc , xc)

        def eq(a, b, thr = 0.0001): 
            if abs(a-b) < thr:return True
            else: return False

        return theta3, theta2, theta1

    def inv_rotate(self):
        """rotate all atom positions by
        1) inverse Z rotation by t1
        2) positive Y rotation by t2
        3) inverse Z rotation by t3
        """

        t1, t2, t3 = self.get_euler()
# Place water molecule with O in origo, and rotate it so hydrogens in xz plane
        H1 = self.h1.get_array() ; H2 = self.h2.get_array() ; O = self.o.get_array()
        TMP = self.o.get_array()
        H1 -= TMP ; H2 -= TMP; O -= TMP

        H1 = np.dot( Rotator.get_Rz_inv(t3) , H1 )
        H1 = np.dot( Rotator.get_Ry(t2) , H1 )
        H1 = np.dot( Rotator.get_Rz_inv(t1) , H1 )

        H2 = np.dot( Rotator.get_Rz_inv(t3) , H2 )
        H2 = np.dot( Rotator.get_Ry(t2) , H2 )
        H2 = np.dot( Rotator.get_Rz_inv(t1) , H2 )

        O = np.dot( Rotator.get_Rz_inv(t3) , O )
        O = np.dot( Rotator.get_Ry(t2) , O )
        O = np.dot( Rotator.get_Rz_inv(t1) , O )
#Put back in oxygen original point
        H1 += TMP ; H2 += TMP; O += TMP

        self.h1.x = H1[0] ;self.h1.y = H1[1] ;self.h1.z = H1[2] 
        self.h2.x = H2[0] ;self.h2.y = H2[1] ;self.h2.z = H2[2] 
        self.o.x  =  O[0] ;  self.o.y = O[1] ;  self.o.z = O[2] 

    def rotate(self, t1, t2, t3):
        """Rotate all coordinates by t1, t2 and t3
        first Rz with theta1, then Ry^-1 by theta2, then Rz with theta 3

        R all in radians

        """
# Place water molecule in origo, and rotate it so hydrogens in xz plane
        self.inv_rotate()

        H1 = self.h1.get_array() ; H2 = self.h2.get_array() ; O = self.o.get_array()
        TMP = self.o.get_array()
        H1 -= TMP ; H2 -= TMP; O -= TMP

# Rotate with angles t1, t2, t3

        H1 = np.dot( Rotator.get_Rz(t1) , H1 )
        H1 = np.dot( Rotator.get_Ry_inv(t2) , H1 )
        H1 = np.dot( Rotator.get_Rz(t3) , H1 )

        H2 = np.dot( Rotator.get_Rz(t1) , H2 )
        H2 = np.dot( Rotator.get_Ry_inv(t2) , H2 )
        H2 = np.dot( Rotator.get_Rz(t3) , H2 )

        O = np.dot( Rotator.get_Rz(t1) , O )
        O = np.dot( Rotator.get_Ry_inv(t2) , O )
        O = np.dot( Rotator.get_Rz(t3) , O )

#Put back in oxygen original point
        H1 += TMP ; H2 += TMP; O += TMP

        self.h1.x = H1[0] ;self.h1.y = H1[1] ;self.h1.z = H1[2] 
        self.h2.x = H2[0] ;self.h2.y = H2[1] ;self.h2.z = H2[2] 
        self.o.x  =  O[0] ;  self.o.y = O[1] ;  self.o.z = O[2] 

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
            self.res_id = atom.res_id
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

class Ethane(list):
    def __init__(self):
        pass

class Cluster(list):
    """
**Molecule container which groups molecules into quantum mechanics, molecular mechanics, and qm/mm regions for easy in generating input files for QM/MM.**
"""
    def __init__(self, *args, **kwargs):
        """ Typical list of molecules """
        self.Property = None
        self.atom_list = []

    def __str__(self):
        return " ".join( [ str(i) for i in self ] )

    def append(self, mol, in_mm = False, in_qm = False,
            in_qmmm = False):
        mol.in_mm = in_mm
        mol.in_qm = in_qm
        mol.in_qmmm = in_qmmm

        super( Cluster, self ).append( mol )

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
    
    @property
    def p(self):
        if self.Property:
            el_dip = np.array([ (at.r-self.coc)*at.Property['charge'] for mol in self for at in mol])
            nuc_dip = np.array([ (at.r-self.coc)*charge_dict[at.element] for mol in self for at in mol])
            dip_lop = np.array([at.Property['dipole'] for mol in self for at in mol])
            dip = el_dip + nuc_dip
            return dip.sum(axis=0)+ dip_lop.sum(axis=0)

        return np.array([at.r*at.q for mol in self for at in mol]).sum(axis=0)



# Specifi

    @property
    def coc(self):
        if self.Property:
            pass
    #obj should be atom
        return sum( [at.r * charge_dict[at.element] for mol in self for at in mol])\
                /sum( map(float,[charge_dict[at.element] for mol in self for at in mol]) )


    def plot(self ):
        """
Plot all the molecule in a 3D frame in the cluster

.. code:: python

    >>> m = Molecule()
    >>> m.append( Atom(element = 'H', x = 1, z = 1) )
    >>> m.append( Atom(element = 'H', x =-1, z = 1) )
    >>> m.append( Atom(element = 'O', z = 0) )
    >>> m.plot()
    
"""

#Plot water molecule in green and  nice xyz axis
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d' )
        ax.plot( [0, 1, 0, 0, 0, 0], [0,0 ,0,1,0,0], [0,0,0,0,0,1] )
        ax.text( 1.1, 0, 0, "X", color = 'red' )
        ax.text( 0, 1.1, 0, "Y", color = 'red' )
        ax.text( 0, 0, 1.1, "Z", color = 'red' )
        x = self.coc[0]
        y = self.coc[1]
        z = self.coc[2]
        p = self.p

        ax.plot( [x,x+p[0]], [y,y+p[1]], [z,z+p[2]], '-', linewidth = 3 )
        for i in self:
            for j in i:
                ax.plot( [j.x], [j.y], [j.z], j.Molecule.style[j.element], linewidth= j.Molecule.linewidth[j.element] )
        ax.set_zlim3d( -5,5)
        plt.xlim(-5,5)
        plt.ylim(-5,5)
        plt.show()



    def get_qm_mol_string(self, basis = ("ano-1 2 1", "ano-1 3 2 1" ) , AA = False):
# If basis len is more than one, treat it like molecular ano type
# Where first element is for first row of elements

        if len( basis ) > 1:
            # Set row index number to periodic table one
            el_to_rowind = {"H" : 0, "C" : 1, "O" : 1, "N" : 1  }
        else:
            # Keep all 0, since basis is only len 1
            el_to_rowind = {"H" : 0, "C" : 0, "O" : 0, "N" : 0 }

        st = ""
        comm1 = "QM: " + " ".join( [ str(m) for m in self if m.in_qm] )[:72]
        comm2 = "MM: " + " ".join( [ str(m) for m in self if m.in_mm] )[:73]
        uni = Molecule.unique([ at.element for mol in self for at in mol if mol.in_qm])
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
    def get_pe_pot_string( self, max_l = 1, pol = 2, hyp = 0, out_AA = False ):
        self.order_mm_atoms()
        st = r'!%s' % (self ) + '\n'
        st += r'@COORDINATES' + '\n'
        st += '%d\n' % sum([len(i) for i in self if i.in_mm ])
        if out_AA:
            st += "AA\n"
        else:
            st += "AU\n"
        #st += '%d\n' %len(mol)
        for mol in [m for m in self if m.in_mm]:
            for at in mol:
                st += "%s %.5f %.5f %.5f\n" % (at.number, \
                        at.x, at.y, at.z )

        st += r'@MULTIPOLES'  + '\n'
        if max_l >= 0:
            st += 'ORDER 0\n'
            st += '%d\n' % sum([len(i) for i in self if i.in_mm ])
            for mol in [m for m in self if m.in_mm]:
                for at in mol:
                    st += "%s %.5f\n" % (tuple( [at.number] ) + tuple( at.Property["charge"] )  )
        if max_l >= 1:
            st += 'ORDER 1\n'
            st += '%d\n' % sum([len(i) for i in self if i.in_mm ])
            for mol in [m for m in self if m.in_mm]:
                for at in mol:
                    st += "%s %.5f %.5f %.5f\n" % ( tuple([at.number]) + tuple(at.Property["dipole"])) 

        st += r'@POLARIZABILITIES' + '\n'
        st += 'ORDER 1 1\n'
        st += '%d\n' % sum([len(i) for i in self if i.in_mm ])
        if pol % 2 == 0:
            for mol in [m for m in self if m.in_mm]:
                #st += 'ORDER 1 1\n'
                #st += '%d\n' % len( mol )
                for at in mol:
                    st += "%s %.5f %.5f %.5f %.5f %.5f %.5f\n" % ( tuple([at.number]) + tuple(at.Property["alpha"])) 

        st += 'EXCLISTS\n%d %d\n' %( sum([len(i) for i in self if i.in_mm ])
 , len(mol))
        for mol in [m for m in self if m.in_mm]:
            for each in mol.exclists():
                ls = ""
                for ind in each:
                    ls += "%s " %ind
                ls += '\n'
                st += ls

        return st
# This is the old *QMMM input style in dalton
    def get_qmmm_pot_string( self, max_l = 1,
            pol = 22,
            hyp = 1,
            in_AA = False,
#Set ignore_qmmm to false to only write qmmm .pot file for molecues in mm region
            ignore_qmmm = True ):
        if in_AA:
            st = "AA\n"
        else:
            st = "AU\n"
# Old qmmm format requires integer at end to properly read charges
        if ignore_qmmm:
            st += "%d %d %d %d\n" % (sum([len(i) for i in self ]), 
                    max_l, pol, 1 )
            st += "".join( [at.potline(max_l, pol, hyp) for mol in self for at in mol ] )
        else:
            st += "%d %d %d %d\n" % (sum([len(i) for i in self if i.in_mm ]), 
                    max_l, pol, 1 )
            st += "".join( [at.potline(max_l, pol, hyp) for mol in self for at in mol if mol.in_mm] )
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
        for mol in [m for m in self if m.in_mm]:
            for at in mol:
                at.number = str(cnt)
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

    >>> c = Cluster.get_water_cluster( 'somefile.mol' , in_AA = False, out_AA = False, N_waaters = 10 )
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
#loop over oxygen and hydrogen and if they are closer than 1 A add them to a water
        waters = []
        cnt = 1
        if fname.endswith( ".xyz" ) or fname.endswith(".mol"):
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
                tmp.res_id = cnt
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

        if in_AA:
            if not out_AA:
                for wat in c:
                    wat.to_AU()
        if not in_AA:
            if out_AA:
                for wat in c:
                    wat.to_AA()
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
                at.Property.transform_ut_properties( t1, t2, t3 )

    def add_mol(self, mol, *args):
        #if type(args):
        self.append( mol )
        mol.cluster = self

    def add_atom(self, *at):
        for i, iat in enumerate(at):
            self.append( iat )
            iat.cluster = self

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
        [tmp_c.add_mol(wat.copy_water()) for wat in self]
        return tmp_c

    @property
    def sum_property(self):
        """
Return the sum properties of all molecules in cluster

.. code:: python
    >>> wat
        """
        p = Property()
        for mol in self:
            for at in mol:
                p += at.Property
        return p

if __name__ == '__main__':
    from use_generator import *
    from gaussian import *
    m1 = Generator().get_mol(model = 'spc' )
    m2 = Generator().get_mol( center = [0,0, 10 ], model = 'spc' )
    m2.res_id = 2
    c = Cluster()
    c.add_mol(m1)
    c.add_mol(m2)
    
    t = Template().get( model = 'SPC' )
#    for m in [m1, m2]:
#        for at in m:
#            Property.add_prop_from_template( at, t )
    c.update_water_props( model = 'SPC', dist = True )
    g = GaussianQuadrupoleList.from_string( c.get_qmmm_pot_string( ignore_qmmm = True ) )
    g.solve_scf()
    print g.beta()

