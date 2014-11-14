#!/usr/bin/env python
#-*- coding: utf-8 -*-

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt

import itertools

import numpy as np
import math as m
import re, os, ut

a0 = 0.52917721092

charge_dict = {"H": 1.0, "C": 6.0, "N": 7.0, "O": 8.0, "S": 16.0}
el_charge_dict = {"H": 0.35,  "O": -0.70, }
mass_dict = {"H": 1.008,  "C": 12.0, "N": 14.01, "O": 15.999, "S": 32.066}

def tensor_to_ut( beta ):
# naive solution, transforms matrix B[ (x,y,z) ][ (xx, xy, xz, yy, yz, zz) ] into array
# Symmtrized UT array    B[ (xxx, xxy, xxz, xyy, xyz, xzz, yyy, yyz, yzz, zzz) ]
    new = np.array( (10) )
    new[ 0 ] = beta[0, 0, 0]
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


class Property( dict ):
    def __init__(self):

        self.max_l = 2
        self.pol = 22
        self.hyper = 2

    def __add__(self, other):
        tmp = {}
        for i, prop in enumerate(self):
            tmp[prop] = np.array( self[prop] ) + np.array(other[prop] )
        return tmp
    def __ladd__(self, other):
        tmp = {}
        for i, prop in enumerate(self):
            tmp[prop] = np.array( self[prop] ) + np.array(other[prop] )
        return tmp
    def __radd__(self, other):
        tmp = {}
        for i, prop in enumerate(self):
            tmp[prop] = np.array( self[prop] ) + np.array(other[prop] )
        return tmp


    def __str__(self):
        return "%.5f %.5f %.5f %.5f" % tuple( self["charge"] + self["dipole"]  )
        #return "%.5f %.5f %.5f %.5f  %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f \n" %(
        #        self["charge"], self["dipole"][0], self["dipole"][1], self["dipole"][2],

        #        self["quadrupole"][0], self["quadrupole"][1], self["quadrupole"][2],
        #        self["quadrupole"][3],self["quadrupole"][4], self["quadrupole"][5],

        #        self["alpha"][0], self["alpha"][1], self["alpha"][2], self["alpha"][3],
        #        self["alpha"][4], self["alpha"][5],
        #        
        #        self["beta"][0], self["beta"][1], self["beta"][2],
        #        self["beta"][3], self["beta"][4], self["beta"][5],
        #        self["beta"][6], self["beta"][7], self["beta"][8],
        #        self["beta"][9]
        #        )

    def potline(self, max_l , pol, hyper):
        string = ""
        if 0  <= max_l :
            string += "%.5f " % tuple(self["charge"] )
        if max_l >= 1 :
            string += "%.5f %.5f %.5f " %( self["dipole"][0], self["dipole"][1], self["dipole"][2] )
        if max_l >= 2 :
            string += "%.5f %.5f %.5f %.5f %.5f %.5f " %( 
                    self["quadrupole"][0], self["quadrupole"][1], self["quadrupole"][2] ,
                    self["quadrupole"][3], self["quadrupole"][4], self["quadrupole"][5] )
        if pol == 2 :
            string += "%.5f %.5f %.5f %.5f %.5f %.5f " %( 
                    self["alpha"][0], self["alpha"][1], self["alpha"][2] ,
                    self["alpha"][3], self["alpha"][4], self["alpha"][5] )

        if pol == 22 :
            string += "%.5f %.5f %.5f %.5f %.5f %.5f " %( 
                    self["alpha"][0], self["alpha"][1], self["alpha"][2] ,
                    self["alpha"][3], self["alpha"][4], self["alpha"][5] )

        if hyper == 1 :
            string += "%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f" %( 
                    self["beta"][0], self["beta"][1], self["beta"][2] ,
                    self["beta"][3], self["beta"][4], self["beta"][5] ,
                    self["beta"][6], self["beta"][7], self["beta"][8] ,
                    self["beta"][9])

        return string

    @staticmethod
    def add_prop_from_template( at, wat_templ ):
        p = Property()
        for i, keys in enumerate( wat_templ ):
            if keys[0] == at.name:
                p[keys[1]] = wat_templ[ keys ]
        at.Property = p

    @staticmethod
    def transform_ut_properties( prop, t1, t2, t3):
        assert isinstance( prop, dict )
        if prop.has_key( "dipole" ):
            prop["dipole"] = Water.transform_dipole( prop["dipole"] , t1, t2, t3 )
        if prop.has_key( "quadrupole" ):
            prop["quadrupole"] = prop.transform_ut_quadrupole( prop["quadrupole"], t1, t2, t3 )

        if prop.has_key( "alpha" ):
            prop["alpha"] = prop.transform_ut_alpha( prop["alpha"],t1, t2, t3 )

        if prop.has_key( "beta" ):
            prop["beta"] = prop.transform_ut_beta( prop["beta"], t1, t2, t3 )

    @staticmethod
    def transform_ut_quadrupole( quad, t1, t2 ,t3 ):
        tmp_Q = Water.ut_2_square( quad )
        tmp_Q = Water.transform_alpha( tmp_Q , t1 ,t2 ,t3 )
        tmp_Q = Water.square_2_ut( tmp_Q )
        return tmp_Q

    @staticmethod
    def transform_ut_alpha( alpha, t1, t2 ,t3 ):
        tmp_a = Water.ut_2_square( alpha )
        tmp_a = Water.transform_alpha( tmp_a , t1 ,t2 ,t3 )
        tmp_a = Water.square_2_ut( tmp_a )
        return tmp_a

    @staticmethod
    def transform_ut_beta( beta, t1, t2 ,t3 ):
        tmp_b = Water.ut_3_square( beta )
        tmp_b = Water.transform_beta( tmp_b, t1 ,t2 ,t3 )
        tmp_b = Water.square_3_ut( tmp_b )
        return  tmp_b




class Atom(object):
    """ By default in Atomic units for coordinates """
    def __init__(self, *args, **kwargs ):
#Element one-key char
        self.element = None

#Name is custom name, for water use O1, H2 (positive x ax), H3
        self.name = None

        self.x = None
        self.y = None
        self.z = None

        self._q = None

        self.cluster = None

        self.res_id = 0
        self.atom_id = None

        self.in_water = False
        self.Water = None
        self.Property = None
        self.AA = False

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

    @property
    def r(self):
        return np.array( [ self.x, self.y, self.z ] )

    @property
    def q(self):
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

    def potline(self, max_l, pol, hyper):
        return  "{0:4}{1:10f}{2:10f}{3:10f} ".format( \
                str(self.res_id), self.x, self.y, self.z ) + self.Property.potline( max_l, pol, hyper ) + "\n"

    def __str__(self):
        return "%s %f %f %f" %(self.name, self.x, self.y, self.z)

    def __sub__(self, other ):
        return self.r - other.r

    def get_array(self):
        return np.array( [self.x , self.y, self.z ] ).copy()

    def dist_to_atom(self, other):
        return np.sqrt( (self.x - other.x)**2 + (self.y -other.y)**2 + (self.z -other.z)**2 )
    def dist_to_point(self, other):
        return np.sqrt( (self.x - other[0])**2 + (self.y -other[1])**2 + (self.z -other[2])**2 )

    def to_au(self):
        self.x /= a0
        self.y /= a0
        self.z /= a0

    def to_AA(self):
        self.x *= a0
        self.y *= a0
        self.z *= a0

class Molecule( list ):
    """General molecule has general methods to obtain euler angles, 
    All molecules inherits from this one"""

    def __init__(self , *args, **kwargs):

#center will be defined for all molecules after all atoms are added
#depends on which molecule
        self.res_id = 0
        self._r = None
        self._com = None

#By default, AU 
        self.AA = False
        self.Property = None
        self.no_hydrogens = True
#if supplied a dictionary with options, gather these in self.info
        self.info = {}
        if kwargs != {} :
            for i in kwargs:
                self.info[ i ] = kwargs[ i ]
#Dipole moment
    @property
    def p(self):
        return np.array([at.r*at.q for at in self]).sum(axis=0)

#Vector pointing to center of atom position
    @property
    def r(self):
        return  np.array([at.r for at in self]).sum(axis = 0) / len(self)

    def translate(self, r):
        vec = r - self.o.r
        for at in self:
            at.x = vec[0] + at.x 
            at.y = vec[1] + at.y 
            at.z = vec[2] + at.z 

    @property
    def com(self):
        if self._com is not None:
            return self._com
        self._com = np.array([at.mass*at.r for at in self]).sum(axis=0) / np.array([at.mass for at in self]).sum()

        return self._com

    @staticmethod
    def get_Rz( theta ):
        vec = np.array(    [[ m.cos(theta),-m.sin(theta), 0],
                            [ m.sin(theta), m.cos(theta), 0],
                            [ 0,    0,  1]])
        return vec
    @staticmethod
    def get_Rz_inv( theta ):
        vec = np.array(     [[ m.cos(theta), m.sin(theta), 0],
                            [ -m.sin(theta), m.cos(theta), 0],
                            [ 0,             0,            1]])
        return vec
    @staticmethod
    def get_Ry( theta ):
        vec = np.array(    [[ m.cos(theta),0, m.sin(theta)],
                            [ 0,    1,  0],
                            [ -m.sin(theta), 0, m.cos(theta)]])
        return vec
    @staticmethod
    def get_Ry_inv( theta ):
        vec = np.array(    [[ m.cos(theta),0, -m.sin(theta)],
                            [ 0,    1,  0],
                            [ m.sin(theta), 0, m.cos(theta)]])
        return vec
    @staticmethod
    def transform_dipole( qm_dipole, t1, t2, t3 ):
        d_new1 = np.zeros([3]) #will be returned
        d_new2 = np.zeros([3]) #will be returned
        d_new3 = np.zeros([3]) #will be returned

        rz  = Water.get_Rz( t1 )
        ryi = Water.get_Ry_inv( t2 )
        rz2 = Water.get_Rz( t3 )

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
    def transform_quadrupole( qm_quadrupole, t1, t2 , t3 ):
        a_new1 = np.zeros([3,3]) #will be calculated
        a_new2 = np.zeros([3,3]) #will be calculated
        a_new3 = np.zeros([3,3]) #will be calculated

        rz = self.get_Rz( t1 )
        ryi = self.get_Ry_inv( t2 )
        rz2 = self.get_Rz( t3 )
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
    def transform_alpha( qm_alpha, t1, t2 , t3 ):
        a_new1 = np.zeros([3,3]) #will be calculated
        a_new2 = np.zeros([3,3]) #will be calculated
        a_new3 = np.zeros([3,3]) #will be calculated

        rz  = Water.get_Rz( t1 )
        ryi = Water.get_Ry_inv( t2 )
        rz2 = Water.get_Rz( t3 )

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
    def transform_beta( qm_beta, t1, t2, t3 ):
        b_new1 = np.zeros([3,3,3]) #will be calculated
        b_new2 = np.zeros([3,3,3]) #will be calculated
        b_new3 = np.zeros([3,3,3]) #will be calculated

        rz =  Water.get_Rz( t1 )
        ryi = Water.get_Ry_inv( t2 )
        rz2 = Water.get_Rz( t3 )

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
        for index, (i, j ) in enumerate( ut.upper_triangular(2) ):
            tmp_a[ index ] = (alpha[i, j] + alpha[ j, i]) / 2
        return tmp_a

    def dist_to_mol(self, other):
        xyz1 = self.com
        xyz2 = other.com
        return m.sqrt( (xyz1[0] - xyz2[0])**2 + \
            (xyz1[1] - xyz2[1])**2 + (xyz1[2] - xyz2[2])**2 )




    @staticmethod
    def square_3_ut(beta):
        assert beta.ndim == 3
        tmp_b = np.zeros( 10 )
        for index, (i, j, k ) in enumerate( ut.upper_triangular(3) ):
            tmp_b[ index ] = ( \
                    beta[i, j, k] + beta[i, k, j] + \
                    beta[j, i, k] + beta[j, k, i] + \
                    beta[k, i, j] + beta[k, j, i] )/ 6
        return tmp_b

    @staticmethod
    def ut_2_square( alpha):
        assert len(alpha) == 6
        tmp_a = np.zeros( (3,3, ))
        for index, val in enumerate( ut.upper_triangular(2) ) :
            tmp_a[ val[0], val[1] ] = alpha[ index ]
            tmp_a[ val[1], val[0] ] = alpha[ index ]
        return tmp_a

    @staticmethod
    def ut_3_square( beta ):
        assert len(beta) == 10
        tmp_b = np.zeros( (3,3,3, ))
        for index, (i, j, k ) in enumerate( ut.upper_triangular(3) ) :
            tmp_b[ i, j ,k] = beta[ index ]
            tmp_b[ i, k ,j] = beta[ index] 
            tmp_b[ j, i, k] = beta [ index ]
            tmp_b[ j, k, i] = beta [ index ]
            tmp_b[ k, i, j] = beta [ index ]
            tmp_b[ k, j, i] = beta [ index ]
        return tmp_b
    def plot(self ):
#Plot water molecule in green and  nice xyz axis
        O1, H1, H2 = self.o, self.h1, self.h2
        fig = plt.figure()
        dip = self.get_dipole()
        ax = fig.add_subplot(111, projection='3d' )
        ax.plot( [0, 1, 0, 0, 0, 0], [0, 0,0,1,0,0], [0,0,0,0,0,1] )
        ax.plot( [O1.x,O1.x + dip[0] ] ,[ O1.y,O1.y+dip[1]],[O1.z,O1.z+dip[2]] ,'-',color="black")
        ax.scatter( [H1.x], [ H1.y] ,[ H1.z], s=25, color='red')
        ax.scatter( [H2.x], [ H2.y] ,[ H2.z], s=25, color='red')
        ax.scatter( [O1.x], [ O1.y] ,[ O1.z], s=50, color='blue')
        ax.set_zlim3d( -5,5)
        plt.xlim(-5,5)
        plt.ylim(-5,5)
        plt.show()
         
 
    @staticmethod
    def transform_dist_quadrupole( qm_quadrupole, t1, t2, t3 ):
        upper0 = np.array( qm_quadrupole[0] )
        upper1 = np.array( qm_quadrupole[1] )
        upper2 = np.array( qm_quadrupole[2] )
        assert upper0.shape == (6,)
        assert upper1.shape == (6,)
        assert upper2.shape == (6,)
        a0 = np.zeros((3, 3))
        a1 = np.zeros((3, 3))
        a2 = np.zeros((3, 3))

        for ij, (i, j) in enumerate(ut.upper_triangular(2)):

            aij0 = upper0[ij]
            aij1 = upper1[ij]
            aij2 = upper2[ij]

            a0[i, j] = aij0
            a0[j, i] = aij0

            a1[i, j] = aij1
            a1[j, i] = aij1

            a2[i, j] = aij2
            a2[j, i] = aij2

        a_new = np.zeros([3,3,3]) #will be returned
        a_new[0, :, :] = self.transform_alpha( a0, t1, t2, t3 )
        a_new[1, :, :] = self.transform_alpha( a1, t1, t2, t3 )
        a_new[2, :, :] = self.transform_alpha( a2, t1, t2, t3 )
        return a_new


    @staticmethod
    def transform_dist_dipole( self, qm_dipole, t1, t2, t3):
#Input qm_dipole is 3 x 3 matrix row is atom col is px, py, pz
        d_new = np.zeros([3,3]) #will be returned
        d_new[0, :] = self.transform_dipole( qm_dipole[0], t1, t2, t3 )
        d_new[1, :] = self.transform_dipole( qm_dipole[1], t1, t2, t3 )
        d_new[2, :] = self.transform_dipole( qm_dipole[2], t1, t2, t3 )
        return d_new

    @staticmethod
    def transform_dist_alpha( qm_alpha, t1, t2, t3 ):
        upper0 = np.array( qm_alpha[0] )
        upper1 = np.array( qm_alpha[1] )
        upper2 = np.array( qm_alpha[2] )
        assert upper0.shape == (6,)
        assert upper1.shape == (6,)
        assert upper2.shape == (6,)
        a0 = np.zeros((3, 3))
        a1 = np.zeros((3, 3))
        a2 = np.zeros((3, 3))

        for ij, (i, j) in enumerate(ut.upper_triangular(2)):

            aij0 = upper0[ij]
            aij1 = upper1[ij]
            aij2 = upper2[ij]

            a0[i, j] = aij0
            a0[j, i] = aij0

            a1[i, j] = aij1
            a1[j, i] = aij1

            a2[i, j] = aij2
            a2[j, i] = aij2

        a_new = np.zeros([3,3,3]) #will be returned
        a_new[0, :, :] = self.transform_alpha( a0, t1, t2, t3 )
        a_new[1, :, :] = self.transform_alpha( a1, t1, t2, t3 )
        a_new[2, :, :] = self.transform_alpha( a2, t1, t2, t3 )
        return a_new


    @staticmethod
    def transform_dist_beta( qm_beta, t1, t2, t3 ):
#Transform upper triangular to 3x3x3x form, rotate it, and transform back to ut style
        upper0 = np.array( qm_beta[0] )
        upper1 = np.array( qm_beta[1] )
        upper2 = np.array( qm_beta[2] )
        assert upper0.shape == (10,)
        assert upper1.shape == (10,)
        assert upper2.shape == (10,)
        b0 = np.zeros( (3, 3, 3))
        b1 = np.zeros( (3, 3, 3))
        b2 = np.zeros( (3, 3, 3))
        for ijk, (i, j, k) in enumerate(ut.upper_triangular(3)):
            bijk0 = upper0[ijk]
            bijk1 = upper1[ijk]
            bijk2 = upper2[ijk]

            b0[i, j, k] = bijk0
            b0[k, i, j] = bijk0
            b0[j, k, i] = bijk0
            b0[i, k, j] = bijk0
            b0[j, i, k] = bijk0
            b0[k, j, i] = bijk0

            b1[i, j, k] = bijk1
            b1[k, i, j] = bijk1
            b1[j, k, i] = bijk1
            b1[i, k, j] = bijk1
            b1[j, i, k] = bijk1
            b1[k, j, i] = bijk1

            b2[i, j, k] = bijk2
            b2[k, i, j] = bijk2
            b2[j, k, i] = bijk2
            b2[i, k, j] = bijk2
            b2[j, i, k] = bijk2
            b2[k, j, i] = bijk2

        b_new0 = np.zeros((3,3,3)) #will be returned
        b_new1 = np.zeros((3,3,3)) #will be returned
        b_new2 = np.zeros((3,3,3)) #will be returned
        b_new0 = self.transform_beta( b0, t1, t2, t3 )
        b_new1 = self.transform_beta( b1, t1, t2, t3 )
        b_new2 = self.transform_beta( b2, t1, t2, t3 )

        b0 = symmetrize_first_beta( b_new0 )
        b1 = symmetrize_first_beta( b_new1 )
        b2 = symmetrize_first_beta( b_new2 )

        return [ b0, b1, b2 ]

    def get_mol_string(self, basis = "cc-pVDZ" ):
        st = ""
        uni = Molecule.unique([ at.element for at in self])
        st += "ATOMBASIS\n\n\nAtomtypes=%d Charge=0 Nosymm\n" %(len(uni))
        for el in uni:
            print el
            st += "Charge=%s Atoms=%d Basis=%s\n" %( str(charge_dict[el]),
                    len( [all_el for all_el in self if (all_el.element == el)] ),
                    basis )
            for i in [all_el for all_el in self if (all_el.element == el) ]:
                st += "%s %.5f %.5f %.5f\n" %(i.element, i.x, i.y, i.z ) 
        return st

    def get_xyz_string(self):
        st = "%d\n\n" % len(self)
        for i in self:
            st += "{0:10s}{1:10f}{2:10f}{3:10f}\n".format(\
                    i.element, i.x,  i.y , i.z )
        return st

    @staticmethod
    def unique(arr):
        tmp = []
        for i in arr:
            if i not in tmp:
                tmp.append(i)
        return tmp

    @staticmethod
    def from_mol_file( molfile, AA = False):
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
    def mollist_to_mol_string( mollist , name ):
        print mollist
        raise SystemExit
    @staticmethod
    def from_xyz( f, in_AA = True, out_AA = True ):

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
        return m

class Water( Molecule ):
    """ Derives all general rotating methods from Molecule
    Specifics here for Water """

    def __init__(self , *args, **kwargs):
        super(Water, self).__init__( **kwargs )
        self.atoms = 0
        self.q = 0.0
        self.r_oh = False
        self.t_hoh = False

        self.euler1 = False
        self.euler2 = False
        self.euler3 = False

        self.no_hydrogens = True
        self.h1 = False
        self.h2 = False
        self.o  = False

        self.atomlist  = []

        self.AA = False
        self.Property = None

        self._coc = None

        self.in_qm = False
        self.in_mm = False
        self.in_qmmm = False

    def center(self):
        tmp = np.array( [0,0,0] )
        self.translate( tmp )

    @staticmethod
    def get_standard():
#Geometrical parameters
        AA = False
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

    @property
    def coc(self):
        if self._coc is not None:
            return _coc
        self._coc = sum( [at.r * charge_dict[at.element] for at in self])\
                /sum( map(float,[charge_dict[at.element] for at in self]) )
        return self._coc

    def append(self, atom):
        if len(self) > 3:
            print "tried to add additional atoms to water, exiting"
            raise SystemExit

        if not isinstance( atom, Atom ):
            print "wront class passed to water append"
            raise SystemExit

        if atom.element == "H":
            if self.no_hydrogens:
                self.h1 = atom
                atom.Water = self
                self.no_hydrogens = False
            else:
                self.h2 = atom
                atom.Water = self
        if atom.element == "O":
            self.o = atom
            atom.Water = self
            atom.name = "O1"
#Add the atom
        super( Water , self).append(atom)

#Define water center, by default set it to center of nuclei

        if (self.h1 and self.h2 and self.o):
            pass
            #self.r = np.array([at.r for at in self]).sum(axis=0) / 3.0

        if self.res_id:
            if self.res_id != atom.res_id:
                print "Tried to add %s to %s, exiting" %(atom, self)
                raise SystemExit
        else:
#Initialize water res_id from atomic res_id
            self.res_id = atom.res_id


#Fix for when arbitrary rotation, let hydrogen closests to upper octant be h1, the other h2
#The hydrogen closest to (1,1,1) gets to be h1, if they are equally close, then the one closest
#to the x axis is h1

#Also calculate center now
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

    def potline(self, max_l, pol, hyper, dist):
        return  "%d %.5f %.5f %.5f " %( 
                self.res_id, self.x, self.y, self.z ) + self.Property.potline( max_l, pol, hyper, dist ) + "\n"

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

    def get_dipole(self):
        hq = 0.25
        oq = -0.5
        return self.h1.get_array() * hq + self.h2.get_array() * hq + self.o.get_array() * oq

    def get_norm(self):
        r1 = self.h1 - self.o
        r2 = self.h2 - self.o
        return np.cross( r1, r2 )

    def dist_to_point( self , point ):
        return m.sqrt( (self.center[0] - point[0])**2 +\
                (self.center[1] - point[1])**2  + ( self.center[2] -point[2])**2 )

    def dist_to_water(self, other):
        xyz1 = self.center
        xyz2 = other.center
        return m.sqrt( (xyz1[0] - xyz2[0])**2 + \
            (xyz1[1] - xyz2[1])**2 + (xyz1[2] - xyz2[2])**2 )

    def set_property_on_each_atom(self):
        for i, prop in enumerate ( self.Property ):
            pass
        o_props =   { prop : self.Property[prop][0] for (key , prop)  in enumerate( self.Property )  }
        h1_props =  { prop : self.Property[prop][1] for (key , prop)  in enumerate( self.Property ) }
        h2_props =  { prop : self.Property[prop][2] for (key , prop)  in enumerate( self.Property ) }

        self.o.Property  = Property.from_template(  **o_props )
        self.h1.Property = Property.from_template( **h1_props )
        self.h2.Property = Property.from_template( **h2_props )

    def get_euler(self):
        """Return euler angles rho1, rho2, rho3 
        required to  water to its default placement
        for which the template properties are calculated """

        H1 = self.h1.r.copy()
        H2 = self.h2.r.copy()
        O1 = self.o.r.copy()

        dip = self.get_dipole()

        origin = O1.copy()
        H1, H2, O1 = H1 - origin, H2 - origin, O1 - origin

        theta1 = m.atan2( dip[1], dip[0])

        H1 =  np.dot( self.get_Rz_inv( theta1 ) , H1 )
        H2 =  np.dot( self.get_Rz_inv( theta1 ) , H2 )
        O1 =  np.dot( self.get_Rz_inv( theta1 ) , O1 )

        dip = np.dot( self.get_Rz_inv( theta1 ) , dip )

#Rotate by theta around y axis so that the dipole is in the z axis 
        theta2 = m.atan2( -dip[0], dip[2] )

        H1 =  np.dot( self.get_Ry( theta2 ) , H1 )
        H2 =  np.dot( self.get_Ry( theta2 ) , H2 )
        O1 =  np.dot( self.get_Ry( theta2 ) , O1 )

        dip = np.dot( self.get_Ry( theta2 ) , dip )

#Rotate around Z axis so that hydrogens are in xz plane.
        if H2[1] >0:
            xc = H2[0]
            yc = H2[1]
        else:
            xc = H1[0]
            yc = H1[1]
        theta3 = m.atan2( yc , xc)

        
        def eq(a, b, thr = 0.001): 
            if abs(a-b) < thr:return True
            else: return False

        #if eq( theta3, -np.pi ):
        #    theta3 = abs( theta3 )

        #if eq( theta2, -np.pi ):
        #    theta2 = abs( theta2 )

        #if eq( theta2, -np.pi ):
        #    theta1 = abs( theta1 )

        #if eq( theta1 , np.pi, ):
        #    theta1 = 0.0
        #if eq( theta1 , np.pi, ):
        #    theta1 = 0.0

        #if eq( theta1 , np.pi, ):
        #    theta1 = 0.0
        #if eq( theta1 , -np.pi, ):
        #    theta1 = 0.0

        return theta3, theta2, theta1

    def rotate(self, t1, t2, t3):
        """Rotate all coordinates by t1, t2 and t3
        first Rz with theta1, then Ry^-1 by theta2, then Rz with theta 3

        R all in radians

        """
        d1, d2, d3 = self.get_euler()

# Place water molecule in origo, and rotate it so hydrogens in xz plane
        H1 = self.h1.get_array() ; H2 = self.h2.get_array() ; O = self.o.get_array()
        TMP = self.o.get_array()
        H1 -= TMP ; H2 -= TMP; O -= TMP

        H1 = np.dot( self.get_Rz_inv(d3) , H1 )
        H1 = np.dot( self.get_Ry(d2) , H1 )
        H1 = np.dot( self.get_Rz_inv(d1) , H1 )

        H2 = np.dot( self.get_Rz_inv(d3) , H2 )
        H2 = np.dot( self.get_Ry(d2) , H2 )
        H2 = np.dot( self.get_Rz_inv(d1) , H2 )

        O = np.dot( self.get_Rz_inv(d3) , O )
        O = np.dot( self.get_Ry(d2) , O )
        O = np.dot( self.get_Rz_inv(d1) , O )

# Rotate with angles t1, t2, t3

        H1 = np.dot( self.get_Rz(t1) , H1 )
        H1 = np.dot( self.get_Ry_inv(t2) , H1 )
        H1 = np.dot( self.get_Rz(t3) , H1 )

        H2 = np.dot( self.get_Rz(t1) , H2 )
        H2 = np.dot( self.get_Ry_inv(t2) , H2 )
        H2 = np.dot( self.get_Rz(t3) , H2 )

        O = np.dot( self.get_Rz(t1) , O )
        O = np.dot( self.get_Ry_inv(t2) , O )
        O = np.dot( self.get_Rz(t3) , O )

#Put back in oxygen original point
        H1 += TMP ; H2 += TMP; O += TMP

        self.h1.x = H1[0] ;self.h1.y = H1[1] ;self.h1.z = H1[2] 
        self.h2.x = H2[0] ;self.h2.y = H2[1] ;self.h2.z = H2[2] 
        self.o.x  =  O[0] ;  self.o.y = O[1] ;  self.o.z = O[2] 

    def to_au(self):
        self.h1.to_au()
        self.h2.to_au()
        self.o.to_au()

    def to_AA(self):
        self.h1.to_aa()
        self.h2.to_aa()
        self.o.to_aa()
# in_AA specifies if input coords are in angstrom

    def get_xyz(self):
        st = "%d\n\n" % len(self)
        for at in self:
            st += "{0:10s}{1:10f}{2:10f}{3:10f}\n".format(\
                    at.element, at.x,  at.y , at.z )
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
                atom.res_id = wat.res_id
        if not out_AA:
            for wat in waters:
                wat.to_au()
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
                string +=  "%d %.5f %.5f %.5f " % tuple(
                        [int(at.res_id)] + at.r)
                string += at.Property.potline( max_l=max_l, pol=pol, hyper= hyper)
                string += '\n'
        return string




class Methanol(Molecule):

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
    def __init__(self, *args, **kwargs):
        """ Typical list of molecules """
        pass

    def __str__(self):
        return " ".join( [ str(i) for i in self ] )

    def append(self, mol, in_mm = False, in_qm = False,
            in_qmmm = False):
        mol.in_mm = in_mm
        mol.in_qm = in_qm
        mol.in_qmmm = in_qmmm

        super( Cluster, self ).append( mol )

    def min_dist(self):
        dist = np.zeros( len(self) )
        for i in range(len(self)):
            for j in range(i ,len(self)):
                if i == j:
                    continue
                dist[i] = ( np.linalg.norm(self[i].center - self[j].center) )
        dist.sort()
        return dist

    def get_qm_mol_string(self, basis = ("cc-pVDZ", ) , AA = False):
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
    def get_qmmm_pot_string( self, max_l = 1, pol = 2, hyp = 0, AA = False ):
        if AA:
            st = "AA\n"
        else:
            st = "AU\n"
# Old qmmm format requires integer at end to properly read charges
        st += "%d %d %d %d\n" % (sum([len(i) for i in self if i.in_mm ]), 
                max_l, pol, 1 )
        st += "".join( [at.potline(max_l, pol, hyp) for mol in self for at in mol if mol.in_mm] )
        return st

    def get_xyz_string(self):
        st = "%d\n\n" % sum([len(i) for i in self ])
        for mol in self:
            for i in mol:
                st += "{0:10s}{1:10f}{2:10f}{3:10f}\n".format(\
                        i.element, i.x,  i.y , i.z )
        return st


    def order_mm_atoms(self):
        cnt = 1
        for mol in [m for m in self if m.in_mm]:
            for at in mol:
                at.number = str(cnt)
                cnt += 1
    @staticmethod
    def get_water_cluster( fname , in_AA = True, out_AA = True , N_waters = 1):
        """From file with name fname, return a Cluster with all waters encountered"""
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
                    if ( i[11:16].strip() == "SW") or (i[11:16] == "DW"):
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
        for wat in c:
            for atom in wat:
                atom.res_id = wat.res_id
        if not out_AA:
            for wat in c:
                wat.to_au()
        return c


    def mol_too_close(self, mol):
        for mols in self:
            for ats in mols:
                for at in mol:
                    if at.dist_to_atom( ats ) < 2.4:
                        return True
        return False


    def add_atom(self, at):
        self.append( at )
        at.cluster = self
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

if __name__ == '__main__':

# Water bonding parameters:
    r_oh = 0.97167 ; theta_hoh = 104.5
# Water rotation parameters
    r = 5.0
    #theta = 45 ; tau = 45
    euler1 = 0.0 ; euler2 = 0.0 ; euler3 = 0.0
    p = m.pi ; c = m.cos ; s = m.sin
    euler1 *= p/180 ;euler2 *= p/180 ;euler3 *= p/180 ;
    theta_hoh *= p/180

#Hardcoded conversions
    for i in np.r_[0 : 180: 8j ]:
        for j in np.r_[0 : 360 : 16j ]:

            theta = i * p/180 
            tau = j * p/180 
            x = r * s(theta) * c(tau )
            y = r * s(theta) * s(tau )
            z = r * c(theta)
            w1 = g.get_water( [ 0, 0, 0], r_oh, theta_hoh )
            w1.res_id = 1 ; w1.r = r ; w1.theta = theta ; w1.tau = tau
            w2 = g.get_water( [ x, y, z], r_oh, theta_hoh )
            w2.res_id = 2 ; w2.r = r ; w2.theta = theta ; w2.tau = tau
            g.writeMol( [ w1, w2 ])
    


