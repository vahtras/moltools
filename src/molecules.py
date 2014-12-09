#!/usr/bin/env python
#-*- coding: utf-8 -*-

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt

import itertools

import numpy as np
import re, os

a0 = 0.52917721092

charge_dict = {"H": 1.0, "C": 6.0, "N": 7.0, "O": 8.0, "S": 16.0}
# from TIP3P charge defs.
el_charge_dict = {"H": .417, "O": -0.834 }
mass_dict = {"H": 1.008,  "C": 12.0, "N": 14.01, "O": 15.999, "S": 32.066}

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

class Property( dict ):
    def __init__(self):

        self["charge"] = np.zeros( 1 )
        self["dipole"] = np.zeros( 3 )
        self["quadrupole"] = np.zeros( 6 )
        self["alpha"] =  np.zeros( 6 ) 
        self["beta"] =  np.zeros( 10 ) 

    def __add__(self, other):
        assert isinstance( other, Property)
        tmp = {}
        for i, prop in enumerate(self):
            tmp[prop] = np.array( self[prop] ) + np.array(other[prop] )
        return tmp
    def __sub__(self, other):
        assert isinstance( other, Property)
        tmp = {}
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
        p = Property()
        for i, keys in enumerate( wat_templ ):
            if keys[0] == at.name:
                p[keys[1]] = wat_templ[ keys ]
        at.Property = p

    def transform_ut_properties( self, t1, t2, t3):
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
    def __init__(self):
        """Container class for all rotation related operations"""
        pass

    @staticmethod
    def transform_1( qm_dipole, t1, t2, t3 ):
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

    def copy_atom(self):
        a = Atom()

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
    def __add__(self, other ):
        return self.r + other.r

    def get_array(self):
        return np.array( self.r ).copy()

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
        self.cluster = None

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
    @property
    def com(self):
        return np.array([at.mass*at.r for at in self]).sum(axis=0) / np.array([at.mass for at in self]).sum()

    @staticmethod
    def dist_to_mol(self, other):
        xyz1 = self.com
        xyz2 = other.com
        return m.sqrt( (xyz1[0] - xyz2[0])**2 + \
            (xyz1[1] - xyz2[1])**2 + (xyz1[2] - xyz2[2])**2 )

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

        for ij, (i, j) in enumerate(upper_triangular(2)):

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
        a_new[0, :, :] = self.transform_2( a0, t1, t2, t3 )
        a_new[1, :, :] = self.transform_2( a1, t1, t2, t3 )
        a_new[2, :, :] = self.transform_2( a2, t1, t2, t3 )
        return a_new

    def get_mol_string(self, basis = ("ano-1 2 1", "ano-1 3 2 1" ) ):
        if len( basis ) > 1:
            el_to_rowind = {"H" : 0, "C" : 1, "O" : 1, "N" : 1  }
        else:
            el_to_rowind = {"H" : 0, "C" : 0, "O" : 0, "N" : 0 }
        st = ""
        s_ = ""
        if self.AA: s_ += "Angstrom"
        uni = Molecule.unique([ at.element for at in self])
        st += "ATOMBASIS\n\n\nAtomtypes=%d Charge=0 Nosymm %s\n" %(len(uni), s_)
        for el in uni:
            st += "Charge=%s Atoms=%d Basis=%s\n" %( str(charge_dict[el]),
                    len( [all_el for all_el in self if (all_el.element == el)] ),
                    basis[ el_to_rowind[el] ])
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

        if in_AA:
            if not out_AA:
                m.to_AU()
        return m

    def to_AU(self):
        for at in self:
            at.x = at.x / a0
            at.y = at.y / a0
            at.z = at.z / a0
    def to_AA(self):
        for at in self:
            at.x *= a0
            at.y *= a0
            at.z *= a0

class Water( Molecule ):
    """ Derives all general methods from Molecule.
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

    def copy_water(self):
        w = Water()
        [w.append(i.copy_atom()) for i in self]

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

#Center of oxygen
    @property
    def coo(self):
        return self.o.r

#Center of charge
    @property
    def coc(self):
        return sum( [at.r * charge_dict[at.element] for at in self])\
                /sum( map(float,[charge_dict[at.element] for at in self]) )

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
# Not accurate magnitude, only direction of dipole vector
        hq = 0.25
        oq = -0.5
        return self.h1.get_array() * hq + self.h2.get_array() * hq + self.o.get_array() * oq

    def get_norm(self):
        r1 = self.h1 - self.o
        r2 = self.h2 - self.o
        return np.cross( r1, r2 )

    def dist_to_point( self , point ):
        return np.sqrt(np.sum((self.coo - np.array(point))**2))

    def dist_to_water(self, other):
        return np.sqrt(np.sum((self.coo - other.coo)**2) )

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

    def update_water_props(self, model = "TIP3P",
            method = "HF", basis = "ANOPVDZ", dist = False,
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
                atom.res_id = wat.res_id
        if in_AA:
            if not out_AA:
                for wat in c:
                    wat.to_AU()
        return c


    def mol_too_close(self, mol):
        for mols in self:
            for ats in mols:
                for at in mol:
                    if at.dist_to_atom( ats ) < 2.0:
                        return True
        return False

    def add_mol(self, at):
        self.append( mol )
        mol.cluster = self

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
    def copy_cluster(self):
        tmp_c = Cluster()
        [tmp_c.append(wat.copy_water()) for wat in self]
        return tmp_c

if __name__ == '__main__':

# Water bonding parameters:
#
    fo = "hfqua_tip3p11_9qm.out"
    fm = "tip3p11_9qm.mol"
    import read_dal
    from gaussian import GaussianQuadrupole, GaussianQuadrupoleList

    at, p, a, b = read_dal.read_beta_hf( fo )
    c = Cluster.get_water_cluster( fm, in_AA = False, out_AA = False, N_waters = 100)
    c.update_water_props( dist = False )
    static_ox = GaussianQuadrupoleList.from_string( Water.get_string_from_waters( c, pol= 22, hyper = 1 ) )
    c.update_water_props( dist = True )
    static_dist = GaussianQuadrupoleList.from_string( Water.get_string_from_waters( c, pol= 22, hyper = 1, dist = True ) )

    print static_ox.total_dipole_moment()
    print static_dist.total_dipole_moment()
    print p



    raise SystemExit
    for i in range(101):
        for j in range(2,10):
            fo = "hfqua_tip3p%d_%dqm.out" %(i,j)
            fm = "tip3p%d_%dqm.mol" %(i,j)
            if not os.path.isfile( os.path.join(os.getcwd(),fo) ):
                continue

            at, p, a, b = read_dal.read_beta_hf( fo )
            c = Cluster.get_water_cluster( fm, in_AA = False, out_AA = False, N_waters = 100)
            c.update_water_props( dist = False )
            static_ox = GaussianQuadrupoleList.from_string( Water.get_string_from_waters( c, pol= 22, hyper = 1 ) )
            c.update_water_props( dist = True )
            static_dist = GaussianQuadrupoleList.from_string( Water.get_string_from_waters( c, pol= 22, hyper = 1, dist = True ) )
            try:
                np.testing.assert_almost_equal( np.linalg.norm(static_ox.total_dipole_moment()) , np.linalg.norm( static_dist.total_dipole_moment()),decimal = 1)
            except:
                print i
                #print np.linalg.norm( static_dist.total_dipole_moment() )
                #raise SystemExit



