#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""
The molecules modules serves as an interface to write water molecule input files using predefined geometries, to be used with the DALTON qm package.
"""

__all__ = [ 'Atom', 'Bond', 'Molecule', 'Water', 'Cluster', ]

import re, os, itertools, functools, warnings, subprocess, shutil, logging, string, tarfile
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt

import numpy as np
import cPickle as pickle
import copy as copymod


from .utilz import scale_vec_to_abs, Rz, Ry, Rz_inv, Ry_inv, au_to_nm, center_and_xz, get_rotation, s2ut, ut2s, reflect, unique, find_dir
from .property import Property
from .template import Template
from .generator import Generator

from loprop.core import MolFrag, penalty_function, shift_function
from applequist import gaussian

a0 = 0.52917721092
au_nm_conv = 45.563352491
elem_array = ['X', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'P']

charge_dict = {"H": 1.0, "He" :2.0, "Li" : 3.0,  "C": 6.0,
        "N": 7.0, "O": 8.0, "S": 16.0, "P" : 15.0, "X" : 0.0 }
# from TIP3P charge defs.
el_charge_dict = {"H": .417, "O": -0.834 , "X" : 0.417 , 'S': -0.25}
mass_dict = {"H": 1.008,  "C": 12.0, "N": 14.01, "O": 15.999, "S": 32.066,
        "X" : 1.008, 'P' : 30.974 }

color_dict = { "X": 'black' ,"H":'brown', "N":'blue',"C":'green',"P":'black', "O":'red', 'S' : 'yellow'}

bonding_cutoff = { 
        ('X','X') : 0.0,
        ('H','X') : 0.0,
        ('H','H') : 0.2,
        ('H','C') : 1.155,
        ('H','N') : 1.155,
        ('H','O') : 1.1,
        ('H','P') : 1.1,
        ('H','S') : 1.3,
        ('C','X') : 0.0,
        ('C','C') : 1.66,
        ('C','N') : 1.60,
        ('C','O') : 1.5,
        ('C','P') : 2.0,
        ('C','S') : 2.0,
        ('N','X') : 0.0,
        ('N','N') : 1.5,
        ('N','O') : 1.5,
        ('N','P') : 1.5,
        ('N','S') : 1.5,
        ('O','X') : 0.0,
        ('O','O') : 1.5,
        ('O','P') : 2.0,
        ('O','S') : 1.75,
        ('P','X') : 0.0,
        ('P','P') : 1.5,
        ('P','S') : 2.0,
        ('S','X') : 0.0,
        ('S','S') : 2.1,
    }

res_dict = {'ALA':'A', 'VAL':'V', 'ILE':'I','LEU':'L','MET':'M',
        'PHE':'F','TYR':'Y','TRP':'W','SER':'S','THR':'T','ASN':'N', 'CRO':'X1',
        'CRO1':'X1', 'CRO2':"X2",'CRO3':'X3','CRO4':'X4',
        'GLN':'Q','CYS':'C','CH6': 'X1' ,'GLY':'G','PRO':'P','ARG':'R','HIS':'H',
        'LYS':'K','ASP':'D','GLU':'E','SEC':'U','PYL':'U', 'HIP':'Z',
        'HIE':'H','CYX':'C','HSE':'H','HID':'H','HSD':'H','HSP':'H2',"TIP3": 'T3',
        'HIP':'H2','HYP':'PX', 'MOL' : 'X', 'WAT' : 'W1', 'SOL' : 'W1' }
chargeDict = {'ARG':1, "LYS" : 1, "ASP":-1, "GLU":-1,
            'R':1 , 'K':1 ,'H':0, 'H2':1 , 'E':-1 , 'D':-1, 'X2': 1}
proline_dict = { "PRO" : "P", "HYP" : "PX" }
custom_dict = { "CRO2" : "X2", "MOL" : "X3" }

#Make permutations of all bonding pair tuples 
for key1, key2 in bonding_cutoff.keys():
    bonding_cutoff[ (key2, key1)] = bonding_cutoff[ (key1, key2) ]

class UnitException( Exception ):
    def __init__( self, unit, label ):
        """unit is bool, labeltom is string label of atom"""
        self.unit = unit
        self.label = label
        
class Bond(object):
    """Docstring for Bond. """
    def __init__(self, at1, at2 ):
        """ first argument is _Atom1, second argument becomes _Atom2 """
        self._r = None
        self._Property = None
        self._label = None
        self._res_id = None
        self._Atom1 = at1
        self._Atom2 = at2
        self._Molecule = None
        self._r = (at2.r - at1.r)/2.0 + at1.r
        self.element = 'X'

        if at1.res_id == at2.res_id:
            self._Molecule = at1.Molecule

    @property
    def res_id(self):
        if self._res_id is not None:
            return self._res_id
        if self._Molecule is not None:
            return self._Molecule.res_id
        if self._Atom1.res_id == self._Atom2.res_id:
            return self._Atom1.res_id
        return None

    @property
    def len(self):
        return self._Atom1.dist_to_atom( self._Atom2 )

    def get_mol_line(self, lab = None):
        if lab is None:
            lab = 'X-%s-%s' %(self._Atom1.pdb_name, self._Atom2.pdb_name, )
        return "{0:15s}{1:10.5f}{2:10.5f}{3:10.5f}\n".format( lab, self.r[0], self.r[1], self.r[2] ) 
    @property
    def r(self):
        if self._r is not None:
            return self._r
    @r.setter
    def r(self, val):
        assert len(val) == 3
        self._r = val

    @property
    def p(self):
        if self._Property is not None:
            return self._Property
        return Property()
    @p.setter
    def p(self, val):
        self._Property = val

#Method of Bond
    def potline(self, max_l=2, pol=22, hyper=1, fmt = "%.5f ",):
        return  "{0:4} {1:10f} {2:10f} {3:10f} ".format( \
                str(self._Atom1.Molecule.cluster_order), *self.r ) + self.p.potline( max_l, pol, hyper, fmt = fmt ) + "\n"

    @property
    def label(self):
        if self._label is None:
            _str = "-".join( map(str, [ self._Atom1.res_id,
                self._Atom2.res_id,
                self._Atom1.pdb_name,
                self._Atom2.pdb_name,
                ] 
                ))
            return _str
        return self._label

    def __str__(self):
        return self.label + " " + " ".join(map(str,self.r))

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

#order in xyz files
        self._order = None

#name is custom name, for water use O1, H2 (positive x ax), H3
        self._name = None

#label is custom name, to identify from different residues
        self._label = None

#cabel is custom name to identify from different chains
        self._clabel = None

        self.x = 0.0
        self.y = 0.0
        self.z = 0.0


# Use populate_bonds in class Molecule to attach all atoms to their neighbours
#bonds is list of class Bond
        self._bonds = None

        self.angles = {}
        self.dihedrals = {}

        self._q = None


        self._number = None
        self._pdb_name = None
        self._res_id = 0
        self._atom_id = None
        self._chain_id = None

        self.in_water = False

#Connectivity
        self._Molecule = None
        self._Cluster = None
        self._System = None
        self._World = None

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
            self.order = kwargs.get( "order", 0 )
            self.in_qm = kwargs.get( "in_qm", False )
            self.in_mm = kwargs.get( "in_mm", False )
            self.in_qmmm = kwargs.get( "in_qmmm", False )
            self._pdb_name = kwargs.get( "pdb_name", None )
            self._res_id = kwargs.get( "res_id", None )
            self._chain_id = kwargs.get( "chain_id", None )
        self._mass = None
#Check if atom has same coordinate and element as other
    def equal(self, other):
        if np.allclose( self.r, other.r, atol = 1e-5 ) and self.element == other.element:
            return True
        return False

    def __eq__(self, other):
        """Just compare coordiante of atom and call it the same if matches,
        skip particular pointers to atoms"""
        return np.allclose( self.r, other.r ) 

    def is_dummy(self):
        if self.label:
            l = self.label.split('-')
            if len( l ) == 4:
                if l[-1][0] == 'X':
                    return True
        return False

    @property
    def number( self ):
        if self._number is not None:
            return self._number
        return None

    @number.setter
    def number( self, val ):
        self._number = val

    @property
    def bonds( self ):
        if self._bonds is not None:
            return self._bonds
        self._bonds = []
        #Can add other stuff here later
        return self._bonds

    @bonds.setter
    def bonds( self, val ):
        self._bonds  = val

#Add a bond for this atom to the other atom
    def add_bond( self, b ):
        if any( np.allclose(b.r, x.r, atol = 1e-7 ) for x in self.bonds):
            #logging.warning( "Tried to add bond which already exsists in %s" %self.pdb_name )
            return
        self.bonds.append( b )

    def t(self, r):
        """translate by vector r"""
        self.x = self.x + r[0]
        self.y = self.y + r[1]
        self.z = self.z + r[2]

#Chain for this atom
    @property
    def Chain(self):
        return self.Cluster

#Chain for this atom
    @property
    def clabel(self):
        if self._clabel:
            return self._clabel
        m = self.Molecule
        if m:
            c = m.Cluster
            if c:
                return "{}-{}-{}-{}".format( c.chain_id, m.res_id,m.res_name,self.pdb_name )
        return None

#Molecule for this atom
    @property
    def Molecule(self):
        if self._Molecule:
            return self._Molecule
        return None

#Molecule for this atom
    @Molecule.setter
    def Molecule(self, val):
        if isinstance( val, Molecule ):
            self._Molecule = val 
        else:
            logging.warning("Tried to set wrong type to atoms _Molecule")

#Cluster for this atom
    @property
    def Cluster(self):
        if self._Cluster:
            return self._Cluster
        elif self.Molecule:
            return self.Molecule.Cluster
        return None

#Molecule for this atom
    @Cluster.setter
    def Cluster(self, val):
        if isinstance( val, Cluster ):
            self._Cluster = val 

#Chain ID property of Atom
    @property
    def chain_id(self):
        if self._chain_id:
            return self._chain_id
        if self.Molecule:
            if self.Molecule.Cluster:
                return self.Molecule.Cluster.chain_id
        return 'X'
    @chain_id.setter
    def chain_id(self, val ):
        self._chain_id = val



# property setters and getters for pdb_name
    @property
    def pdb_name(self):
        if self._pdb_name is None:
            return '%s' %(self.element)
        return self._pdb_name
    @pdb_name.setter
    def pdb_name(self, val):
        self._pdb_name = val

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
        if self.Molecule and self._element:
            return self._element + str(self.order)
        return 'X'
    @name.setter
    def name(self,val):
        self._name = val
    @property
    def label(self):
        if self._label is not None:
            return self._label
        if self.Molecule:
            m = self.Molecule
            return '-'.join( [ str(m.res_id), m.res_name, self.pdb_name] )
        return None
    @label.setter
    def label(self, val):
        self._label = val

#Atoms
    @property
    def order(self):
        if self._order:
            return self._order
        if self.Molecule:
            return self.Molecule.index(self) + 1
        return 0
    @order.setter
    def order(self, val):
        self._order = val

    def disconnect(self, other):
        if other.name in self.bonds:
            del self.bonds[ other.name ]

#Atoms
    def transfer_props(self, level = 1, own_res = True ):
        """
        Transfer own props completely based on level.

        level = 1: Transfer to bonds evenly
        level = 2: Transfer own, and bonds, to nearest bonded atoms.
        """
        bonds = [b for b in self.bonds if b._Atom2.res_id == self.res_id]
        sites = len( bonds )
        p = Property()
        if level == 1:
            for own_bond in bonds:
                try:
                    other_bond = [b for b in own_bond._Atom2.bonds if b._Atom2 == self ][0]
#The other bond has no bond to this one! its fine in some cases
                except IndexError:
                    continue
                own_bond.p += self.p/sites
                other_bond.p += self.p/sites
                p += self.p/sites

        elif level == 2:
            for own_bond in bonds:
                try:
                    other_bond = [b for b in own_bond._Atom2.bonds if b._Atom2 == self ][0]
                except IndexError:
                    print 'no bonds in other atom'
                    continue
#If for some reason we forgot to put properties on this bond but not for the other

                own_bond._Atom2.p += own_bond.p
                own_bond._Atom2.p += self.p / sites
                own_bond.p = Property()
        self.p = Property()
        

    @property
    def com(self):
        return self.r
    def get_mol_line(self):
        return "{0:15s}{1:10.5f}{2:10.5f}{3:10.5f}\n".format( self.label, self.x, self.y, self.z ) 

    def in_mm(self):
        return self.Molecule.in_mm
    def in_qm(self):
        return self.Molecule.in_qm

    def xyz_string(self):
        x, y, z = map( lambda x: string.rjust( "%.3f"%x, 8)[:8], self.r )
        return "{0:10s} {1:10f} {2:10f} {3:10f}\n".format(self.element, self.x,  self.y, self.z )

    def pdb_string(self):
        x, y, z = map( lambda x: string.rjust( "%.3f"%x, 8)[:8], self.r )
        """Return pdb line as specified by the PDB standard"""
        st = "{0:6s}{1:5s}{2:1s}{3:4s}{4:1s}{5:4s}{6:1s}{7:4s}{8:1s}{9:3s}{10:8s}{11:8s}{12:8s}{13:22s}{14:2s}\n".format( "ATOM", 
                str(self.atom_id),
                "",
                self.pdb_name,
                "",
                self.Molecule.res_name,
                self.chain_id,
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
        if self.Molecule:
            self._atom_id = self.Molecule.index( self )
        return None

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
        new = scale_vec_to_abs( vec, value = final_value )
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
    def get_async_bond(self, alist = [], first = None ):
        """By specifying the first atom, will recursively
        grab all atoms in that direction until the tree stops
        
        useful when rotating coordinates around the cross vector of an angle
        in order to rotate whole molecule about that angle"""

        if first:
            self.Molecule.populate_bonds()
            alist = []
        if first is not None:
            alist.append( self )
            return first.get_async_bond( alist = alist )
        else:
            visit = []
            for a in self.bonds.values():
                if a in alist:
                    continue
                visit.append( a )
            alist.append( self )
            if visit == []:
                return alist
            for a in visit:
                alist = a.get_async_bond( alist = alist )
            return alist
            

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
        Rz, Rzi = Rz, Rz_inv
        Ry, Ryi = Ry, Ry_inv

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
            self.x, self.y, self.z = self.bonds[at].r + (self.r - self.bonds[at].r)*scale

#Atom method
    def copy(self):
        return self.copy_atom()
    def copy_self(self):
        return self.copy_atom()

    def copy_atom(self):
#type(self) returns this particular atom type, which can be inherited
        a = type(self)(**{'x':self.x, 'y':self.y, 'z':self.z,'AA':self.AA,
            'element':self.element,'name':self.name,'number':self.number,
            'pdb_name':self.pdb_name} )

        for key, val in self.__dict__.iteritems():
#Do not copy pointers
            if isinstance( val, list ):
                continue
            if isinstance( val, Property ):
                a.p = self.p.copy_property()
            a.__dict__[ key ] = val

        return a

    @property
    def r(self):
        return np.array( [ self.x, self.y, self.z ] )
    @r.setter
    def r(self, val):
        assert len(val) == 3
        self.x, self.y, self.z = val[0], val[1], val[2]

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

#Atom method
    @property
    def res_id(self):
        if self._res_id:
            return self._res_id
        if self.Molecule:
            return self.Molecule.res_id
        return 0

#Method of Atom
    def potline(self, max_l=2, pol=22, hyper=1, fmt = "%.5f ",):
        return  "{0:4} {1:10f} {2:10f} {3:10f} ".format( \
                str(self.Molecule.cluster_order), self.x, self.y, self.z ) + self.Property.potline( max_l = max_l, pol = pol, hyper = hyper , fmt = fmt ) + "\n"

#Atom string method
    def __str__(self):
        return "%s %f %f %f" %(self.label, self.x, self.y, self.z)

    #def __repr__(self):
    #    st = '"A%d' %self.order
    #    if self.Molecule:
    #        st = st.rstrip('"') +  '-M%d"' %self.Molecule.cluster_order
    #        if self.Molecule.Cluster:
    #            st = st.rstrip('"') + '-C%d"' %self.Molecule.Cluster.system_order
    #            if self.Molecule.Cluster.System:
    #                st.rstrip('"')
    #                st += '-S"'
    #    return st

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
        return np.sqrt( ((self.r - other.r)**2).sum() )



    def dist_to_point(self, other):
        """
Return the distance to a point

.. code:: python

   >>> a = Atom( z = 0 )
   >>> print H1.dist_to_point( [0, 3, 4] )
   5.0

"""
        return np.sqrt( ((self.r - other)**2).sum() )

    def to_AU(self):
        if self.AA:
            self.x /= a0
            self.y /= a0
            self.z /= a0
            for b in self.bonds:
                b.r /= a0
            self.AA = False

    def to_AA(self):
        if not self.AA:
            self.x *= a0
            self.y *= a0
            self.z *= a0
            for b in self.bonds:
                b.r *= a0
            self.AA = True

class Molecule( list ):
    """
**Inherits list methods, specific molecules will inherit from this class.**
"""

    def __init__(self , *args, **kwargs):
        super( Molecule, self).__init__()

# Dictionary with bonds
        self.bond_dict = {}

#center will be defined for all molecules after all atoms are added
#depends on which molecule
        self._chain_id = None
        self._res_name = None
        self._freq = None
        self._res_id = None
        self._r = None
        self._com = None
        self._Cluster = None
        self._order  = None
        self._cluster_order  = None
        self.no_hydrogens = True
        self.is_concap = False
        self.is_ready = False
        self.is_bridge = False
        self.n_term = False
        self.c_term = False
        self._level = None


# This will be set True if Property is represented by a point on molecule
        self._is_Property = False
        self._property_r = None

# This will be set True if attaching LoProp properties
        self._LoProp = False

# For plotting different elements:
        self.style = { "X": 'ko' ,"H":'wo', "N":'bo',"C":'go',"P":'ko', "O":'ro',
                'S' : 'yo'}
        self.linewidth = {"X":25,"H":25, "N": 30, "C": 30, "O":40, "P" : 40,
                'S' : 45 }

# Make empty, beware that this causes molecules to give zero dipole momnet
# before template is loaded
        self._Property = None

#Access snapshot of cluster if exists
        self._snap = None

#By default, in no region
        self._in_qm = False
        self._in_mm = False
        self._in_qmmm = False

#By default, AU 
        #self.AA = False

#if supplied a dictionary with options, gather these in self.info
        self.info = {}

        if type(args) == tuple:
            if len(args) == 1:
                if type(args[0]) == list:
                    for i in args[0]:
                        self.initiate( i )
                else:
                    self.initiate( args[0] )
            else:
                for at in args:
                    self.initiate( at )

        if kwargs != {} :
            for i in kwargs:
                self.info[ i ] = kwargs[ i ]
            self.AA = kwargs.get( "AA" , False )

    def atom_in_other_mol_closer_than( self, other, dist = 1.0, AA=True ):
        """Given other molecule, return True if other molecule
        has an Atom with distance less than 1.0 to this molecule"""
        if (not self.AA) and AA:
            dist /= a0
        for i in self:
            for j in other:
                if i.dist_to_atom( j ) < dist:
                    return True
        return False


    @property
    def nuc_charge(self ):
        return sum( [ charge_dict[ i.element ] for i in self] )

#To return all atoms and bonds in molecule
    def get_ats_and_bonds(self):
        """Important not to overwrite bonds which has properties, just return
        
        """
        tot = []
        bond_visited = []
        for at in self:
            tot.append( at )
            for b in at.bonds:
                if any( np.allclose(b.r, x.r, atol = 1e-7) for x in bond_visited ):
                    continue
                bond_visited.append( b )
                tot.append( b )
        return tot

#Unique identifier which will produce file name string unique to this residue
    def file_label( self, freq = True ):
        if freq:
            f = au_to_nm( self.freq )
            return "_".join( map(str, [self.res_name, self.res_id, f]) )
        else:
            return "_".join( map(str, [self.res_name, self.res_id ]) )



#Unique identifier which will produce file name string unique to this residue
    def molfile_label( self, freq = True ):
        if freq:
            f = au_to_nm( self.freq )
            return "_".join( map(str, [self.res_name, self.res_id, f]) )
        else:
            return "_".join( map(str, [self.res_name, self.res_id ]) )


    @property
    def bonds(self):
        b = reduce(lambda a, x: a + x, [a.bonds for a in self], [] )
        return b

#Molecule chain_id
    @property
    def chain_id(self):
        if self._chain_id is not None:
            return self._chain_id
        if self.Cluster:
            return self.Cluster.chain_id
        tmp_ch = self[0].chain_id
        for at in self:
            try:
                assert tmp_ch == at.chain_id
            except AssertionError:
                logging.error( "No Chain object or _chain_id in NewResidue and not all atoms have same chain_id")
        return tmp_ch

#Molecule
    @property
    def Cluster(self):
        if self._Cluster:
            return self._Cluster
        return None
#Molecule
    @Cluster.setter
    def Cluster(self, val):
        if isinstance( val, Cluster ):
            self._Cluster = val 

#Molecule method to transfer props from all atoms in at_list evenly to neighbours
#In order to remove props to closest bonded neighbours
    def transfer_props(self, at_list, 
            transfer = { 'charge' : 0,
                'quadrupole' : 0,
                'dipole' : 0,
                'alpha' : 0,
                'beta' : 0 }, 
            center = None ):
        """
        The algorithm checks for edges, takes properties from there 
        if they are also in the list, and then
        distributes evenly to other neighbours
        """

        if center is not None:
            if center.shape == (3,):
                at_list = sorted( at_list, key = lambda x: (len(x.bonds), x.dist_to_point(center),) )
            else:
                logging.error('Provided center keyword in transfer props\
                        without proper numpy array of shape 3!')
        else:
            at_list = sorted( at_list, key = lambda x: len(x.bonds) )

#atoms and their neighbours should have a tmp copy of bonds
        for each in at_list:
            each.tmp_bonds = each.bonds.copy()
            for at in each.bonds.values():
                at.tmp_bonds = at.bonds.copy()

#
        props = [k for k, v in transfer.iteritems() if v == 1]
        for each in at_list:
            if len( each.tmp_bonds ) == 0:
                logging.error('Tried to transfer props from non bonded atom')
                logging.error('Offending atom is : %s' %each )
                raise SystemExit
            p = each.p / len( each.tmp_bonds )
            for at in each.tmp_bonds.values():
                for prop in props:
                    at.p[prop] += p[prop]
                    each.p[prop] -= p[prop]
                del at.tmp_bonds[ each.name ]
            
#So the side groups dont transfer props back to this atom


#Properties for grabbing uniquely PDB named atoms
    @property
    def N(self):
        return self.get_atom_by_pdbname( 'N' )
    @property
    def CA(self):
        return self.get_atom_by_pdbname( 'CA' )
    @property
    def HA(self):
        return self.get_atom_by_pdbname( 'HA' )
    @property
    def HA1(self):
        return self.get_atom_by_pdbname( 'HA1' )
    @property
    def HA2(self):
        return self.get_atom_by_pdbname( 'HA2' )
    @property
    def C(self):
        return self.get_atom_by_pdbname( 'C' )
    @property
    def O(self):
        return self.get_atom_by_pdbname( 'O' )
    @property
    def OG(self):
        return self.get_atom_by_pdbname( 'OG' )
    @property
    def OG2(self):
        return self.get_atom_by_pdbname( 'OG2' )
    @property
    def CA(self):
        return self.get_atom_by_pdbname( 'CA' )
    @property
    def CB(self):
        return self.get_atom_by_pdbname( 'CB' )
    @property
    def CB2(self):
        return self.get_atom_by_pdbname( 'CB2' )
    @property
    def CD(self):
        return self.get_atom_by_pdbname( 'CD' )
    @property
    def HD1(self):
        return self.get_atom_by_pdbname( 'HD1' )
    @property
    def HD2(self):
        return self.get_atom_by_pdbname( 'HD2' )
    @property
    def CG(self):
        return self.get_atom_by_pdbname( 'CG' )
    @property
    def H(self):
        return self.get_atom_by_pdbname( 'H' )
    @property
    def HN(self):
        return self.get_atom_by_pdbname( 'HN' )


    @property
    def X1(self):
        return self.get_atom_by_pdbname( 'X1' )

    @property
    def LoProp(self):
        return self._LoProp
    @LoProp.setter
    def LoProp(self, val):
        assert type(val) is bool or val is None
        if val:
            self.is_Property = False
            self.property_r = np.zeros( 3 )
        self._LoProp = val

#Properties relating to whether molecule has fixed point for Properties instead of LoProp
    @property
    def is_Property(self):
        return self._is_Property
    @is_Property.setter
    def is_Property(self, val):
        assert type(val) is bool or val is None
        if val:
            self._is_Property = val
            self.LoProp = False
    @property
    def property_r(self):
        return self._property_r
    @property_r.setter
    def property_r(self,val):
        self._property_r = val
#// end of properties

#Molecule freq method
    @property
    def freq(self):
        if self._freq is not None:
            return self._freq
        if self.Cluster:
            return self.Cluster.freq
        return "0.0000000"

    def to_Property(self, r = None ):
        """Will remove the property of individual atoms and put them in center of 
        mass of molecule instead"""
        if not self.LoProp and self.Property:
            warnings.warning("Can't convert to Property, Property is set")
        if r is None:
            r = self.com
        self.Property = self.sum_property
        for at in self:
            at.p = Property()
        self.is_Property = True



    @freq.setter
    def freq(self, val):
        self._freq = val

#Molecule snap method
    @property
    def snap(self):
        if self._snap is not None:
            return self._snap
        if self.Cluster:
            return self.Cluster.snap
        return "0"

    @snap.setter
    def snap(self, val):
        self._snap = val


#Scale bond
    def scale_bond(self, at1, at2, scale = 1.0):
        """Bend bond from at1 to at2 by a factor scale = 1.0"""
        self.populate_bonds()
        assert at2 in at1.bonds.values()
        ats = at1.get_async_bond( first = at2 )[1:]
        dr = ( at2.r - at1.r ) * scale - ( at2.r - at1.r )
        for at in ats:
            at.x, at.y, at.z = at.r + dr

#Bend angle in Molecule
    def bend(self, at1, at2, at3, theta = 0.0):
        """Bend angle between at1 and at2 by theta"""
        self.populate_angles()
        
#So we only rotate atoms, not the one in middle of angle bond
        ats = at2.get_async_bond( first = at1 )[1:]

        trans, r3, r2, r1 = center_and_xz( at2.r, at3.r, at1.r )

        for at in ats:
            at.x, at.y, at.z = at.r + trans
            at.x, at.y, at.z = np.einsum('ab,bc,cd,d', Rz_inv(r3), Ry(r2), Rz_inv(r1), at.r )
            at.x, at.y, at.z = np.einsum('ab,b', Ry_inv(theta), at.r )
            at.x, at.y, at.z = np.einsum('ab,bc,cd,d', Rz(r1), Ry_inv(r2), Rz(r3), at.r )
            at.x, at.y, at.z = at.r - trans

#Method of Molecule
    def potline(self, max_l = 2 , pol = 22, hyper=1, fmt = "%.5f ",
            prop_point = None,
            ):
        string = ""
        if self.is_Property:
#can override molecules property location if we want by prop_point keyword
            center = self.com
            if prop_point is not None:
                center = prop_point
            tmp_atom = Atom()
            tmp_atom.Molecule = self
            tmp_atom.x, tmp_atom.y, tmp_atom.z = center
            tmp_atom.Property = self.sum_property
            string += tmp_atom.potline( max_l, pol, hyper, fmt = fmt )
        else:
            for at in self:
                string += at.potline( max_l, pol, hyper, fmt = fmt )
        return string

    @property
    def cluster_order(self):
        if self._cluster_order is not None:
            return self._cluster_order
        return 0

    @cluster_order.setter
    def cluster_order(self, val):
        self._cluster_order = val


    def connect_everything(self):
        for a in self:
            a.Molecule = self

    def reflect_mol_plane(self, key = lambda x:(x[0].r, x[1].r, x[2].r),
            plane = 'xz'):
        """docstring for reflect"""
        p1, p2, p3 = key( self )
        origin = p1.copy()
        t, r1, r2, r3 = center_and_xz( p1, p2, p3 )
        self.t( -origin )

        R1 = Rz( r1 )
        R2_inv = Ry_inv( r2 )
        R3 = Rz( r3 )
        S = reflect( plane )
        R1_inv = np.einsum( 'ij->ji', R3 )
        R2 = np.einsum( 'ij->ji', R2_inv )
        R3_inv = np.einsum( 'ij->ji', R1 )

# Rotate each atom to fit key orientation, followed by plane reflection and rotation back
        for at in self:
            at.x, at.y, at.z = np.einsum( 'ab, bc, cd, de, ef, fg, gi, i', R3, R2_inv, R1, S, R3_inv, R2, R1_inv, at.r )
            at.p = at.p.transform_by_matrix( R1_inv )
            at.p = at.p.transform_by_matrix( R2 )
            at.p = at.p.transform_by_matrix( R3_inv )
            at.p = at.p.transform_by_matrix( S )
            at.p.rotate( r1, r2, r3)
        self.t( origin )

    def reflect(self, plane = 'xz' ):
        S = reflect( plane )
        for at in self:
            at.x, at.y, at.z = np.einsum( 'ij,j', S,  at.r)
            at.p = at.p.transform_by_matrix( S )

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
            try:
                return super(Molecule,self).__getitem__( item )
#Return False if looking up molecule that is empty
            except IndexError:
                warnings.warn("Looking up molecule that is empty")
                return None
#Will return different type if reloaded import in ipython
            #except TypeError:
            #    return self[ item ]
                
    def rotate_around(self, p1, p2, theta = 0.0):
        """Rotate All aomts clockwise around line formed from point p1 to point p2 by theta"""
        for at in self :
            at.x, at.y, at.z = rotate_point_by_two_points( at.r, p1, p2, theta)

    def inv_rotate(self, t1, t2, t3):
        """rotate all atoms and properties as
        1) inverse Z rotation by t1
        2) positive Y rotation by t2
        3) inverse Z rotation by t3
        """
#Put back in original point
        r1 = Rz_inv(t1)
        r2 = Ry(t2)
        r3 = Rz_inv(t3)
        for at in self:
            at.x, at.y, at.z = np.einsum('ab,bc,cd,d', r3, r2, r1, at.r )
            at.p = at.p.inv_rotate( t1, t2, t3 )

    def rotate(self, t1, t2, t3, inplace = False):
        """Rotate atomss and properties by t1, t2, t3
        t1 positive rotation around Z-axis
        t2 negative rotation around Y-axis
        t3 positive rotation around Z-axis
        """
        r1 = Rz(t1)
        r2 = Ry_inv(t2)
        r3 = Rz(t3)

        com = self.com.copy()
        if inplace:
            self.t( -com )
        if self.LoProp:
            for at in self:
                at.x, at.y, at.z = np.einsum('ab,bc,cd,d', r3, r2, r1, at.r )
                at.p = at.p.rotate( t1, t2, t3 )
        elif self.is_Property:
            for at in self:
                at.x, at.y, at.z = np.einsum('ab,bc,cd,d', r3, r2, r1, at.r )
            self.Property = self.Property.rotate( t1, t2, t3 )
        else:
            for at in self:
                at.x, at.y, at.z = np.einsum('ab,bc,cd,d', r3, r2, r1, at.r )
        if inplace:
            self.t( com )

    def template(self, max_l = 0, pol = 1, hyp = 0,
            label_func = lambda x: x.pdb_name,
            centered = None,
            decimal = 7):
        """Write out the template with properties for molecule
        if centered, will put the property on position given as array"""
        if len(label_func.func_code.co_varnames) != 1:
            print "Unsupported multi key function"
            raise SystemExit
        st_label = "_".join( label_func.func_code.co_names )
        st = "{\n"
        st += "'meta' : { 'label' : '%s', },\n" % st_label

        if centered:
            p = self.p
            st += ("( 'X', {0:8s}) : " + "{1:%df},\n"%decimal).format( "'charge'", p.q )
            tmp = "( 'X', {0:8s}) : [%s],\n"%(reduce(lambda a,x:a+x,map(lambda x: " {%d:.%df}, " %(x,decimal), range(1,4) )))

            st += tmp.format( "'dipole'", *(p.d.tolist()) )
            tmp = "( 'X', {0:8s}) : [%s],\n"%(reduce(lambda a,x:a+x,map(lambda x: " {%d:.%df}, " %(x,decimal), range(1,7) )))

            st += tmp.format(  "'quadrupole'", *p.Q.tolist() )
            tmp = "( 'X', {0:8s}) : [%s],\n"%(reduce(lambda a,x:a+x,map(lambda x: " {%d:.%df}, " %(x,decimal), range(1,7) )))

            st += tmp.format( "'alpha'", *p.a )
            tmp = "( 'X', {0:8s}) : [%s],\n"%(reduce(lambda a,x:a+x,map(lambda x: " {%d:.%df}, " %(x,decimal), range(1,11) )))

            st += tmp.format( "'beta'", *p.b )
            st += '}'
            return st

        for at in self:
            st += ("( {0:5s}, {1:8s}) : "+"{2:.%df},\n" %decimal).format( "'" +label_func(at)  +"'" , "'charge'", at.p.q )
            tmp = "( {0:5s}, {1:8s}) : [%s],\n"%(reduce(lambda a,x:a+x,map(lambda x: " {%d:.%df}, " %(x,decimal), range(2,5) )))

            st += tmp.format( "'" + label_func(at) + "'", "'dipole'", *(at.p.d.tolist()) )
            tmp = "( {0:5s}, {1:8s}) : [%s],\n"%(reduce(lambda a,x:a+x,map(lambda x: " {%d:.%df}, " %(x,decimal), range(2,8) )))

            st += tmp.format(  "'" + label_func(at) + "'", "'quadrupole'",*at.p.Q.tolist() )
            tmp = "( {0:5s}, {1:8s}) : [%s],\n"%(reduce(lambda a,x:a+x,map(lambda x: " {%d:.%df}, " %(x,decimal), range(2,8) )))

            st += tmp.format( "'" + label_func(at) + "'", "'alpha'",*at.p.a )
            tmp = "( {0:5s}, {1:8s}) : [%s],\n"%(reduce(lambda a,x:a+x,map(lambda x: " {%d:.%df}, " %(x,decimal), range(2,12) )))

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
            return 0

    @property
    def res_name(self):
        if self._res_name is not None:
            return self._res_name
        return "MOL"

    @res_name.setter
    def res_name(self, val):
        self._res_name = val
    
    
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
            item.Molecule = self
        else:
            logging.warning('Passed %s instance to Molecule' %type(item) )

    def initiate(self, item):
        item = item.copy()
        if isinstance( item, Atom ):
#Small hack to set molecule to LoProp if initiated with atom that has Property on it
            if item.p.q != 0.0 or ( item.p.d != np.zeros(3)).all():
                self.LoProp = True
            self.append( item )
            item.Molecule = self
        else:
            logging.warning('Passed %s instance to Molecule' %type(item) )



    def plot_2d(self, copy = True, smart = False, key = None, ):
        """Plots the 2D projetion of the molecule onto the y plane,
        PLAN: TODO:
        Later implement internal plane projection on arbitray two-vector plane
        """
        fig, ax = plt.subplots()
        #norm  = self.get_internal_plane()
        if copy:
            copy = self.copy()
        else:
            copy = self
        copy.populate_bonds()

        if smart:
            p1, p2, p3 = largest_triangle( [at.r for at in copy] )
            t, r1, r2, r3 = center_and_xz( p1, p2, p3 )
            v = np.cross( p3 - p1, p2 - p1 )
            norm = v / np.linalg.norm( v )
        else:
            norm = np.array( [0, 0, 1] )

        copy.t( t )
        copy.inv_rotate( r1, r2, r3 )


        if key is None:
            key = lambda x: x.element

        for at in copy:
            x, y = (at.r - np.einsum('i,i', norm, at.r)*norm)[:-1]
            ax.scatter( x, y, color = color_dict[at.element], linewidth=4 )
            ax.annotate(  key( at ) , (x, y ), (x, y+0.1) )
        
        x = np.linspace( -5, 5, 100 )
        y = np.linspace( -5, 5, 100 )

        #plot the bonds
#Plot bonds
        for each in copy:
            for key in each.bonds.values():
                ax.plot( [key.x, each.x],
                         [key.y, each.y],
                         color = 'black', linewidth = 0.25 )

        for at in copy:
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

    def save(self, fname = "molecule.p"):
        pickle.dump( self, open( fname, 'wb' ), protocol = 2 )

    @staticmethod
    def load(fname = 'molecule.p'):
        if not os.path.isfile( fname):
            raise IOError
        return pickle.load( open(fname, 'rb' ) )
    
    @property
    def b_proj(self):
        return b_para( self.b, self.p )

    def attach_properties(self, 
            model = "TIP3P_PDB",
            method = "B3LYP",
            basis = "ANO631",
            loprop = True,
            freq = "0.0",
            euler_key = lambda x: (x[0].r, x[1].r, x[2].r),
            template_key = lambda x: x.pdb_name,
            force_template = False,
            centered = None,
            ):
        """Attach property for Molecule method, by default TIP3P/HF/ANOPVDZ, static
        
        If the attachment is not LoProp, the default placement of the property
        is at the center-of-mass

        the centered argument can be used to put it elsewhere, i.e. center-of-charge
        """
        self.Property = None
        if centered is None:
            centered = self.com
        if isinstance(self, Water) and not force_template:
            template_key = lambda x: x.element + str(x.order)

        if isinstance(self, Water):
            R = np.einsum('ij->ji', get_rotation( self.o.r, (self.h1.r - self.h2.r)/2 + self.h2.r, self.h1.r) )
        else:
            R = np.einsum('ij->ji', get_rotation( *euler_key( self ) ))

        templ = Template().get( *(model, method, basis, loprop, freq) )
        if loprop:
            for at in self:
                at.p = Property.from_template( template_key(at), templ )
        else:
            self.Property = Property.from_template( 'X', templ )

        #t1, t2, t3 = self.get_euler( key = euler_key )
        for at in self:
            at.p.d = np.einsum('ij,j', R, at.p.d )
            at.p.a = s2ut( np.einsum('ai,bj,ij', R, R, ut2s(at.p.a) ) )
            at.p.Q = s2ut( np.einsum('ai,bj,ij', R, R, ut2s(at.p.Q) ) )
            at.p.b = s2ut( np.einsum('ai,bj,ck,ijk', R, R, R, ut2s(at.p.b) ) )

        if loprop:
            self.is_Property = False
            self.LoProp = True
        else:
            self.is_Property = True
            self.LoProp = False
            self.property_r = centered

            self.Property.d = np.einsum('ij,j', R, self.Property.d )
            self.Property.a = s2ut( np.einsum('ai,bj,ij', R, R, ut2s(self.Property.a) ) )
            self.Property.Q = s2ut( np.einsum('ai,bj,ij', R, R, ut2s(self.Property.Q) ) )
            self.Property.b = s2ut( np.einsum('ai,bj,ck,ijk', R, R, R, ut2s(self.Property.b) ) )

        self._freq = freq




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
        if isinstance( self, Water ):
            key = lambda x: (x.o.r, x.h2.r + (x.h1.r-x.h2.r)/2, x.h1.r)
        if rot_type == 'water':
            key = lambda x: (x.o.r, x.h2.r + (x.h1.r-x.h2.r)/2, x.h1.r)
        try:
            p1, p2, p3 = key( self )
        except IndexError:
            logging.error('Tried to get euler on Molecule with too few atoms')
        t, r1, r2, r3 = center_and_xz( p1, p2, p3 )
        return r1, r2, r3

    def props_from_targz(self,
            f_ = None,
            tmpdir = None,
            maxl = 2,
            pol = 22,
            hyp = 2,
            decimal = 7,
            freq = None,
            bonds = None,
            ):
        if tmpdir is None:
            tmpdir = "/tmp"
        if f_ is None:
            print "Supply .tar.gz file from dalton quadratic .QLOP calculation"
            return
        if freq == None:
            freq = 0.0
        else:
            freq = float(freq)
 
#Using Olavs external scripts
        tarfile.open( f_, 'r:gz' ).extractall( tmpdir )
        try:
            outpot = MolFrag( tmpdir = tmpdir,
                    max_l = maxl,
                    pol = pol,
                    pf = penalty_function( 2.0 ),
                    freqs = (freq,)
                    ).output_potential_file(
                            maxl = maxl,
                            pol = pol,
                            hyper = hyp,
                            decimal = decimal,
                            bond_centers = bonds,
                            )
        except IOError:
            print tmpdir

        f_at = lambda x: map(float, x.get_mol_line().split()[1:] )
        f_prop = lambda x: map(float,x.split()[1:4])

        self.populate_bonds()
        if bonds:
            relevant = sorted( self.get_ats_and_bonds(), key = f_at)
        else:
            relevant = sorted( self, key = f_at)
        for center in relevant:
            center.p = Property()

        lines = [ " ".join(l.split()) for l in outpot.split('\n') if len(l.split()) > 4 ]
        try:
            assert len( relevant ) == len( lines )
        except:
            logging.error("Some error went undetected despite creating .tar.gz and .out files")
            return


        for center, prop in zip( relevant, sorted( lines, key = f_prop )):
            prop = Property.from_propline( prop ,
                    maxl = maxl,
                    pol = pol,
                    hyper = hyp )
            center.p = prop
            center.p.freq = freq
            if isinstance( center, Bond ):
                my_label = center._Atom1.label
                other_bond = [b for b in center._Atom2.bonds if b._Atom2.label == my_label ][0]
                other_bond.p = prop

        self.LoProp = True

    def props_from_qm(self,
            tmpdir = None,
            dalpath = None,
            procs = 4,
            decimal = 7,
            maxl = 2,
            pol = 22,
            hyper = 2,
            method = 'hflin',
            basis = ['ano-1 2 0 0 0', 'ano-1 4 3 1 0', 'ano-2 5 4 1 0' ],
            dalexe = None,
            basdir = '/home/x_ignha/repos/dalton/basis',
            log = None,
            keep_outfile = False,
            keep_targzfile = False,
            freq = None,
            bonds = False,
            ):
        """
        Will generate a .mol file of itself, run a DALTON calculation as a
        childed process, get the properties back and put them on all atoms.

        Might take long time for large residues.
        """
        if freq == None:
            freq = 0.0
        else:
            freq = float(freq)
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

        dal = method + '.dal'
        mol = self.molfile_label() + '.mol'
        dal_full, mol_full = map( lambda x: os.path.join( tmpdir, x ), [dal,mol])
        if method == 'hfqua':
            open( dal, 'w').write( Generator.get_hfqua_dal( ) )
        elif method == 'b3lypqua':
            open( dal, 'w').write( Generator.get_b3lypqua_dal( ) )
        elif method == 'camb3lypqua':
            open( dal, 'w').write( Generator.get_camb3lypqua_dal( ) )
        elif method == 'ccsdlin':
            open( dal, 'w').write( Generator.get_ccsdlin_dal( ) )
        elif method == 'hflin':
            if freq:
                open( dal, 'w').write( Generator.get_hflin_freq_dal( freq = freq ) )
            else:
                open( dal, 'w').write( Generator.get_hflin_dal( ) )
        elif method == 'b3lyplin':
            if freq:
                open( dal, 'w').write( Generator.get_b3lyplin_freq_dal( freq = freq ) )
            else:
                open( dal, 'w').write( Generator.get_b3lyplin_dal() )
        elif method == 'camb3lyplin':
            if freq:
                open( dal, 'w').write( Generator.get_camb3lyplin_freq_dal( freq = freq ) )
            else:
                open( dal, 'w').write( Generator.get_camb3lyplin_dal() )
        else:
            print "wrong calculation type specified"
            return
        open( mol, 'w').write( self.get_mol( basis = basis) )

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
            print """The DALTON runscript can not be found. Please set the env variable DALTON to the dalton runscript or supply the script to props_from_qm directly as dalpath = <path-to-dalscript> """
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
                '-N', str(procs), '-D', '-t', tmpdir,
                dal, mol
                ], stdout = subprocess.PIPE,
                )
            p.stdout
            lines_iterator = iter(p.stdout.readline, b"")
            for line in lines_iterator:
                logging.debug( line ) 
            out, err = p.communicate()


            wrkdir = os.getcwd()
            out_wrk = os.path.join( wrkdir, "%s_%s.out" % tuple(map(lambda x:x.split('.')[0], [dal,mol])))
            tar_wrk = os.path.join( wrkdir, "%s_%s.tar.gz" % tuple(map(lambda x:x.split('.')[0], [dal,mol])))
            mol_wrk, dal_wrk = map( lambda x: os.path.join( wrkdir, x), [mol, dal] )
            mol_tmp, dal_tmp = map( lambda x: os.path.join( tmpdir, x ), [mol, dal])

            tar_tmp = "%s_%s.tar.gz" % tuple(map(lambda x:x.split('.')[0], [dal,mol]))
            of_tmp = "DALTON.OUT"


#Need to do this since late dalton scripts appends the tmp with separate PID
            real_tmp = find_dir( of_tmp, tmpdir )

#If smth happend to the dalton subprocess, the of will not exist and throw exception
            try:
                of_tmp, tar_tmp = map( lambda x: os.path.join( real_tmp, x ), [of_tmp, tar_tmp ] )
            except AttributeError:
                logging.error( "Some internal HPC specific error occured" )
                logging.error( "Will dump the .mol file of botched calculation" )
                logging.error( self.get_mol( basis = basis) )
                return

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
                            bond_centers = bonds, 
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


        try:
            assert len( self ) == len( lines )
        except:
            logging.error("Some error went undetected despite creating .tar.gz and .out files")
            return

        for at, prop in zip(sorted(self, key = f_at), sorted( lines, key = f_prop )):
            at.Property = Property.from_propline( prop ,
                    maxl = maxl,
                    pol = pol,
                    hyper = hyper )


#So that we do not pollute current directory with dalton outputs
#Also remove confliction inter-dalton calculation files
# For triolith specific calculations, remove all files in tmp
        for f_ in [os.path.join(real_tmp, f) for f in os.listdir(real_tmp) if os.path.isfile( os.path.join( real_tmp, f)) ]:
            os.remove( f_ )

        try:
            for f_ in [f for f in os.listdir(real_tmp) if os.path.isfile(f) ]:
                print os.path.join(real_tmp, f_)
                raise SystemExit
                os.remove( os.path.join(real_tmp, f_))
            for f_ in [mol, dal]:
                os.remove( f_ )
        except OSError:
            logging.error('Something wrint in trying to remove files in real_tmp')
        if not keep_outfile:
            os.remove( out_wrk )
        if not keep_targzfile:
            os.remove( tar_wrk )

        self.is_Property = False
        self.LoProp = True

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

    def populate_bonds(self, cluster = False):
        for at in self:
            at.bonds = []
#Implement later that it can only be called once
        if self.AA:
            conv = 1.0
        else:
            conv = 1/a0

#Populate bonds on cluster level

        if cluster:
            ats = self.Cluster.atoms
            for i1, at1 in enumerate( self ):
                for i2, at2 in enumerate( ats ):
                    if at1 == at2:
                        continue
                    if self[i1].dist_to_atom( ats[i2] ) < conv*bonding_cutoff[(self[i1].element, ats[i2].element)]:
                        b1 = Bond( self[i1], ats[i2] )
                        b2 = Bond( ats[i2], self[i1] )
                        self[i1].add_bond( b1 )
                        ats[i2].add_bond( b2 )
        else:
            for i2 in range( 1, len(self) ):
                for i1 in range( i2 ):
                    if self[i1].dist_to_atom( self[i2] ) < conv*bonding_cutoff[(self[i1].element, self[i2].element)]:

                        b1 = Bond( self[i1], self[i2] )
                        b2 = Bond( self[i2], self[i1] )
                        self[i1].add_bond( b1 )
                        self[i2].add_bond( b2 )

    def populate_angles(self):
        pass
    #    self.populate_bonds()
    #    for at1 in self:
    #        for at2 in [at2 for at2 in at1.bonds.values()]:
    #            for at3 in [at3 for at3 in at2.bonds.values() if at3 is not at1 ]:
    #                at1.angles[(at2.name,at3.name)] = at3

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

    def populate_dihedrals(self):

        self.populate_bonds()
        dihed = []
        for at1 in self:
            if at1.bonds == []:
                continue
            for at2 in at1.bonds.values():
                if at2.bonds == []:
                    continue
                for at3 in [a for a in at2.bonds.values() if a != at1]:
                    if at3.bonds == []:
                        continue
                    for at4 in [a for a in at3.bonds.values() if a != at2]:
                        at1.dihedrals[ (at2.name, at3.name, at4.name) ] = (at1,at2,at3,at4)

#For molecular properties on whole molecule, instead of atomic ones
    @property
    def Property(self):
        return self._Property
    @Property.setter
    def Property(self, val):
        if val is None or val is False:
            self._Property = val
        #if self.LoProp is None:
        #    warnings.warn("Setting property despite LoProp being False")
        else:
            self.LoProp = False
            self._Property = val
            self.is_Property = True

#Wrapper func for Molecule
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
        if self.Property:
            return self.Property
        coc = self.coc
        conv = 1.0
        p = Property()

        if self.AA:
            conv = 1/a0
        el_dip = np.array([ conv*(center.r-coc)*center.p.q for center in self.get_ats_and_bonds() ])

        nuc_dip = np.array([ conv*(center.r-coc)*charge_dict[center.element] for center in self.get_ats_and_bonds()])

        dip_lop = np.array([center.p.d for center in self.get_ats_and_bonds()])
        dip = el_dip + nuc_dip
        d = (dip + dip_lop).sum(axis=0)
        p.d = d
        for center in self.get_ats_and_bonds():
            p.q += center.p.q
            p.a += center.p.a
            p.b += center.p.b
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
        if type(args[0]) is int or type(args[0]) is float:
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


    def plot(self, copy = True, center = False, d = False, names = False, attr = False, key = False ):
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
        for at in copy:
            for bond in at.bonds:
                ax.plot( [at.x, bond._Atom2.x],
                         [at.y, bond._Atom2.y],
                         [at.z, bond._Atom2.z], color = 'black' )

        ax.plot( [0, 1, 0, 0, 0, 0], [0,0 ,0,1,0,0], [0,0,0,0,0,1] )
        ax.text( 1.1, 0, 0, "X", color = 'red' )
        ax.text( 0, 1.1, 0, "Y", color = 'red' )
        ax.text( 0, 0, 1.1, "Z", color = 'red' )
        if d:
            x = copy.coc[0]
            y = copy.coc[1]
            z = copy.coc[2]
            p = copy.p.d
            ax.plot( [x,x+p[0]], [y,y+p[1]], [z,z+p[2]], 'k-', linewidth = 3 )
            #ax.plot( [p[0]],[p[1]],[p[2]],'ko', markersize = 5, linewidth = 5 )
        for i in copy:
            ax.plot( [i.x], [i.y], [i.z], copy.style[i.element], linewidth= copy.linewidth[i.element] )
        if names:
            for i in copy:
                ax.text( i.x, i.y, i.z, i.name)
        if attr:
            for i in copy:
                ax.text( i.x, i.y, i.z, getattr(i, attr))

        ax.set_zlim3d( -5,5)
        plt.xlim(-5,5)
        plt.ylim(-5,5)
        plt.show()

    def get_mol(self, basis = ("ano-1 2 0 0 0", "ano-1 4 3 1 0",
        "ano-2 5 4 1 0" ), terminal = False):
        return self.get_mol_string( basis = basis, 
                terminal = terminal)

    def get_mol_string(self, basis = ("ano-1 2 0 0 0", "ano-1 4 3 1 0",
        "ano-2 5 4 1 0" ), terminal = False  ):
        if len( basis ) > 1:
            el_to_rowind = {"H" : 0, "He" : 0, "Li" : 1, "C" : 1, "O" : 1, "N" : 1,
                    "S" : 2, "P" : 2 }
        else:
            el_to_rowind = {"H" : 0, "C" : 0, "O" : 0, "N" : 0, "S" : 0, "P" :0 }
        st = ""
        s_ = ""
        if self.AA: s_ += " Angstrom"

        charge = 0.0
# Start Residue specifics that take terminal charge into account
        if terminal:
            if not self.is_concap:
                if self.n_term:
                    charge += 1
                elif self.c_term:
                    charge -= 1
                if self.res_name in res_dict:
                    if res_dict[ self.res_name ] in chargeDict:
                        charge += chargeDict[ res_dict[ self.res_name] ]
            if self._level == 3:
                charge = 0
                for at in self:
                    if at.pdb_name == "H1":
                        charge += 1
                        break
                    if at.pdb_name == "OC1":
                        charge -= 1
                        break
##  // end of Residue
        ats = sorted( self, key = lambda x: (x.element,) + (x.x, x.y, x.z) ) 
        uni = sorted(unique([ at.element for at in ats ]), key = lambda x: charge_dict[x] )

        st += "ATOMBASIS\n\n\nAtomtypes=%d Charge=%d Nosymm%s\n" %(len(uni), charge,  s_)
        for el in uni:
            st += "Charge=%s Atoms=%d Basis=%s\n" %( str(charge_dict[el]),
                    len( [all_el for all_el in ats if (all_el.element == el)] ),
                    basis[ el_to_rowind[el] ])
            for i in [all_el for all_el in ats if (all_el.element == el) ]:
                st += i.get_mol_line()
        return st

    def get_bond_and_xyz(self, ):
        st = ''
        ats = sorted( self, key = lambda x: (x.element,) + (x.x, x.y, x.z) ) 
        uni = sorted(unique([ at.element for at in ats ]), key = lambda x: charge_dict[x] )

        self.populate_bonds()
        bonds_outputted = []
        for el in uni:
            for at in [all_el for all_el in ats if (all_el.element == el) ]:
                st += at.get_mol_line()
                for b in at.bonds:
                    if any( np.allclose(b.r == x, atol = 1e-7) for x in bonds_outputted ):
                        continue
                    st += b.get_mol_line( lab = 'XX-%s-%s' %(at.pdb_name, b._Atom2.pdb_name))
                    bonds_outputted.append( b.r )
                         

        return st



    def get_pdb_string(self):
        st = """"""
        for at in self:
            st += at.pdb_string()
        return st

#Molecule copy method
    def copy_self(self):
        return self.copy()

    def copy(self, shallow = False):
        if shallow:
            return copymod.copy(self)
        else:
            return copymod.deepcopy(self)

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
    def from_mol_file( molfile, in_AA = True, out_AA = True):
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
                    pdb_name = lab.split('-')[3]
                else:
                    element = lab.split('-')[2][0]
                    pdb_name = lab.split('-')[2]
                kwargs = { "AA": in_AA,
                        "element" :  element,
                        "x" : matched[1],
                        "y" : matched[2],
                        "z" : matched[3],
                        "pdb_name" : pdb_name
                        }
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

    @property
    def AA(self):
        AA = self[0].AA
        for at in self:
            try:
                assert AA == at.AA
            except AssertionError:
                logging.error("Not all atoms same unit in molecule %s" %self)
        return AA

    @AA.setter
    def AA(self, val):
        assert type(val) == bool
        self._AA = val


    def to_AU(self):
        for at in self:
            at.to_AU()

    def to_AA(self):
        for at in self:
            at.to_AA()

    def byname(self,label ):
        return self.get_atom_by_name( label.upper() )
    def get_atom_by_name(self, label ):
        for at in self:
            if at.name == label:
                return at
        return None

    def abl(self, label, dup = False):
        """Wrapper for get_atom_by_label"""
        return self.get_atom_by_label( label = label, dup = dup )

    def abp(self, label, dup = False):
        """Wrapper for get_atom_by_pdbname"""
        return self.get_atom_by_pdbname( label = label, dup = dup )

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

    def get_atom_by_label(self, label, dup = False):
        for i in self:
            if i.label == label:
                return i
        logging.warning( "No %s in %s" %(label, self) )
        return None



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

        self._coc = None

        self.in_qm = False
        self.in_mm = False
        self.in_qmmm = False

#Since we added the atoms already, we need to define .o, .h1, and .h2
#the element + index of occurance of atom in molecule must match here for proper
#Water initialization from list of atoms
        if len(self) > 0:
            for atom in self: 
                if atom.name == 'O1':
                    self.o = atom
                if atom.name == 'H2':
                    self.h1 = atom
                if atom.name == 'H3':
                    self.h2 = atom



#Water - Properties for grabbing uniquely PDB named atoms
    @property
    def OW(self):
        return self.get_atom_by_pdbname( 'OW' )
    @property
    def HW1(self):
        return self.get_atom_by_pdbname( 'HW1' )
    @property
    def HW2(self):
        return self.get_atom_by_pdbname( 'HW2' )

#If initializing water from list of 3 atoms, set H, O for given atoms

    def copy(self, shallow = False):
        return super(Water, self).copy( shallow = shallow )

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

other models:
    sign_change : 
    r = 0.958019
    theta = 104.50


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
        elif model == 'sign_change':
            r_oh = 0.958019
            a_hoh =  104.50
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
        w = Water( )
        w.add_atom( o )
        w.add_atom( h1 )
        w.add_atom( h2 )
        w.o.pdb_name = 'OW'
        w.h1.pdb_name = 'HW1'
        w.h2.pdb_name = 'HW2'
        if worst:
            w.populate_bonds()
            w.populate_angles()
            w.h1.scale_angle( 0.988 )
            w.h1.scale_bond( 0.985 )
            w.h2.scale_bond( 1.015 )
            w.inv_rotate( *w.get_euler() )
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
    def add_atom(self, atom):
        """
Override list append method, will add up to 3 atoms,
1 must be oxygen, 2 must be hydrogens.

.. code:: python

    >>> m = Water()
    >>> m.add_atom( Atom( z = 0.11, element = 'H' ) )
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
            d1 = hyd1.dist_to_point( np.array([1,1,1]) )
            d2 = hyd2.dist_to_point( np.array([1,1,1]) )
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


class Cluster(list):
    """
**Molecule container which groups molecules into quantum mechanics, molecular mechanics, and qm/mm regions for easy in generating input files for QM/MM.**
"""
    def __init__(self, *args, **kwargs):
        """ Can be a list of molecules or a list of clusters"""

        super(Cluster, self).__init__()

        copy = False
        if kwargs != {}:
            copy = kwargs.get( "copy", False )

        self._chain_id = None
        self._freq = None
        self._snap = None
        self._system_order = None
        self._System = None
        self.Property = None
        self.atom_list = []

        if type(args) == tuple:
            if len(args) == 1:
                if type(args[0]) == list:
                    for i in args[0]:
                        self.add( i, copy = copy )
                else:
                    self.add( args[0], copy = copy )
            else:
                for item in args:
                    self.add( item, copy = copy )

    def get_distance_matrix(self):

        N = len(self.atoms)
        mat = np.zeros( (N,N ) )
        for i in range(1, self.atoms):
            for j in range( i ):
                mat[i, j] = self.atoms[i].dist_to_atom( self.atoms[j] )
        return mat


    @staticmethod
    def rand_water_cluster(N = 10):
        c = Cluster()
        cnt = 0 
        while cnt < N:
            w = Water.get_standard()
            w.rotate( *np.random.uniform( 0, 2*np.pi, 3) )
            w.t( np.random.uniform( -25, 25, 3) )
            if c.mol_too_close( w ):
                continue
            c.add( w )
            cnt += 1
        return c

#New implementation for Angstroms in clusters
    @property
    def AA(self):
        AA = [at for res in self for at in res ][0].AA
        for each in [at for res in self for at in res ]:
            try:
                assert each.AA == AA
            except AssertionError:
                logging.error("All objects in cluster are not of same unit")
        return AA
    @AA.setter
    def AA(self, val):
        assert type(val) == bool
        self._AA = val

    def dump_xyz(self, s = 'tmp' ):
        """Quickl;y wrote xyz file for fancier visualization in avogadro/vmd"""
        if s.endswith( '.xyz' ):
            s = s[:-4]
        open( s + '.xyz', 'w' ).write( self.get_xyz_string() )

#Cluster freq method
    @property
    def system_order(self):
        if self._system_order:
            return self._system_order
        return 0
#Cluster 
    @system_order.setter
    def system_order(self, val):
        self._system_order = val

#Cluster
    @property
    def System(self):
        return self._System
#Cluster
    @System.setter
    def System(self, val):
        self._System = val

#Cluster freq method
    @property
    def freq(self):
        if self._freq is not None:
            return self._freq
        if self.System:
            return self.System.freq
        return "0.0000000"

    @freq.setter
    def freq(self, val):
        self._freq = val
#Cluster snap method
    @property
    def snap(self):
        if self._snap is not None:
            return self._snap
        if self.System:
            return self.System.snap
        return "0"

    @snap.setter
    def snap(self, val):
        self._snap = val

#Cluster method
    @property
    def chain_id(self):
        if self._chain_id:
            return self._chain_id
        return 'X'

    @chain_id.setter
    def chain_id(self, val ):
        self._chain_id = val


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


#Cluster level
    def transfer_props(self, at_list, 
            transfer = { 'charge' : 0,
                'quadrupole' : 0,
                'dipole' : 0,
                'alpha' : 0,
                'beta' :0},
            center = None ):
        """
        Convenient wrapper for clusters with several interconnected
        Molecules.

        splitter puts atoms in seperate lists belonging to same molecules

        Make sure that these are not interconnected

        Note: the populate_bonds method must be run before this 
        in order to get correct transfering behaviour
        """
        

        mol_arr = splitter( at_list, key = lambda x: x.Molecule )

        for ind, ats in enumerate( mol_arr ):
            mol_arr[ ind ][0].Molecule.transfer_props( ats,
                    transfer = transfer,
                    center = center )

    @property
    def density(self):
        """Return the density in SI units kg/m^3"""
        N_A = 6.02214129e+23

#g/mol
        M = sum( [at.mass for mol in self for at in mol] )
#AU**3 or AA**3
        v = convex_hull_volume( np.array([at.r for mol in self for at in mol]))
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

    def get_pdlist(self, 
            max_l = 1,
            pol = 22,
            hyp = 1,
            model = 'pointdipole',
            bonds = False,
            cell = False,
            cell_cutoff = 25.0, 
            ):
        """Given cutoff in Angstromgs, will return a GassuanQuadrupoleList
        where atomic within AA_cutoff between different interacting segments
        
        has a damped gaussian """
        from applequist.particles import PointDipoleList
        from applequist.gaussian import GaussianQuadrupoleList
        from applequist.thole import TholeList

        opts = { 'pointdipole' : PointDipoleList,
                'point' :PointDipoleList,
                'thole' : TholeList,
                'gaussian' :GaussianQuadrupoleList,
                'p' : PointDipoleList,
                't' : TholeList,
                'g' :GaussianQuadrupoleList,
                'add' :PointDipoleList,
                }

        aa = self.AA
        self.to_AU()
        if cell:
            init_f = functools.partial( opts[model].cell_from_string, co = float(cell_cutoff) ) 
        else:
            init_f = opts[model].from_string
        g = init_f( self.get_qmmm_pot_string(max_l = max_l, pol = pol, hyp = hyp, bonds = bonds))

#Set all particles to group 1, meaning no interaction will take place
        if model == 'add':
            for particle in g:
                particle.group = 1
        if aa:
            self.to_AA()
        return g

    def g_list_from_damped(self, 
            max_l = 1,
            pol = 22,
            hyp = 1,
            rq = 1e-9,
            rp = 1e-9,
            AA_cutoff = 1.5,
            nullify = False,
            model = 'thole',
            cell = False,
            cell_cutoff = 25.0, 
            ):
        """Given cutoff in Angstromgs, will return a GassuanQuadrupoleList
        where atomic within AA_cutoff between different interacting segments
        
        has a damped gaussian """
        from .pd.particles import PointDipoleList
        from .pd.gaussian import GaussianQuadrupoleList
        from .pd.thole import TholeList

        opts = { 'pointdipole' : PointDipoleList,
                'point' :PointDipoleList,
                'thole' : TholeList,
                'gaussian' :GaussianQuadrupoleList,
                'p' : PointDipoleList,
                't' : TholeList,
                'g' :GaussianQuadrupoleList,
                }

        aa = self.AA
        self.to_AU()
        if cell:
            init_f = functools.partial( opts[model].cell_from_string, co = float(cell_cutoff) ) 
        else:
            init_f = opts[model].from_string
        g = init_f( self.get_qmmm_pot_string(max_l = max_l, pol = pol, hyp = hyp))

        for atom, res in map( lambda x: [x, x.Molecule], self.min_dist_atoms_separate_res(AA_cutoff) ):
            ind = reduce( lambda a, x: a + len(x), res.Chain[ : res.Chain.index(res) ] , 0) + atom.order
            g[ ind ]._R_q = rq
            g[ ind ]._R_p = rp
            if nullify:
                g[ ind ]._q = 0.0
                g[ ind ]._p0 = np.zeros( (3,) )
                g[ ind ]._Q0 = np.zeros( (3,3,) )
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
        basis = ("ano-1 2 0 0 0", "ano-1 4 3 1 0", "ano-2 5 4 1 0" ) ):
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

        """Given pdb/mol/xyz  file return a Cluster with all separate molecules"""
        if fil.endswith('.xyz'):
            with open(fil,'r') as f_:
                pass
    def min_dist_atoms_separate_res(self, AA_cutoff = 1.5 ):
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
        return unique( tmp )



     
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

        return unique(min_ats)


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
        """Docstring for __eq__ in cluster method, will equate to """
        return id(self) == id(other)
        
        
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
# Atoms list for cluster
    @property
    def atoms(self):
        return [at for mol in self for at in mol]
   
# Center of charge of Cluster
    @property
    def coc(self):
    #obj should be atom
        return sum( [at.r * charge_dict[at.element] for mol in self for at in mol])\
                /sum( map(float,[charge_dict[at.element] for mol in self for at in mol]) )


    def plot(self, copy = True, center = False, d = False, names = False, attr = False, key = False ):
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


#Plot in nice xyz axis
        copy.populate_bonds()
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d' )

#Plot bonds
        for at in copy.atoms:
            for bond in at.bonds:
                ax.plot( [at.x, bond._Atom2.x],
                         [at.y, bond._Atom2.y],
                         [at.z, bond._Atom2.z], color = 'black' )

        ax.plot( [0, 1, 0, 0, 0, 0], [0,0 ,0,1,0,0], [0,0,0,0,0,1] )
        ax.text( 1.1, 0, 0, "X", color = 'red' )
        ax.text( 0, 1.1, 0, "Y", color = 'red' )
        ax.text( 0, 0, 1.1, "Z", color = 'red' )

        for i in copy:
            for j in i:
                ax.plot( [j.x], [j.y], [j.z], j.Molecule.style[j.element], linewidth= j.Molecule.linewidth[j.element] )
        if names:
            for mol in copy:
                for at in mol:
                    ax.text( at.x, at.y, at.z, at.name)
        if attr:
            for mol in copy:
                for at in mol:
                    ax.text( at.x, at.y, at.z, getattr(at, attr ) )

        ax.set_zlim3d( -5,5)
        plt.xlim(-5,5)
        plt.ylim(-5,5)
        plt.show()


    def get_mol_string(self, basis = ("ano-1 2 0 0 0", "ano-1 4 3 1 0", "ano-2 5 4 1 0" ) ):
# If basis len is more than one, treat it like molecular ano type
# Where first element is for first row of elements
        if len( basis ) > 1:
            # Set row index number to periodic table one
            el_to_rowind = {"H" : 0, "C" : 1, "O" : 1, "N" : 1 , "S" : 2  }
        else:
            # Keep all 0, since basis is only len 1
            el_to_rowind = {"H" : 0, "C" : 0, "O" : 0, "N" : 0, "S" : 0 }

        st = ""
        comm1 = "Comment 1"
        comm2 = "Comment 2"
        uni = unique([ at.element for mol in self for at in mol])
        s_ = ""
        if self.AA: 
            s_ += "Angstrom"

        st += "ATOMBASIS\n%s\n%s\nAtomtypes=%d Charge=0 Nosymm %s\n" %( \
                comm1,
                comm2,
                len(uni),
                s_)
        for el in uni:
            st += "Charge=%s Atoms=%d Basis=%s\n" %( str(charge_dict[el]),
                    len( [all_el for mol in self for all_el in mol if (all_el.element == el)] ),
                     basis[ el_to_rowind[el] ] )
            for i in [all_el for mol in self for all_el in mol if (all_el.element == el) ]:
                st += "{0:5s}{1:10.5f}{2:10.5f}{3:10.5f}\n".format( i.element, i.x, i.y, i.z )
        return st




    def get_qm_mol_string(self, basis = ("ano-1 2 0 0 0", "ano-1 4 3 1 0", "ano-2 5 4 1 0" ) ):
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
        uni = unique([ at.element for mol in self for at in mol if mol.in_qm])
        s_ = ""
        if self.AA: s_ += "Angstrom"

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
#set center to true to avoid emply loprop atom lines in output, good to make matrices smaller in pointdipole list
            center = False,
#Set ignore_qmmm to false to only write qmmm .pot file for molecues in mm region
            ignore_qmmm = True,
            bonds = False,
            fmt = "%.7f "
            ):

#make sure that each unique residue is in separate residue in POT output
        self.order_mm_mols()

        if self.AA:
            st = "AA\n"
        else:
            st = "AU\n"

# Old qmmm format requires integer at end to properly read charges ?
#Disable as it introduce a bug in reading zero beta inputs
        #if hyp == 0:
        #    hyp_int = 1

#If bonds is true, put all bond points in
        if bonds:
            st += "%d %d %d %d\n" % (sum([len(center.get_ats_and_bonds()) for center in self if center.LoProp]) + len([m for m in self if m.is_Property]), max_l, pol, 1 )

            st += "".join( [center.potline(max_l, pol, hyp, fmt = fmt) for mol in self for center in mol.get_ats_and_bonds() if mol.LoProp] )

            st += "".join( [mol.potline(max_l, pol, hyp, fmt = fmt) for mol in self if mol.is_Property ] )

        else:
            if ignore_qmmm:
                st += "%d %d %d %d\n" % (sum([len(center) for center in self if center.LoProp]) + len([m for m in self if m.is_Property]), max_l, pol, 1 )

                st += "".join( [center.potline(max_l, pol, hyp, fmt = fmt) for mol in self for center in mol if mol.LoProp] )

                st += "".join( [mol.potline(max_l, pol, hyp, fmt = fmt) for mol in self if mol.is_Property ] )
            else:
                st += "%d %d %d %d\n" % (sum([len(i) for i in self if i.in_mm ]), 
                        max_l, pol, 1 )
                st += "".join( [at.potline(max_l, pol, hyp, fmt = fmt) for mol in self for at in mol if mol.in_mm and mol.LoProp] )
                st += "".join( [mol.potline(max_l, pol, hyp, fmt = fmt) for mol in self if mol.in_mm and mol.is_Property ] )
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
        self.sort( key = lambda x:x.is_Property )
        for mol in self:
            mol.cluster_order = cnt
            cnt += 1

    def update_water_props(self, model = "TIP3P",
            method = "HF", basis = "ANOPVDZ", dist = True,
            freq = "0.0"):

        kwargs_dict = Template().get( *(model, method, basis,
            dist , freq ))
        for wat in self:
            t1, t2, t3 = wat.get_euler()
            for at in wat:
                Property.add_prop_from_template( at, kwargs_dict )
                at.Property.rotate( t1, t2, t3)

    @staticmethod
    def get_water_cluster_from_string( _str, in_AA,
            out_AA, N_waters = 1000, file_ending = '.mol',
            md_type = 'xinli', ):

        """
Return a cluster of water molecules given file.

.. code:: python

    >>> c = Cluster.get_water_cluster( 'somefile.mol' , in_AA = False, out_AA = False, N_waters = 10 )
    >>> print len( c )
    10

"""
        atoms = []
        c = Cluster()
        lines = _str.split( '\n' )

        if file_ending == '.mol' or file_ending == '.xyz':
            pat_xyz = re.compile(r'^\s*(\w+)\s+(-*\d*.+\d+)\s+(-*\d*.+\d+)\s+(-*\d*.+\d+) *$')
            for i in lines:
                if pat_xyz.match(i):
                    f = pat_xyz.match(i).groups()
                    matched = pat_xyz.match(i).groups()
                    kwargs = { "element" :  matched[0],
                            "AA": in_AA,
                            "x" : matched[1],
                            "y" : matched[2], "z" : matched[3] }
                    tmpAtom = Atom( **kwargs )
                    atoms.append( tmpAtom )

        elif file_ending == '.pdb':
            pat1 = re.compile(r'^(ATOM|HETATM)')
#Temporary atom numbering so that it is compatible with PEQM reader in dalton
            for i in lines:
                if pat1.search(i):
                    n = i[11:16].strip()
                    allowed_res_names = ['MOL', 'SOL', 'WAT']
#Ignore residue names that are not in allowed ater residue names
                    if i[17:20].strip() not in allowed_res_names:
                        continue
                    if n in [ "SW", "DW", "MW" ]:
                        continue
                    kwargs = {
                            "AA" : in_AA,
                            "x" : float(i[30:38].strip()),
                            "y" : float(i[38:46].strip()),
                            "z" : float(i[46:54].strip()),
                            "pdb_name" : i[11:16].strip(),
                            "element": i[11:16].strip()[0],
                            "number" : i[6:11].strip()  }
                    tmpAtom = Atom( **kwargs )
                    atoms.append( tmpAtom )
        elif file_ending == '.out':
            pat_xyz = re.compile(r'^(\w+)\s+(-*\d*.+\d+)\s+(-*\d*.+\d+)\s+(-*\d*.+\d+) *$')
            for i in lines:
                if pat_xyz.match(i):
                    f = pat_xyz.match(i).groups()
                    tmpAtom = Atom(f[0][0], float(f[1]), float(f[2]), float(f[3]), 0)
                    atoms.append( tmpAtom )
        elif file_ending == '.log':
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
#Find the box encompassing all atoms
        xmin, ymin, zmin = np.inf, np.inf, np.inf
        xmax, ymax, zmax = -np.inf, -np.inf, -np.inf
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
#Dead center of the box
        center = np.array(map( lambda x: (x[1] - x[0])/2.0+x[0], zip([xmin, ymin, zmin], [xmax, ymax, zmax] )) )

        from .dstruct import Cell
        cell = Cell( [xmin, ymin, zmin], [xmax, ymax, zmax], 3 )
        for at in atoms:
            cell.add_atom(at )

        
        wlist = []
        for i in cell:
            if i.element == "H":
                continue
            if i.in_water:
                continue
            tmp = Water( AA = in_AA )
#__Water__.append() method will update the waters residue number and center coordinate
#When all atoms are there
#Right now NOT center-of-mass
            i.in_water= True
            tmp.add_atom(i)
            for j in cell.get_closest(i):
                if j.element == "O":
                    continue
                if j.in_water:
                    continue
#1.05 because sometimes spc water lengths can be over 1.01
                if in_AA:
                    if i.dist_to_atom(j) <= 1.05:
                        j.in_water = True
                        tmp.add_atom( j )
                else:
                    if i.dist_to_atom(j) <= 1.05/a0:
                        j.in_water = True
                        tmp.add_atom( j )
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
#If in cartesian:
        for wat in c:
            for atom in wat:
                atom._res_id = wat.res_id

        if in_AA or file_ending == '.log':
            if not out_AA:
                for wat in c:
                    wat.to_AU()
        if not in_AA and out_AA:
            c.to_AA()
        for wat in c:
            wat.o.order = 1
            wat.h1.order = 2
            wat.h2.order = 3
        c.set_qm_mm(100)
        if md_type == 'xinli':
            for wat in c:
                wat.o.pdb_name = 'OW'
                wat.h1.pdb_name = 'HW1'
                wat.h2.pdb_name = 'HW2'
        return c

    @staticmethod
    def get_water_cluster( fname , in_AA = False, out_AA = False , N_waters = 1000, md_type = 'xinli' ):
        """ md_type = 'xinli' assigns the standard pdb names to oxygen and hydrogen
        o.pdb_name = OW
        o.h1 = HW1
        o.h2 = HW2
        """
        file_ending = '.' + fname.split('.')[-1]
        return Cluster.get_water_cluster_from_string( open(fname).read(),
                in_AA = in_AA, out_AA = out_AA, N_waters = N_waters,
                file_ending = file_ending, md_type = md_type)

    def mol_too_close(self, mol, dist = 2.5):
        for mols in self:
            for ats in mols:
                for at in mol:
                    if at.dist_to_atom( ats ) < dist:
                        return True
        return False

    def attach_properties(self, 
            model = "TIP3P_PDB",
            method = "B3LYP",
            basis = "ANO631",
            loprop = True,
            freq = "0.0",
            euler_key = lambda x: (x[0].r, x[1].r, x[2].r),
            template_key = lambda x: x.pdb_name,
            force_template = False,
            centered = None,
            ):
        """Attach properties to all molecules in this cluster"""
        for mol in self:
            mol.attach_properties( model = model,
                    method = method,
                    basis = basis,
                    loprop = loprop,
                    freq = freq,
                    euler_key = euler_key,
                    template_key = template_key,
                    force_template = force_template,
                    centered = centered )

#Add method for cluster
    def add(self, item, copy = False ):
        if isinstance( item , Molecule ):
            self.add_mol( item, copy = copy )
        elif isinstance( item, Atom) :
            self.add_atom( item, copy = copy )
        else:
            logging.warning( 'Tried to pass other instance than Atom or Molecule to Cluster' )

    def add_mol(self, mol, copy = False ):
        if isinstance( mol , Molecule ):
            if copy:
                mol = copymod.deepcopy( mol )
            if mol.chain_id is not 'X':
                self._chain_id = mol.chain_id
            mol.Cluster = self
            super( Cluster, self ).append( mol )
        elif type( mol ) == list:
            for each in mol:
                if copy:
                    each = copymod.deepcopy( each )
                if each.chain_id is not 'X':
                    self._chain_id = each.chain_id
                super( Cluster, self ).append( each )
                each.Cluster = each

    def add_atom(self, at, copy = False ):
        for i, iat in enumerate(at):
#If we create a slice of a cluster, retain the chain ID in the new one,solves a lot of problems
            if copy:
                iat = copymod.deepcopy( iat )
            if iat.chain_id is not 'X':
                print iat
                self._chain_id = iat.chain_id
                print self.chain_id
            self.append( iat )
            iat.Cluster = self

#Special cluster method when dealing with different molecules in a clusters
# By defalt only connect atoms in seperate peptides, meaning carbons
    def populate_bonds(self, cluster = False):
        """ given argument cluster = True, will populate bonds"""
        for at in self.atoms:
            at.bonds = []
#Implement later that it can only be called once
        if self.AA:
            conv = 1.0
        else:
            conv = 1/a0

#Populate bonds on cluster level
        for i2 in range( 1, len(self.atoms) ):
            for i1 in range( i2 ):
                if self.atoms[i1].dist_to_atom( self.atoms[i2] ) < conv*bonding_cutoff[(self.atoms[i1].element, self.atoms[i2].element)]:
                    b = Bond( self.atoms[i1], self.atoms[i2] )
                    self.atoms[i1].add_bond( b )
                    self.atoms[i2].add_bond( b )

#Cluster method for angles
    def populate_angles(self):
        self.populate_bonds()
        for at1 in self.atoms:
            for at2 in [at2 for at2 in at1.bonds.values()]:
                for at3 in [at3 for at3 in at2.bonds.values() if at3 is not at1 ]:
                    at1.angles[(at2.name,at3.name)] = at3

#Bend angle in Cluster, should be fine as long as bonds are correctly connected
    def bend(self, at1, at2, at3, theta = 0.0):
        """Bend angle between at1 and at2 by theta"""
        self.populate_angles()
        
#So we only rotate atoms, not the one in middle of angle bond
        ats = at2.get_async_bond( first = at1 )[1:]
        trans, r3, r2, r1 = center_and_xz( at2.r, at3.r, at1.r )
        for at in ats:
            at.x, at.y, at.z = at.r + trans
            at.x, at.y, at.z = np.einsum('ab,bc,cd,d', Rz_inv(r3), Ry(r2), Rz_inv(r1), at.r )
            at.x, at.y, at.z = np.einsum('ab,b', Ry_inv(theta), at.r )
            at.x, at.y, at.z = np.einsum('ab,bc,cd,d', Rz(r1), Ry_inv(r2), Rz(r3), at.r )
            at.x, at.y, at.z = at.r - trans



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
        return copymod.deepcopy( self )

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
            decimal = 7,
            maxl = 2,
            pol = 22,
            hyper = 2,
            method = 'hflin',
            env = os.environ,
            basis = ['ano-1 2 0 0', 'ano-1 4 3 1 0', 'ano-2 5 4 1 0' ],
            dalexe = None,
            basdir = '/home/x_ignha/repos/dalton/basis',
            log = None,
            keep_outfile = False,
            freq = None,
            bonds = False,
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
                    bonds = bonds,
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

#Function for cluster
    @property
    def p(self):
        return self.sum_property

    @property
    def alpha(self):
        """
        Return the sum properties of all molecules in cluster
        Now it is dead wrong, need to adjust dipoles and quadrupoles to coc
        """
        p = Property()
        alphas =  [mol.p.a for mol in self ]
        for alpha in alphas:
            p.a += alpha
        return p.a

    @property
    def beta(self):
        """
        asdf
        """
        p = Property()
        betas =  [mol.p.b for mol in self ]
        for beta in betas:
            p.b += beta
        return p.b

    @property
    def sum_property(self):
        """
        Return the sum properties of all molecules in cluster
        Now it is dead wrong, need to adjust dipoles and quadrupoles to coc

        Important bugfix: if molecules have molecular centered properties,
        this sum_property will return zero as all atoms and bonds have zero

        """
        coc = self.coc
        conv = 1.0
        p = Property()

        if self.AA:
            conv = 1/a0

#el_dip and nuc_dip will add the property differently depending on if it 
#is a LoProp placement on the molecule, or if the whole molecule is represented
#by one property at a selected center
        el_dip = np.array( [ conv*(center.coc - coc)*center.p.q for center in self ] )
        dip_lop = [mol.p.d for mol in self ]

        d = ( el_dip + dip_lop).sum( axis = 0 )

        p.d = d

        c_mols =  [mol.p for mol in self ]

        for center in c_mols:
            p.q += center.q
            p.a += center.a
            p.b += center.b
        return p

    def to_AA(self):
        if not self.AA:
            for i in self:
                i.to_AA()
        return self
    def get_atom_by_label(self, val):
        for at in self.atoms:
            if at.label == val:
                return at
        return None

    def get_atom_by_clabel(self, val):
        for at in self.atoms:
            if at.clabel == val:
                return at
        return None
     

    def to_AU(self):
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
