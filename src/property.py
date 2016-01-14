__all__ = [ 'Property' ]

import numpy as np
import logging
import warnings

from .utilz import ut2s, s2ut, Ry, Rz, Ry_inv, Rz_inv


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

    def __getitem__(self, item):


        new = { 'c' : 'charge', 'd' : 'dipole', 'q' :'quadrupole',
                'a' : 'alpha', 'b' : 'beta' }
        try:
            key = new[ item[0].lower() ]
        except KeyError:
            logging.error("unknown command for getting item of Property")
            raise SystemExit
        return super(Property, self).__getitem__( key )


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

    def __div__(self, other):
        tmp = Property()
        for i, prop in enumerate(self):
            tmp[prop] = np.array( self[prop] )/float(other)
        return tmp

    def is_null(self):
        empty = True
        for key, val in self.iteritems():
            if not np.allclose( np.zeros( np.array((val,)).shape ), np.array((val,)) , atol = 1e-14):
                empty = False
        return empty



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
        return np.einsum('ijj,i', ut2s(self.b), self.d)/np.linalg.norm(self.d) #* #self.d / np.linalg.norm( self.d )

#Method of Property
    def potline(self, max_l =2 , pol= 22, hyper=1, fmt = "%.7f "):
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
            p['alpha'][3] = iso
            p['alpha'][5] = iso
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
                print "'( %s, %s )' not found in provided template" %(at_string, p)
                warnings.warn("Could not find in provided template" )
                return Property()

        p = Property() 
        for key in template:
            if key[0] == at_string:
                for each in all_props:
                    p[ each ] = template[ (at_string, each ) ]
        return p

    def inv_rotate( self, t1, t2, t3, plane = 'xz' ):
        """Rotate all properties by t1, t2, t3
        t1 negative rotation around Z-axis
        t2 positiv rotation around Y-axis
        t3 negative rotation around Z-axis
        """
        rots = { 'xz' : (Rz_inv(t1), Ry(t2), Rz_inv(t3) ),
                 'xy' : (Ry_inv(t1), Rz(t2), Ry_inv(t3) )
            }

        p = Property()
        r1 = rots[ plane ][0]
        r2 = rots[ plane ][1]
        r3 = rots[ plane ][2]
        p.q = self.q
        p.d = np.einsum('ab,bc,cd,d', r3, r2, r1, self.d )
        p.a = s2ut( np.einsum('ec,fd,ca,db,ai,bj,ij', r3, r3, r2, r2, r1, r1, ut2s(self.a) ) )
        p.Q = s2ut( np.einsum('ec,fd,ca,db,ai,bj,ij', r3, r3, r2, r2, r1, r1, ut2s(self.Q) ) )
        p.b = s2ut( np.einsum('Id,Je,Kf,da,eb,fc,ai,bj,ck,ijk', r3, r3, r3, r2, r2, r2, r1, r1, r1, ut2s(self.b) ) )
        return p

    def rotate( self, t1, t2, t3 ):
        """Rotate all properties by t1, t2, t3
        t1 positive rotation around Z-axis
        t2 negative rotation around Y-axis
        t3 positive rotation around Z-axis
        """
        p = Property()
        r1 = Rz(t1)
        r2 = Ry_inv(t2)
        r3 = Rz(t3)
        p.q = self.q
        p.d = np.einsum('ab,bc,cd,d', r3, r2, r1, self.d )
        p.a = s2ut( np.einsum('ec,fd,ca,db,ai,bj,ij', r3, r3, r2, r2, r1, r1, ut2s(self.a) ) )
        p.Q = s2ut( np.einsum('ec,fd,ca,db,ai,bj,ij', r3, r3, r2, r2, r1, r1, ut2s(self.Q) ) )
        p.b = s2ut( np.einsum('Id,Je,Kf,da,eb,fc,ai,bj,ck,ijk', r3, r3, r3, r2, r2, r2, r1, r1, r1, ut2s(self.b) ) )
        return p

    def transform_by_matrix(self, matrix):
        """docstring for by_matrix"""
        assert matrix.shape == (3,3,)
        p = Property()
        p.q = self.q
        p.d = np.einsum( 'ij,j', matrix, self.d )
        p.a = s2ut(np.einsum( 'ai,bj,ij', matrix, matrix, ut2s(self.a) ))
        p.Q = s2ut(np.einsum( 'ai,bj,ij', matrix, matrix, ut2s(self.Q) ))
        p.b = s2ut(np.einsum( 'ai,bj,ck,ijk', matrix, matrix, matrix, ut2s(self.b) ))
        return p

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



