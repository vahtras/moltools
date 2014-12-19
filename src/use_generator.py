#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""
Purpose
------------------

Automate creating water dimers for DALTON input in order to study the properties dependence on various parameters.


Usage
------------------

This script creates a water dimer with one water with its oxygen in origo and the other water with its oxygen at a separation :math:`r` from the first oxygen with spheciral coordinate angles :math:`\\tau` and :math:`\\theta`.

.. code:: bash

    $ use_generator.py -param -r_min 5 -r_max 10 -r_points 10

Will create 10 .mol files where the distance between the water molecules goes from 5 to 10 atomic units.

To analyze the results from the output of the files generated, see :ref:`use_calculator`.

-----------

Further parameters that can be varied are :math:`\\tau`, :math:`\\theta`, :math:`\\rho_1`, :math:`\\rho_2` and :math:`\\rho_3`.

:math:`\\rho_1` - :math:`\\rho_3` are described in :ref:`Water`.


By default :math:`r=5.0` in atomic units, and :math:`\\tau` = :math:`\\theta` = :math:`\\rho_{1}` = :math:`\\rho_2` = :math:`\\rho_3` = 0

.. note:
    1 atonic unit = 0.529 Å.

"""

import argparse, re, fractions

import numpy as np
import math as m

from template import Template
from molecules import Molecule, Water, Methanol, Atom, Property, Cluster

a0 = 0.52917721092

class Generator( dict ):
    """
    Used to create molecules, write dalton .mol files 
    using -param for study with use_calculator.py

    water currently implemented only

    plans to implement methanol

    """
    def __init__(self, *args, **kwargs):

#This waater is TIP3P model, generalzie later
        self[ ("water","a_hoh", "degree") ] = 104.5
        self[ ("water","r_oh", "AA") ] = 0.9572

        self[ ("methanol", "r_oh", "AA" ) ] = 0.967
        self[ ("methanol", "r_co", "AA" ) ] = 1.428
        self[ ("methanol", "r_ch", "AA" ) ] = 1.098

        self[ ("methanol", "a_coh", "degree" ) ] = 107.16
        self[ ("methanol", "a_hch", "degree" ) ] = 109.6
        self[ ("methanol", "a_hco", "degree" ) ] = 109.342

        self[ ("methanol", "d_hcoh", "h4", "degree" ) ] =  60.0
        self[ ("methanol", "d_hcoh", "h5", "degree" ) ] = -60.0
        self[ ("methanol", "d_hcoh", "h6", "degree" ) ] =  180.0

        
#Default options for water
        for val in ["r", "tau", "theta", "rho1", "rho2", "rho3", ]:
            self[ ( val, 'min') ]    = 0.0
            self[ ( val, 'max') ]    = 0.0
            self[ ( val, 'points') ] = 1
        self[ ( 'r', 'min') ]    = 5.0
        self[ ( 'r', 'max') ]    = 10.0
        self[ ( 'r', 'points') ] = 1

# Set by default all parameters to False
        for val in ["r", "tau", "theta", "rho1", "rho2", "rho3", ]:
            self[ ( val, "active" ) ]  = False

    @staticmethod
    def get_hfqua_dal():
        return """**DALTON INPUT
.RUN RESPONSE
.DIRECT
.PARALLELL
**WAVE FUNCTION
.HF
**RESPONSE
*QUADRATIC
.DIPLEN
**END OF DALTON INPUT"""

    def gen_mols_param(self, mol = "water", 
            basis = ["ano-1 2 1", "ano-1 3 2 1"],
            AA = True):
        r = np.linspace( self[ ('r', 'min')] , self[ ('r', 'max')] ,
            self[ ('r', 'points' ) ]  )
        tau = np.linspace( self[ ('tau', 'min')] , self[ ('tau', 'max')] ,
            self[ ('tau', 'points' ) ] )
        theta = np.linspace( self[ ('theta', 'min')] , self[ ('theta', 'max')] ,
            self[ ('theta', 'points' )  ] )
        rho1 = np.linspace( self[ ('rho1', 'min')], self[ ('rho1', 'max')],
            self[ ('rho1', 'points' )  ] )
        rho2 = np.linspace( self[ ('rho2', 'min')], self[ ('rho2', 'max')],
            self[ ('rho2', 'points' )  ] )
        rho3 = np.linspace( self[ ('rho3', 'min')], self[ ('rho3', 'max')],
            self[ ('rho3', 'points' )  ] )

        r_oh = self[ ("water","r_oh", "AA") ]
        a_hoh = np.pi * self[ ("water", "a_hoh", "degree" )] / 180.0

        for i in r:
            for j in tau:
                for k in theta:
                    for l in rho1:
                        for m in rho2:
                            for n in rho3:
                                c= Cluster()
                                w1 = self.get_mol( [0, 0, 0], mol , AA = AA)
                                c.append( w1, in_qm = True )
                                x, y, z = self.polar_to_cartesian( i, j, k )
                                w2 = self.get_mol( [x,y,z], mol, AA = AA)
                                w2.rotate( l, m, n )

                                c.append( w2, in_qm = True )
                                name = ""
                                name += "-".join( map( str, ["%3.2f"%i, "%3.2f"%j, "%3.2f"%k, "%3.2f"%l, "%3.2f"%m, "%3.2f"%n] ) )
                                name += ".mol"

                                m = c.get_qm_mol_string( AA = AA,
                                        basis = tuple(basis),
                                        )
                                f_ = open(name, 'w')
                                f_.write( m )
        return 0

    def vary_parameters( self, opts ):
        """Given two parameters, e.g. r and theta, keeps all other static
        param_list should be list of strings of parameters
        ["r":{"min": 2, "max":5, "points": 10}, "rho1" , ... ]

        Has sane defaults, but can be overrided by passing arguments to 
        main program as:

        -r_min 5
        -r_max 10
        -r_points 10

        Which overrides defaults 

        """
        for val in opts:
            self[ (val, 'active') ] = True
            self[ (val, 'min') ] = opts[val][ "min" ]
            self[ (val, 'max') ] = opts[val][ "max" ]
            self[ (val, 'points') ] = opts[val][ "points" ]

    def get_mol( self, center = [0,0,0], mol = "water", AA = False ):
        """return molecule in center, all molecules have different definition
        of euler angles

        for water place O in origo
        for methanol place C=O bond in origo
        
        """

        if mol == "water":
#Geometrical parameters
            r_oh = self[ ("water","r_oh", "AA") ]
            a_hoh = self[ ("water","a_hoh","degree") ]
            if not AA:
                r_oh = r_oh / a0

            d = (90 - a_hoh/2 ) * np.pi / 180


            xo = center[0]
            yo = center[1]
            zo = center[2] 

            xh1 = (center[0] + r_oh * np.cos(d))
            yh1 =  center[1] 
            zh1 = (center[2] + r_oh* np.sin(d))

            xh2 = (center[0] - r_oh * np.cos(d)) 
            yh2 = center[1] 
            zh2 = (center[2] + r_oh* np.sin(d))

            h1 = Atom( **{ "AA" : AA,
                "x" : xh1,
                "y" : yh1,
                "z" : zh1,
                "element" : "H"} )
            h2 = Atom( **{ "AA" : AA,
                "x" : xh2,
                "y" : yh2,
                "z" : zh2,
                "element" : "H"} )
            o = Atom( **{ "AA" : AA,
                "x" : xo,
                "y" : yo,
                "z" : zo,
                "element" : "O"} )

            w = Water( AA = AA)
            w.append( o )
            w.append( h1 )
            w.append( h2 )
            
            return w

        elif mol == "methanol":

            r_co = self[ ("methanol", "r_co", "AA" )]
            r_oh = self[ ("methanol", "r_oh", "AA" )]
            r_ch = self[ ("methanol", "r_ch", "AA" )]

            a_coh = self[ ("methanol", "a_coh", "degree" ) ]
            #a_hch = self[ ("methanol", "a_hch", "degree" ) ]
            a_hco = self[ ("methanol", "a_hco", "degree" ) ]

            a_coh *= np.pi / 180
            a_hco *= np.pi / 180

            d_hcoh_4 = self[ ("methanol", "d_hcoh", "h4", "degree" ) ]
            d_hcoh_4 *= np.pi / 180
            d_hcoh_5 = self[ ("methanol", "d_hcoh", "h5", "degree" ) ]
            d_hcoh_5 *= np.pi / 180
            d_hcoh_6 = self[ ("methanol", "d_hcoh", "h6", "degree" ) ]
            d_hcoh_6 *= np.pi / 180

            if not AA:
                r_co, r_oh, r_ch = r_co/a0, r_oh/a0, r_ch/a0

            c1 = Atom( **{"x":0, "y":0, "z":-r_co/2, "AA": AA, "element":"C" } )
            o2 = Atom( **{"x":0, "y":0, "z": r_co/2, "AA": AA, "element":"O" } )

            h3 = Atom( **{"x":r_oh*np.cos( a_coh-np.pi/2),
                "y":0,
                "z":r_oh*np.sin( a_coh-np.pi/2) + r_co/2,
                "AA": AA, "element":"H" } )

            h4 = Atom( **{"x": r_ch*np.sin( a_hco ) * np.cos( d_hcoh_4 ),
                "y": r_ch*np.sin( a_hco) * np.sin( d_hcoh_4 ),
                "z": r_ch*np.cos( a_hco) - r_co/2 ,
                "AA": AA, "element":"H" } )
            h5 = Atom( **{"x": r_ch*np.sin( a_hco ) * np.cos( d_hcoh_5 ),
                "y": r_ch*np.sin( a_hco) * np.sin( d_hcoh_5 ),
                "z": r_ch*np.cos( a_hco) - r_co/2 ,
                "AA": AA, "element":"H" } )
            h6 = Atom( **{"x": r_ch*np.sin( a_hco ) * np.cos( d_hcoh_6 ),
                "y": r_ch*np.sin( a_hco) * np.sin( d_hcoh_6 ),
                "z": r_ch*np.cos( a_hco) - r_co/2 ,
                "AA": AA, "element":"H" } )

            m = Methanol()
            m.append(c1)
            m.append(o2)
            m.append(h3)
            m.append(h4)
            m.append(h5)
            m.append(h6)

            return m

    def polar_to_cartesian(self, r, tau, theta):
        x, y, z = r* m.sin( theta )*m.cos( tau ) \
               , r* m.sin(  theta )*m.sin( tau )  \
               , r* m.cos(  theta ) 

        return x , y , z

    def build_pna( self,  xyz = "tmp.xyz", waters = 0,
            min_r = 2.0,
            mult_r = 10,
            seed = 111 ):
        pna = Molecule.from_xyz( xyz )
        freqs = [ "0.0", "0.0238927", "0.0428227", "0.0773571" ] 

        np.random.seed( seed )

        c = Cluster()
        c.append(pna, in_qm = True)
        cnt = 0
        while cnt < waters:
# Random rotation angles
            t1 = np.random.uniform( 0, np.pi/2 )
            t2 = np.random.uniform( 0, np.pi   )
            t3 = np.random.uniform( 0, np.pi/2 )

# random length, rho and tau 
            r =  np.random.uniform( min_r , min_r * mult_r)
            tau =  np.random.uniform( 0, np.pi*2)
            theta =  np.random.uniform( 0,np.pi)

            center = self.polar_to_cartesian( r, tau, theta )

            wat = self.get_mol( center = pna.com + center,
                    mol = "water")

            wat.rotate( t1, t2, t3 )
            wat.res_id = cnt

            if c.mol_too_close( wat ):
                continue

#We are satisfied with this position, add properties to the water, and rotate them according to t1, t2, t3 so they match the water orientation
            c.append( wat, in_mm = True )
            cnt += 1

        for f_mm in freqs:
            for dist in ["nodist", "dist"]:
                for wat in [ m for m in c if m.in_mm ]:
                    t1, t2, t3 =  wat.get_euler()
                    kwargs_dict = Template().get( *("TIP3P", "HF", "ANOPVDZ",
                        dist == "dist",f_mm ) )
                    for at in wat:
                        Property.add_prop_from_template( at, kwargs_dict )
                    Property.transform_ut_properties( wat.h1.Property, t1,t2,t3 )
                    Property.transform_ut_properties( wat.h2.Property, t1,t2,t3 )
                    Property.transform_ut_properties( wat.o.Property,  t1,t2,t3 )
#Write out QM and MM region separately with properties
                open("pna.mol" ,'w').write(c.get_qm_mol_string(
                    basis= ("ano-1 2 1", "ano-1 3 2 1"),
                    AA = True))
                open("%dmm_%s_%s.pot" %(waters, f_mm, dist ),'w').write(c.get_qmmm_pot_string( in_AA = True ))
                open("tmp.xyz", 'w').write( c.get_xyz_string() )

if __name__ == '__main__':

    A = argparse.ArgumentParser( add_help= True)

#Related to generating molecules with specified parameters
# Right now 2 waters implemented with 6 parameters

    A.add_argument( "-param", action = 'store_true', default = False )
    A.add_argument( "-param_mol", type = str, default = 'water' )



    A.add_argument( "-r_min"     ,   type = float , default = 3.00  ) 
    A.add_argument( "-theta_min" ,   type = float , default = 0.00  ) 
    A.add_argument( "-tau_min"   ,   type = float , default = 0.00  ) 
    A.add_argument( "-rho1_min"  ,   type = float , default = 0.00  ) 
    A.add_argument( "-rho2_min"  ,   type = float , default = 0.00  ) 
    A.add_argument( "-rho3_min"  ,   type = float , default = 0.00  ) 

    A.add_argument( "-r_max"     ,   type = float , default = 10.0 ) 
    A.add_argument( "-theta_max" ,   type = float , default = np.pi/2    ) 
    A.add_argument( "-tau_max"   ,   type = float , default = np.pi/2    ) 
    A.add_argument( "-rho1_max"  ,   type = float , default = np.pi/2    ) 
    A.add_argument( "-rho2_max"  ,   type = float , default = np.pi      ) 
    A.add_argument( "-rho3_max"  ,   type = float , default = np.pi/2    ) 

    A.add_argument( "-r_points"     ,   type = int , default = 1) 
    A.add_argument( "-theta_points" ,   type = int , default = 1) 
    A.add_argument( "-tau_points"   ,   type = int , default = 1) 
    A.add_argument( "-rho1_points"  ,   type = int , default = 1) 
    A.add_argument( "-rho2_points"  ,   type = int , default = 1) 
    A.add_argument( "-rho3_points"  ,   type = int , default = 1) 

#Related to generating/getting one water
    A.add_argument( "-get_mol", 
            type = str, choices = ["water",
                "methanol" , "ethane"] )
    A.add_argument( "-basis", type = str, nargs = "*",
            default =["ano-1 2 1", "ano-1 3 2 1" ] )
    A.add_argument( "-AA" ,  default = False, action = 'store_true' )

#########################
#      PNA RELATED
#########################

    A.add_argument( '-seed', type = int, default = 111 )
    A.add_argument( '-pna', action = 'store_true', default = False )
    A.add_argument( '-pna_waters', type = int, default = 10 )
    A.add_argument( '-pna_min_r', type = float, default = 0 )
    A.add_argument( '-pna_mult_r', type = float, default = 10 )

########################################################################

    args = A.parse_args()

#Obtain variable parameters from program call and set them
    g = Generator( )
    opts =  {
       "r"    :{"min":args.r_min,    "max": args.r_max,    "points":args.r_points,    },
       "tau"  :{"min":args.tau_min,  "max": args.tau_max,  "points":args.tau_points,  },
       "theta":{"min":args.theta_min,"max": args.theta_max,"points":args.theta_points,},
       "rho1" :{"min":args.rho1_min, "max": args.rho1_max, "points":args.rho1_points, },
       "rho2" :{"min":args.rho2_min, "max": args.rho2_max, "points":args.rho2_points, },
       "rho3" :{"min":args.rho3_min, "max": args.rho3_max, "points":args.rho3_points, },
    }
    g.vary_parameters( opts )


    if args.param:
        g.gen_mols_param( 
                mol = args.param_mol ,
                basis = args.basis,
                AA = False )

    if args.get_mol:
        mol = g.get_mol( center = [0, 0, 0 ], mol = args.get_mol , AA = args.AA )

    if args.pna:
        g.build_pna( xyz = "pna_opt.xyz", 
                waters = args.pna_waters,
                min_r = args.pna_min_r,
                mult_r = args.pna_mult_r,
                seed = args.seed )
                
 
