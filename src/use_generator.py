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
from molecules import Generator, Molecule, Water, Methanol, Atom, Property, Cluster

a0 = 0.52917721092
if __name__ == '__main__':

    A = argparse.ArgumentParser( add_help= True)

# Related to generating molecules with specified parameters
# Right now 2 waters implemented with 6 parameters
# TWO MOLS GEN RELATED

    A.add_argument( "-two_mols_gen", action = 'store_true', default = False,
            help = "Generate two molecule at certain distance/parameters relative to each other")
    A.add_argument( "-two_mols_mol", type = str, default = 'water',
            help = "Which specific molecules to generate")
    A.add_argument( "-two_mols_model", type = str, default = 'tip3p',
            help = "Which specific model of the molecule to generate")
    A.add_argument( "-worst", action = 'store_true', default = False,)

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


# ONE MOL GEN RELATED
    A.add_argument( "-one_mol_gen", action = 'store_true', default = False,)
    A.add_argument( "-one_mol_mol", type = str, default = 'water',)
    A.add_argument( "-one_mol_model", type = str, default = 'tip3p',)
    A.add_argument( "-r_oh_dev", type = float, default = 5,)
    A.add_argument( "-theta_hoh_dev", type = float, default = 5,)

    A.add_argument( "-r_oh_points", type = int, default = 1,)
    A.add_argument( "-theta_hoh_points", type = int, default = 1,)

#Related to generating/getting one water
    A.add_argument( "-get_mol", 
            type = str, choices = ["water",
                "methanol" , "ethane"] )
    A.add_argument( "-model", type = str, default = 'tip3p' )
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
       "r_oh_dev" :{"min":0, "max": args.r_oh_dev, "points": args.r_oh_points },
       "theta_hoh_dev" :{"min":0, "max": args.theta_hoh_dev, "points": args.theta_hoh_points },
    }
    g.vary_parameters( opts )


    if args.two_mols_gen:
        g.gen_mols_param( 
                mol = args.two_mols_mol ,
                model = args.two_mols_model,
                basis = args.basis,
                AA = False,
                worst = args.worst )
    if args.one_mol_gen:
        g.one_mol_gen(
            mol = args.one_mol_mol,
            model = args.one_mol_model,
            )

    if args.get_mol:
        mol = g.get_mol( center = [0, 0, 0 ], mol = args.get_mol , AA = args.AA )

    if args.pna:
        g.build_pna( xyz = "pna_opt.xyz", 
                waters = args.pna_waters,
                min_r = args.pna_min_r,
                mult_r = args.pna_mult_r,
                seed = args.seed )
                
 
