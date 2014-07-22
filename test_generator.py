#!/usr/bin/env python
#-*- coding: utf-8 -*-

from generator import *
import numpy as np
import argparse

if __name__ == '__main__':

    A = argparse.ArgumentParser( add_help= True)
    A.add_argument( "-l", default = "hyper")
    A.add_argument( "-p", default = "dipole")
    A.add_argument( "-c", default = "Z")

    A.add_argument( "-vary_r"   , default = False , action = 'store_true' ) 
    A.add_argument( "-vary_theta" , default = False , action = 'store_true' )
    A.add_argument( "-vary_tau"   , default = False , action = 'store_true' ) 
    A.add_argument( "-vary_rho1", default = False , action = 'store_true' ) 
    A.add_argument( "-vary_rho2", default = False , action = 'store_true' ) 
    A.add_argument( "-vary_rho3", default = False , action = 'store_true' ) 

    A.add_argument( "-r"   , type = str ) 
    A.add_argument( "-theta" , type = str ) 
    A.add_argument( "-tau"   , type = str ) 
    A.add_argument( "-rho1", type = str ) 
    A.add_argument( "-rho2", type = str ) 
    A.add_argument( "-rho3", type = str ) 
    args = A.parse_args()
    g = Generator()

    opts =  {
             "r"        : {"min":3.15     , "max": 5.00,  "points":1   } ,
             "theta"    : {"max": np.pi   , "min": 0.00 , "points":1   } ,
             "tau"      : {"max": np.pi/2 , "min": 0.00 , "points":1   } ,
             "rho1"     : {"max": np.pi/2   , "min": 0.00 , "points":10   }  ,
             "rho2"     : {"max": np.pi   , "min": 0.00 , "points":10   }  ,
             "rho3"     : {"max": np.pi*2   , "min": 0.00 , "points":10   }  ,
             }
    g.varyParameters( opts )
    g.genMols()
