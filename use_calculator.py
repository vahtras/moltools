#!/usr/bin/env python
#-*- coding: utf-8 -*-

import os, re, subprocess
import numpy as np
import math as m
import argparse

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt

from particles import *
from gaussian import *

from template import Template
from water import Atom, Water
from dic import Dic
from calculator import Calculator
from generator import Generator
from xms import Xms

a0 = 0.52917721092
charge_dic = {"H": 1.0, "C": 6.0, "N": 7.0, "O": 8.0, "S": 16.0}
mass_dict = {"H": 1.008,  "C": 6.0, "N": 7.0, "O": 15.999, "S": 16.0}
index_dict = { 
        "static" : {"dipole":{"X": [0,0] , "Y": [0,1], "Z": [0,2]}},

        "polar" : {"dipole":{"X": [1,0] , "Y": [1,1], "Z": [1,2]},
                   "alpha":{"X": [3,0] , "Y": [3,1], "Z": [3,2]} },
        "hyper" : { "dipole":{"X": [2,0], "Y": [2,1], "Z": [2,2]},
                   "alpha":{"X": [4,0] , "Y": [4,1], "Z": [4,2] } ,
                   "beta":{"X": [5,0] , "Y": [5,1], "Z": [5,2]} }}
line_thick_dict = { 
        "static" : {"dipole":{"X": 0, "Y": 0, "Z": 0}},

        "polar" : {"dipole":{"X": 2 , "Y": 2, "Z": 2},
                   "alpha":{"X":  2 , "Y": 2, "Z": 2} },
        "hyper" : { "dipole":{"X": 3 , "Y": 3, "Z": 3 },
                   "alpha":{  "X": 3 , "Y": 3, "Z": 3 },
                   "beta":{   "X": 3 , "Y": 3, "Z": 3 }}}
line_style_dict = { 
        "static" : {"dipole":{"X": 1, "Y": 1, "Z": 1}},

        "polar" : {"dipole":{"X": 2 , "Y": 2, "Z": 2},
                   "alpha":{"X":  2 , "Y": 2, "Z": 2} },
        "hyper" : { "dipole":{"X": 3, "Y": 3, "Z": 3 },
                   "alpha":{  "X": 3, "Y": 3, "Z": 3 },
                   "beta":{   "X": 3, "Y": 3, "Z": 3 }}}

if __name__ == '__main__':

    A = argparse.ArgumentParser( add_help= True)

    A.add_argument( "-model"  , type = str, default = "gaussian", choices = [ "pointdipole", "quadrupole", "gaussian"] )
    A.add_argument( "-Rp"  , type = str, default = "0.00001" )
    A.add_argument( "-Rq"  , type = str, default = "0.00001" )
    A.add_argument( "-param"  , default = True )
    A.add_argument( "-dist"   , action = 'store_true', default = False )
    A.add_argument( "-no_subtitle"  , default = False, action = 'store_true' )

    A.add_argument( "-qm"   , type = str,  default = "hfqua", dest = "qm_method" )
    A.add_argument( "-rel", action = 'store_true' , default = False)

    A.add_argument( "-l", nargs = "*", default = ["static"] )
    A.add_argument( "-p", nargs = "*", default = ["dipole"] )
    A.add_argument( "-c", nargs = "*", default = ["Z"] )

    A.add_argument( "-max_l", type = int, default = 2, choices = [1 , 2] )
    #A.add_argument( "-pol", type = int, default = 22, choices = [1 , 2, 22] )
    #A.add_argument( "-hyper", type = int , default = 1, choices = [1] )


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
    c = Calculator()


    if args.r:
        c.opts[ "r" ] = { "constant" : args.r }
    if args.theta:
        c.opts[ "theta" ] = { "constant" : args.theta }
    if args.tau:
        c.opts[ "tau" ] = { "constant" : args.tau }
    if args.rho1:
        c.opts[ "rho1" ] = { "constant" : args.rho1 }
    if args.rho2:
        c.opts[ "rho2" ] = { "constant" : args.rho2 }
    if args.rho3:
        c.opts[ "rho3" ] = { "constant" : args.rho3 }

    if args.vary_r:
        c.opts[ "r" ] = { "vary" : True }
        args.var = "r"
    if args.vary_theta:
        c.opts[ "theta" ] = { "vary" : True }
        args.var = "theta"
    if args.vary_tau:
        c.opts[ "tau" ] = { "vary" : True }
        args.var = "tau"
    if args.vary_rho1:
        c.opts[ "rho1" ] = { "vary" : True }
        args.var = "rho1"
    if args.vary_rho2:
        c.opts[ "rho2" ] = { "vary" : True }
        args.var = "rho2"
    if args.vary_rho3:
        c.opts[ "rho3" ] = { "vary" : True }
        args.var = "rho3"

    if args.rel:
        c.get_rel_error( args )
        string = c.get_xvg_string_rel( args )
    else:
        c.get_abs_value( args )
        string = c.get_xvg_string_abs( args )
    open('tmp.xvg', 'w').write( string )
