#!/usr/bin/env python
#-*- coding: utf-8 -*-

from calculator import *
from water import *
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
        var = "r"
    if args.vary_theta:
        c.opts[ "theta" ] = { "vary" : True }
        var = "theta"
    if args.vary_tau:
        c.opts[ "tau" ] = { "vary" : True }
        var = "tau"
    if args.vary_rho1:
        c.opts[ "rho1" ] = { "vary" : True }
        var = "rho1"
    if args.vary_rho2:
        c.opts[ "rho2" ] = { "vary" : True }
        var = "rho2"
    if args.vary_rho3:
        c.opts[ "rho3" ] = { "vary" : True }
        var = "rho3"


    c.getRelError()
    c.writeLog( level = args.l , prop = args.p,  component = args.c, variable = var )


