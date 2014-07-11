#!/usr/bin/env python
#-*- coding: utf-8 -*-

from calculator import *
from water import *
import argparse

if __name__ == '__main__':

    A = argparse.ArgumentParser( add_help= True)
    A.add_argument( "-l", default = "polar")
    A.add_argument( "-p", default = "dipole")
    A.add_argument( "-c", default = "Z")

    args = A.parse_args()

    c = Calculator()
    c.getRelError()
    c.writeLog( level = args.l , prop = args.p,  component = args.c )
