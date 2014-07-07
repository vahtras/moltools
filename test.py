#!/usr/bin/env python
#-*- coding: utf-8 -*-

import os, re, subprocess
import numpy as np
import  math as m
from particles import *

from water import *
from generator import *
from calculator import *
from template import *
from R import *


def main():
    """Main function for all classes"""
    c = Calculator()
    g = Generator()
    r = R()
    w1 = g.getWater( [0, 0, 0] , 1.0 , 104.5/180 * m.pi )

    refDipole, refAlpha, refBeta = Template().getReference("TIP3P", "HF","PVDZ")

    for i in c.getOutFiles():
        for j in c.getMolFiles():
            if i.rstrip(".out").lstrip("hfqua_") == j.rstrip(".mol"):
                base = j.rstrip(".mol")
                print base

                for k in g.readWaters( j ):
                    d = c.getQmDipole( i )
                    a = c.getQmAlpha( i )
                    b = c.getQmBeta( i )
                    t1, t2, t3 =  k.getEuler()
                    print k.transformDipole( d, t1, t2, t3)
                    print k.transformAlpha( a, t1, t2, t3)
                    print k.transformBeta( b, t1, t2, t3)



if __name__ == '__main__':
    main()


