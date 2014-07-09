#!/usr/bin/env python
#-*- coding: utf-8 -*-

from generator import *
import numpy as np

if __name__ == '__main__':
    g = Generator()

    opts = { "r" : {"min":2.30, "max":5.00,  "points":10} ,
             "theta" : {"max": np.pi , "min": 0.00 , "points":10},
             "tau" : {"max": np.pi/2 , "min": 0.00 , "points":10},
             "rho1" : {"max": np.pi , "min": 0.00 , "points":1},
             "rho2" : {"max": np.pi , "min": 0.00 , "points":1},
             "rho3" : {"max": np.pi , "min": 0.00 , "points":1},
             }

    g.varyParameters( opts )
    g.genMols()
