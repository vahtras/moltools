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



if __name__ == '__main__':
    """Main function for all classes"""
    c = Calculator()
    c.getRelError()

    for i in c.getMatchingOutAndMol():
        r, theta, tau, rho1, rho2, rho3 = i.split('-')
        print "\n" , theta, tau
        print c.Dict.getVal( r, theta, tau, rho1, rho2, rho3 )[2][2]
