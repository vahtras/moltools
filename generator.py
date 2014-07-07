#!/usr/bin/env python
#-*- coding: utf-8 -*-

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt
from water import *

import numpy as np
import math as m

class Generator:
    """
    General class to manipulate water molecules, write and read dalton .mol files and .out files
    Also to plot water molecules for testing rotations.

    And to write xmgrace data from dalton output
    """

    def __init__(self, *args, **kwargs):
        if kwargs is not None:
            self.options = kwargs
        else:
            self.options = { "isAA" : True  }

    def getWater(self, origin, r, theta ):
        h1 = Atom() ; h2 = Atom() ; o = Atom()
        d = (m.pi/2 - theta/2)
        o.element = "O" ; h1.element = "H" ; h2.element = "H"
        o.x = origin[0] ; o.y = origin[1] ; o.z = origin[2] 
        h1.x = (origin[0] + r * m.cos(d)) ; h1.y = origin[1] ; h1.z = (origin[2] + r*m.sin(d))
        h2.x = (origin[0] - r * m.cos(d)) ; h2.y = origin[1] ; h2.z = (origin[2] + r*m.sin(d))
        w = Water(); w.addAtom( o) ;w.addAtom( h2 ) ;w.addAtom( h1 ) 
        w.theta_hoh = theta
        w.r_oh = r
        w.center = origin
        w.euler1 = 0.00
        w.euler2 = 0.00
        w.euler3 = 0.00

        if "isAU" in self.options:
            w.AA = False
            w.h1.AA = False
            w.h2.AA = False
            w.o.AA  = False
        else:
            w.AA = True
            w.h1.AA = True
            w.h2.AA = True
            w.o.AA  = True
        return w

    def plotWater(self, water ):
#Plot water molecule in green and  nice xyz axis
        O1, H1, H2 = water.o, water.h1, water.h2
        fig = plt.figure()
        dip = water.getDipole()
        ax = fig.add_subplot(111, projection='3d' )
        ax.plot( [0, 1, 0, 0, 0, 0], [0, 0,0,1,0,0], [0,0,0,0,0,1] )
        ax.plot( [O1.x,O1.x + dip[0] ] ,[ O1.y,O1.y+dip[1]],[O1.z,O1.z+dip[2]] ,'-',color="black")
        ax.scatter( [H1.x], [ H1.y] ,[ H1.z], s=25, color='red')
        ax.scatter( [H2.x], [ H2.y] ,[ H2.z], s=25, color='red')
        ax.scatter( [O1.x], [ O1.y] ,[ O1.z], s=50, color='blue')
        ax.set_zlim3d( -5,5)
        plt.xlim(-5,5)
        plt.ylim(-5,5)
        plt.show()
    def writeMol(self, wlist ):
        b = ""
        for i in wlist:
            if i.resId == 1:
                continue
            b += "-".join( map( str, [i.r, "%3.2f"%i.theta, "%3.2f"%i.tau, i.euler1, i.euler2, i.euler3] ) )
            b += ".mol"
        f_ = open (b, 'w')
        f_.write("ATOMBASIS\n\n\nAtomtypes=2 Charge=0 Angstrom Nosymm\n")
        f_.write("Charge=1.0 Atoms=4 Basis=cc-pVDZ\n")
        for i in wlist:
            for j in i.atomlist:
                if j.element != "H":
                    continue
                f_.write( "%s %.5f %.5f %.5f\n" %(j.element, j.x, j.y, j.z ) )
        f_.write("Charge=8.0 Atoms=2 Basis=cc-pVDZ\n")
        for i in wlist:
            for j in i.atomlist:
                if j.element != "O":
                    continue
                f_.write( "%s %.5f %.5f %.5f\n" %(j.element, j.x, j.y, j.z ) )
        f_.close()


if __name__ == '__main__':
    g = Generator()
    w1 = g.getWater( [1,1,1], 2.3, 101.4 )
    g.plotWater( w1 )
