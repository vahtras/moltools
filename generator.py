#!/usr/bin/env python
#-*- coding: utf-8 -*-

import re

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt

import math as m
import numpy as np

from water import *
from owndict import *

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
            self.options = { "isAA" : False  }

        self.dic = Dic()

        self.r_oh = 1.83681
        self.theta_hoh = np.pi * 104.5/ 180.0

        self.varyR =     False
        self.varyTheta = False
        self.varyTau =   False
        self.varyRho1 =  False
        self.varyRho2 =  False
        self.varyRho3 =  False

        self.optionsR =     {  "min": 10.00, "max" : 10.00, "points": 1 }
        self.optionsTheta = {  "min": 0.00, "max" : 0.00, "points": 1 }
        self.optionsTau =   {  "min": 0.00, "max" : 0.00, "points": 1 }
        self.optionsRho1 =  {  "min": 0.00, "max" : 0.00, "points": 1 }
        self.optionsRho2 =  {  "min": 0.00, "max" : 0.00, "points": 1 }
        self.optionsRho3 =  {  "min": 0.00, "max" : 0.00, "points": 1 }

        opts = { "r" : {"min": 2.30, "max": 10.00,  "points":50} ,
             "theta" : {"max": np.pi , "min": 0.00 , "points":1 },
             "tau"  : {"max": np.pi/2 , "min": 0.00 , "points":1 },
             "rho1" : {"max": np.pi , "min": 0.00 , "points":1},
             "rho2" : {"max": np.pi , "min": 0.00 , "points":1},
             "rho3" : {"max": np.pi , "min": 0.00 , "points":1},
             }

        self.vary_parameters()

    def gen_mols(self, AA = False,):
        r = np.r_[ self.optionsR[ "min" ] : self.optionsR[ "max" ] : \
                complex( "%sj"%self.optionsR[ "points" ] ) ]
        theta = np.r_[ self.optionsTheta[ "min" ] : self.optionsTheta[ "max" ] : \
                complex( "%sj"%self.optionsTheta[ "points" ] ) ]
        tau = np.r_[ self.optionsTau[ "min" ] : self.optionsTau[ "max" ] : \
                complex( "%sj"%self.optionsTau[ "points" ] ) ]
        rho1 = np.r_[ self.optionsRho1[ "min" ] : self.optionsRho1[ "max" ] : \
                complex( "%sj"%self.optionsRho1[ "points" ] ) ]
        rho2 = np.r_[ self.optionsRho2[ "min" ] : self.optionsRho2[ "max" ] : \
                complex( "%sj"%self.optionsRho2[ "points" ] ) ]
        rho3 = np.r_[ self.optionsRho3[ "min" ] : self.optionsRho3[ "max" ] : \
                complex( "%sj"%self.optionsRho3[ "points" ] ) ]
        for i in r:
            for j in theta:
                for k in tau:
                    for l in rho1:
                        for m in rho2:
                            for n in rho3:
                                w1 = self.get_water( [0, 0, 0], self.r_oh, self.theta_hoh, AA = AA)
                                x, y, z = self.getCartesianFromDegree( i, j, k )
                                w2 = self.get_water( [x,y,z], self.r_oh, self.theta_hoh)
                                w2.rotate( l, m, n )
                                name = ""
                                name += "-".join( map( str, ["%3.2f"%i, "%3.2f"%j, "%3.2f"%k, "%3.2f"%l, "%3.2f"%m, "%3.2f"%n] ) )
                                name += ".mol"
                                self.write_mol( [w1, w2], name )

    def vary_parameters( self, *args ):

        """Given two parameters, for e.g. r and theta, keeps all other static"""

        if args:
            for j in args:
                for i in j:
                    if i == "r":
                        self.varyR = True
                        self.optionsR = j[i]
                    if i == "theta":
                        self.varyTheta = True
                        self.optionsTheta = j[i]
                    if i == "tau":
                        self.varyTau = True
                        self.optionsTau = j[i]
                    if i == "rho1":
                        self.varyRho1 = True
                        self.optionsRho1 = j[i]
                    if i == "rho2":
                        self.varyRho2 = True
                        self.optionsRho2 = j[i]
                    if i == "rho3":
                        self.varyRho3 = True
                        self.optionsRho3 = j[i]

    def get_water(self, origin, r, theta, AA = True ):
        h1 = Atom() ; h2 = Atom() ; o = Atom()
        d = (m.pi/2 - theta/2)
        o.element = "O" ; h1.element = "H" ; h2.element = "H"
        o.x = origin[0] ; o.y = origin[1] ; o.z = origin[2] 
        h1.x = (origin[0] + r * m.cos(d)) ; h1.y = origin[1] ; h1.z = (origin[2] + r*m.sin(d))
        h2.x = (origin[0] - r * m.cos(d)) ; h2.y = origin[1] ; h2.z = (origin[2] + r*m.sin(d))
        w = Water(); w.append( o) ;w.append( h2 ) ;w.append( h1 ) 
        w.theta_hoh = theta
        w.r_oh = r
        w.center = origin
        w.euler1 = 0.00
        w.euler2 = 0.00
        w.euler3 = 0.00

        if AA:
            w.AA = True
            w.h1.AA = True
            w.h2.AA = True
            w.o.AA  = True
        else:
            w.AA = False
            w.h1.AA = False
            w.h2.AA = False
            w.o.AA  = False
        return w


    def readWaters(self, fname):
        """From file with name fname, return a list of all waters encountered"""
#If the file is plain xyz file

        atoms = []
        if fname.endswith( ".xyz" ) or fname.endswith(".mol"):
            pat_xyz = re.compile(r'^\s*(\w+)\s+(-*\d*.+\d+)\s+(-*\d*.+\d+)\s+(-*\d*.+\d+) *$')
            for i in open( fname ).readlines():
                if pat_xyz.match(i):
                    f = pat_xyz.match(i).groups()
                    tmpAtom = Atom()
                    tmpAtom.AA = True
                    tmpAtom.x = float(f[1])
                    tmpAtom.y = float(f[2])
                    tmpAtom.z = float(f[3])
                    tmpAtom.element = f[0][0]
                    atoms.append( tmpAtom )

        elif fname.endswith( ".pdb" ):
            pat1 = re.compile(r'^(ATOM|HETATM)')
            for i in open( fname ).readlines():
                if pat1.search(i):
                    #Ignore charge centers for polarizable water models
                    if ( i[11:16].strip() == "SW") or (i[11:16] == "DW"):
                        continue
                    tmpAtom = Atom(i[11:16].strip()[0], \
                            float(i[30:38].strip()), \
                            float(i[38:46].strip()), \
                            float(i[46:54].strip()), \
                            int(i[22:26].strip()) )

                    if fnameAAorAU == "AU":
                        if args.opAAorAU == "AA":
                            tmpAtom.toAA()
                    elif fnameAAorAU == "AA":
                        if args.opAAorAU == "AU":
                            tmpAtom.toAU()
                    atoms.append( tmpAtom )
        elif fname.endswith( ".out" ):
            pat_xyz = re.compile(r'^(\w+)\s+(-*\d*.+\d+)\s+(-*\d*.+\d+)\s+(-*\d*.+\d+) *$')
            for i in open( fname ).readlines():
                if pat_xyz.match(i):
                    f = pat_xyz.match(i).groups()
                    tmpAtom = Atom(f[0][0], float(f[1]), float(f[2]), float(f[3]), 0)
                    if fnameAAorAU == "AU":
                        if args.opAAorAU == "AA":
                            tmpAtom.toAA()
                    elif fnameAAorAU == "AA":
                        if args.opAAorAU == "AU":
                            tmpAtom.toAU()
                    atoms.append( tmpAtom )
#loop over oxygen and hydrogen and if they are closer than 1 A add them to a water
        waters = []
        cnt = 1

        if fname.endswith( ".xyz" ) or fname.endswith(".mol"):
            for i in atoms:
                if i.element == "H":
                    continue
                if i.inWater:
                    continue
                tmp = Water(  )
                i.inWater = True
                tmp.append( i )
                for j in atoms:
                    if j.element == "O":
                        continue
                    if j.inWater:
                        continue
#If in cartesian:
                    if j.AA:
                        if i.distToAtom(j) < 1.1:
                            tmp.append ( j )
                            j.inWater = True
                    else:
                        if i.distToAtom(j) < 1.1/a0:
                            tmp.append ( j )
                            j.inWater = True
                tmp.number = cnt
                cnt += 1
                waters.append( tmp )
        elif fname.endswith( ".pdb" ):
#Find out the size of the box encompassing all atoms
            xmin = 10000.0; ymin = 10000.0; zmin = 10000.0; 
            xmax = -10000.0; ymax = -10000.0; zmax = -10000.0; 
            for i in atoms:
                if i.x < xmin:
                    xmin = i.x
                if i.y < ymin:
                    ymin = i.y
                if i.z < zmin:
                    zmin = i.z
                if i.x > xmax:
                    xmax = i.x
                if i.y > ymax:
                    ymax = i.y
                if i.z > zmax:
                    zmax = i.z
            center = np.array([ xmax - xmin, ymax -ymin, zmax- zmin]) /2.0
            wlist = []
            for i in atoms:
                if i.element != "O":
                    continue
                tmp = Water()
                i.inWater= True
#__Water__.append() method will update the waters residue number and center coordinate
#When all atoms are there
#Right now NOT center-of-mass
                tmp.append(i)
                for j in atoms:
                    if j.element != "H":
                        continue
                    if j.inWater:
                        continue
#1.05 because sometimes spc water lengths can be over 1.01
                        
                    if args.opAAorAU == "AA":
                        if i.dist(j) <= 1.05:
                            j.inWater = True
                            tmp.append( j )
                            if len(tmp.atomlist) == 3:
                                break
                    elif args.opAAorAU == "AU":
                        if i.dist(j) <= 1.05/a0:
                            j.inWater = True
                            tmp.append( j )
                            if len(tmp.atomlist) == 3:
                                break
                wlist.append( tmp )
            wlist.sort( key = lambda x: x.distToPoint( center ))
            center_water = wlist[0]
            cent_wlist = wlist[1:]
            cent_wlist.sort( key= lambda x: x.distToWater( center_water) )
            waters = [center_water] + cent_wlist[ 0:args.waters - 1 ]
        elif fname.endswith( ".out" ):
            for i in atoms:
                if i.element == "H":
                    continue
                if i.inWater:
                    continue
                tmp = Water(  )
                i.inWater = True
                tmp.append( i )
                for j in atoms:
                    if j.element == "O":
                        continue
                    if j.inWater:
                        continue
#If in cartesian:
                    if i.AA:
                        if i.dist(j) < 1.0:
                            tmp.append ( j )
                            j.inWater = True
                    else:
                        if i.dist(j) < 1.0/a0:
                            tmp.append ( j )
                            j.inWater = True
                tmp.number = cnt
                cnt += 1
                waters.append( tmp )
        return waters

    def write_mol(self, wlist, name = "tmp.mol" ):
        f_ = open (name, 'w')
        f_.write("ATOMBASIS\n\n\nAtomtypes=2 Charge=0 Nosymm\n")
        f_.write("Charge=1.0 Atoms=4 Basis=cc-pVDZ\n")
        for i in wlist:
            for j in i:
                if j.element != "H":
                    continue
                f_.write( "%s %.5f %.5f %.5f\n" %(j.element, j.x, j.y, j.z ) )
        f_.write("Charge=8.0 Atoms=2 Basis=cc-pVDZ\n")
        for i in wlist:
            for j in i:
                if j.element != "O":
                    continue
                f_.write( "%s %.5f %.5f %.5f\n" %(j.element, j.x, j.y, j.z ) )
        f_.close()

    def getCartesianFromDegree(self, r, theta, tau):
        return r* m.sin( m.pi*theta/180.0 )*m.cos( m.pi*tau/180.0) \
               , r* m.sin( m.pi*theta/180.0 )*m.sin( m.pi*tau/180.0)  \
               , r* m.cos( m.pi*theta/180.0 )

if __name__ == '__main__':
    g = Generator()
    w1 = g.get_water( [0,1,1], 1.0, 101.4 )
    g.write_mol( [w1] )
    #print w1.get_mol()
    #w1.plotWater()


