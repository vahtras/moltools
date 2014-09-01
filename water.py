#!/usr/bin/env python
#-*- coding: utf-8 -*-

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt

import numpy as np
import math as m

a0 = 0.52917721092
charge_dic = {"H": 1.0, "C": 6.0, "N": 7.0, "O": 8.0, "S": 16.0}
mass_dict = {"H": 1.008,  "C": 6.0, "N": 7.0, "O": 15.999, "S": 16.0}

class Dic:
    def __init__(self):
        rho3 = { 0.00 : [] }
        rho2 = { 0.00 : rho3 }
        rho1 = { 0.00 : rho2 }
        tau = { 0.00 : rho1}
        theta = { 0.00 : tau }
        r = { 0.00 : theta }
        self.dic = r
    def getVal(self, r, theta, tau, rho1, rho2, rho3):
        if self.dic.has_key( r ):
            if self.dic[r].has_key( theta ):
                if self.dic[r][theta].has_key(tau):
                    if self.dic[r][theta][tau].has_key(rho1):
                        if self.dic[r][theta][tau][rho1].has_key(rho2):
                            if self.dic[r][theta][tau][rho1][rho2].has_key(rho3):
                                return self.dic[r][theta][tau][rho1][rho2][rho3]
    def setVal(self, r, theta, tau, rho1, rho2, rho3, val ): #, rho1, rho2, rho3, val):
        tmpr = r
        tmptheta = theta
        tmptau = tau
        tmprho1 = rho1
        tmprho2 = rho2
        tmprho3 = rho3
        if self.dic.has_key( r ):
            if self.dic[r].has_key( theta ):
                if self.dic[r][theta].has_key(tau):
                    if self.dic[r][theta][tau].has_key(rho1):
                        if self.dic[r][theta][tau][rho1].has_key(rho2):
                            if self.dic[r][theta][tau][rho1][rho2].has_key(rho3):
                                self.dic[r][theta][tau][rho1][rho2][rho3] = val
                            else:
                                self.dic[r][theta][tau][rho1][rho2][rho3] = val
                        else:
                            rho3 = { tmprho3 : val }
                            self.dic[r][theta][tau][rho1][rho2] = rho3
                    else:
                        rho3 = { tmprho3 : val }
                        rho2 = { tmprho2 : rho3 }
                        self.dic[r][theta][tau][rho1] = rho2
                else:
                    rho3 = { tmprho3 : val }
                    rho2 = { tmprho2 : rho3 }
                    rho1 = { tmprho1 : rho2 }
                    self.dic[r][theta][tau] = rho1
            else:
                rho3 = { tmprho3 : val }
                rho2 = { tmprho2 : rho3 }
                rho1 = { tmprho1 : rho2 }
                tau = { tmptau : rho1 }
                self.dic[r][theta] = tau
        else:
            rho3 = { tmprho3 : val }
            rho2 = { tmprho2 : rho3 }
            rho1 = { tmprho1 : rho2 }
            tau = { tmptau : rho1 }
            theta = { tmptheta : tau }
            self.dic[r] = theta

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

        self.dic = Dic()

        self.r_oh = 0.97167
        self.theta_hoh = np.pi * 104.5/ 180.0

        self.varyR =     False
        self.varyTheta = False
        self.varyTau =   False
        self.varyRho1 =  False
        self.varyRho2 =  False
        self.varyRho3 =  False

        self.optionsR =     {  "min": 5.00, "max" : 5.00, "points": 1 }
        self.optionsTheta = {  "min": 0.00, "max" : 0.00, "points": 1 }
        self.optionsTau =   {  "min": 0.00, "max" : 0.00, "points": 1 }
        self.optionsRho1 =  {  "min": 0.00, "max" : 0.00, "points": 1 }
        self.optionsRho2 =  {  "min": 0.00, "max" : 0.00, "points": 1 }
        self.optionsRho3 =  {  "min": 0.00, "max" : 0.00, "points": 1 }

        opts = { "r" : {"min":2.30, "max":5.00,  "points":100} ,
             "theta" : {"max": np.pi , "min": 0.00 , "points":10},
             "tau"  : {"max": np.pi/2 , "min": 0.00 , "points":10},
             "rho1" : {"max": np.pi , "min": 0.00 , "points":1},
             "rho2" : {"max": np.pi , "min": 0.00 , "points":1},
             "rho3" : {"max": np.pi , "min": 0.00 , "points":1},
             }
        self.varyParameters()

    def genMols(self):
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
                                w1 = self.getWater( [0, 0, 0], self.r_oh, self.theta_hoh)
                                x, y, z = self.getCartesianFromDegree( i, j, k )
                                w2 = self.getWater( [x,y,z], self.r_oh, self.theta_hoh)
                                w2.rotate( l, m, n )
                                name = ""
                                name += "-".join( map( str, ["%3.2f"%i, "%3.2f"%j, "%3.2f"%k, "%3.2f"%l, "%3.2f"%m, "%3.2f"%n] ) )
                                name += ".mol"
                                self.writeMol( [w1, w2], name )

    def varyParameters( self, *args ):
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

    def getWater(self, origin, r, theta, AA = True ):
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
                tmp.addAtom( i )
                for j in atoms:
                    if j.element == "O":
                        continue
                    if j.inWater:
                        continue
#If in cartesian:
                    if j.AA:
                        if i.distToAtom(j) < 1.1:
                            tmp.addAtom ( j )
                            j.inWater = True
                    else:
                        if i.distToAtom(j) < 1.1/a0:
                            tmp.addAtom ( j )
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
#__Water__.addAtom() method will update the waters residue number and center coordinate
#When all atoms are there
#Right now NOT center-of-mass
                tmp.addAtom(i)
                for j in atoms:
                    if j.element != "H":
                        continue
                    if j.inWater:
                        continue
#1.05 because sometimes spc water lengths can be over 1.01
                        
                    if args.opAAorAU == "AA":
                        if i.dist(j) <= 1.05:
                            j.inWater = True
                            tmp.addAtom( j )
                            if len(tmp.atomlist) == 3:
                                break
                    elif args.opAAorAU == "AU":
                        if i.dist(j) <= 1.05/a0:
                            j.inWater = True
                            tmp.addAtom( j )
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
                tmp.addAtom( i )
                for j in atoms:
                    if j.element == "O":
                        continue
                    if j.inWater:
                        continue
#If in cartesian:
                    if i.AA:
                        if i.dist(j) < 1.0:
                            tmp.addAtom ( j )
                            j.inWater = True
                    else:
                        if i.dist(j) < 1.0/a0:
                            tmp.addAtom ( j )
                            j.inWater = True
                tmp.number = cnt
                cnt += 1
                waters.append( tmp )
        return waters
    def writeMol(self, wlist, name = "tmp.mol" ):
        f_ = open (name, 'w')
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

    def getCartesianFromDegree(self, r, theta, tau):
        return r* m.sin( m.pi*theta/180.0 )*m.cos( m.pi*tau/180.0) \
               , r* m.sin( m.pi*theta/180.0 )*m.sin( m.pi*tau/180.0)  \
               , r* m.cos( m.pi*theta/180.0 )





class Atom:
    def __init__(self, *args, **kwargs ):

        self.element = None
        self.polar = False
        self.cartesian = False
        self.theta = False
        self.tau = False
        self.r = False
        self.x = None
        self.y = None
        self.resid = None
        self.inWater = False

        self.AA = True

    def __str__(self):
        return "%s %f %f %f" %(self.element, self.x, self.y, self.z)
    def __sub__(self, other ):
        """return numpy array between this and other atom"""
        return self.getArray() - other.getArray()
    def getArray(self):
        return np.array( [self.x , self.y, self.z ] ).copy()
    def distToAtom(self, other):
        return np.sqrt( (self.x - other.x)**2 + (self.y -other.y)**2 + (self.z -other.z)**2 )
    def toAU(self):
        if self.AA:
            self.x /= a0
            self.y /= a0
            self.z /= a0
            self.AA = False
    def toAA(self):
        if not self.AA:
            self.x *= a0
            self.y *= a0
            self.z *= a0
            self.AA = True

class Water:
    def __init__(self , *args, **kwargs):

        self.atoms = 0
        self.q = 0.0
        self.r_oh = False
        self.theta_hoh = False
#Parameters describing water location and rotation in space
        self.r = False
        self.theta = False
        self.tau = False
        self.euler1 = False
        self.euler2 = False
        self.euler3 = False

        self.bNoHydrogens = True
        self.h1 = False
        self.h2 = False
        self.o  = False
#By default put center on oxygen, put center of mass if args.com
        self.center  = False
        self.resId = 0
        self.atomlist  = []

        self.AA = True
    def addAtom(self, atom):
        if self.atoms > 3:
            print "tried to add additional atoms to water, exiting"
            raise SystemExit
        self.atoms += 1
        if atom.element == "H":
            if self.bNoHydrogens:
                self.h1 = atom
                self.bNoHydrogens = False
            else:
                self.h2 = atom
        if atom.element == "O":
            self.o = atom
#Keep track of atomlist for easier looping
        self.atomlist.append(atom)
#Define water center, implement center-of mass later
        if (self.h1 and self.h2 and self.o):
            self.center = np.array([self.h1.x + self.h2.x + self.o.x,  \
                    self.h1.y + self.h2.y + self.o.y , \
                    self.h1.z + self.h2.z + self.o.z ]) / 3.0
            hm = mass_dict[ self.h1.element ]
            om = mass_dict[ self.h1.element ]
            self.com = np.array([ self.h1.x * hm  + self.h2.x *hm + self.o.x *om,  \
                self.h1.y *hm + self.h2.y *hm + self.o.y *om , \
                self.h1.z *hm + self.h2.z *hm + self.o.z *om ]) /( 2*hm +om)

        if self.resId:
            if self.resId != atom.resid:
                print "Tried to add %s to %s, exiting" %(atom, self)
                raise SystemExit
        else:
#Initialize water resId from atomic resid
            self.resId = atom.resid

    def __str__(self):
        return str("WAT") + str(self.resId) 

    def getAngleRho(self, other):
        d1 = self.getDipole()
        d2 = other.getDipole()
        return np.arccos( np.dot( d1, d2)/ ( np.linalg.norm(d1) * np.linalg.norm(d2) ) )

    def getAngleTau(self, other):
        r1= self.getNorm()
        r2= other.getNorm()
        return np.arccos( np.dot( r1, r2 ) / (np.linalg.norm( r1 ) * np.linalg.norm( r2 )))

    def getDipole(self):
        hq = 0.25
        oq = -0.5
        return self.h1.getArray() * hq + self.h2.getArray() * hq + self.o.getArray() * oq

    def getNorm(self):
        r1 = self.h1 - self.o
        r2 = self.h2 - self.o
        return np.cross( r1, r2 )

    def distToPoint( self , point ):
        return m.sqrt( (self.center[0] - point[0])**2 +\
                (self.center[1] - point[1])**2  + ( self.center[2] -point[2])**2 )

    def distToWater(self, other):
        xyz1 = self.center
        xyz2 = other.center
        return m.sqrt( (xyz1[0] - xyz2[0])**2 + \
            (xyz1[1] - xyz2[1])**2 + (xyz1[2] - xyz2[2])**2 )

    def getEuler(self):
        """Return euler angles required to rotate water in oxygen at origo to current"""

        H1 = self.h1.getArray()
        H2 = self.h2.getArray()
        O1 = self.o.getArray()

        dip = self.getDipole()

        origin = O1.copy()
        H1, H2, O1 = H1 - origin, H2 - origin, O1 - origin

        theta1 = m.atan2( dip[1], dip[0])
        H1 =  np.dot( self.getRzinv( theta1 ) , H1 )
        H2 =  np.dot( self.getRzinv( theta1 ) , H2 )
        O1 =  np.dot( self.getRzinv( theta1 ) , O1 )
        dip = np.dot( self.getRzinv( theta1 ) , dip )
#Rotate by theta around y axis so that the dipole is in the z axis 
        theta2 = m.atan2( -dip[0], dip[2])
        H1 =  np.dot( self.getRy( theta2 ) , H1 )
        H2 =  np.dot( self.getRy( theta2 ) , H2 )
        O1 =  np.dot( self.getRy( theta2 ) , O1 )
        dip = np.dot( self.getRy( theta2 ) , dip )
#Rotate around Z axis so that hydrogens are in xz plane.
        if H2[1] >0:
            xc = H2[0]
            yc = H2[1]
        else:
            xc = H1[0]
            yc = H1[1]
        theta3 = m.atan2( yc , xc)

        return theta3, theta2, theta1

    def rotate(self, t1, t2, t3):
        """Rotate self by t1, t2 and t3
        first Rz with theta1, then Ry^-1 by theta2, then Rz with theta 3"""
        d1, d2, d3 = self.getEuler()
        
# Place water molecule in origo, and rotate it so hydrogens in xz plane
        H1 = self.h1.getArray() ; H2 = self.h2.getArray() ; O = self.o.getArray()
        TMP = self.o.getArray()
        H1 -= TMP ; H2 -= TMP; O -= TMP

        H1 = np.dot( self.getRzinv(d3) , H1 )
        H1 = np.dot( self.getRy(d2) , H1 )
        H1 = np.dot( self.getRzinv(d1) , H1 )

        H2 = np.dot( self.getRzinv(d3) , H2 )
        H2 = np.dot( self.getRy(d2) , H2 )
        H2 = np.dot( self.getRzinv(d1) , H2 )

        O = np.dot( self.getRzinv(d3) , O )
        O = np.dot( self.getRy(d2) , O )
        O = np.dot( self.getRzinv(d1) , O )

# Rotate with angles t1, t2, t3

        H1 = np.dot( self.getRz(t1) , H1 )
        H1 = np.dot( self.getRyinv(t2) , H1 )
        H1 = np.dot( self.getRz(t3) , H1 )

        H2 = np.dot( self.getRz(t1) , H2 )
        H2 = np.dot( self.getRyinv(t2) , H2 )
        H2 = np.dot( self.getRz(t3) , H2 )

        O = np.dot( self.getRz(t1) , O )
        O = np.dot( self.getRyinv(t2) , O )
        O = np.dot( self.getRz(t3) , O )

#Put back in oxygen original point
        H1 += TMP ; H2 += TMP; O += TMP

        self.h1.x = H1[0] ;self.h1.y = H1[1] ;self.h1.z = H1[2] 
        self.h2.x = H2[0] ;self.h2.y = H2[1] ;self.h2.z = H2[2] 
        self.o.x  =  O[0] ;  self.o.y = O[1] ;  self.o.z = O[2] 
    def getRz( self, theta ):
        vec = np.array( [[ m.cos(theta), -m.sin(theta), 0],
                            [ m.sin(theta), m.cos(theta), 0],
                            [ 0,    0,  1]])
        return vec
    def getRzinv( self, theta ):
        vec = np.array(     [[ m.cos(theta), m.sin(theta), 0],
                            [ -m.sin(theta), m.cos(theta), 0],
                            [ 0,             0,            1]])
        return vec
    def getRy( self, theta ):
        vec = np.array( [[ m.cos(theta),0, m.sin(theta)],
                            [ 0,    1,  0],
                            [ -m.sin(theta), 0, m.cos(theta)]])
        return vec
    def getRyinv( self, theta ):
        vec = np.array( [[ m.cos(theta),0, -m.sin(theta)],
                            [ 0,    1,  0],
                            [ m.sin(theta), 0, m.cos(theta)]])
        return vec
    def plotWater(self ):
#Plot water molecule in green and  nice xyz axis
        O1, H1, H2 = self.o, self.h1, self.h2
        fig = plt.figure()
        dip = self.getDipole()
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

    def transformDipole( self, qmdipole, t1, t2, t3 ):

        d_new1 = np.zeros([3]) #will be returned
        d_new2 = np.zeros([3]) #will be returned
        d_new3 = np.zeros([3]) #will be returned

        rz  = self.getRz( t1 )
        ryi = self.getRyinv( t2 )
        rz2 = self.getRz( t3 )

        for i in range(3):
            for x in range(3):
                d_new1[i] += rz[i][x] * qmdipole[x]
        for i in range(3):
            for x in range(3):
                d_new2[i] += ryi[i][x] * d_new1[x]
        for i in range(3):
            for x in range(3):
                d_new3[i] += rz2[i][x] * d_new2[x]
        self.qmDipole = d_new3
        return d_new3

    def transformAlpha( self, qmalpha, t1, t2 , t3 ):

        a_new1 = np.zeros([3,3]) #will be calculated
        a_new2 = np.zeros([3,3]) #will be calculated
        a_new3 = np.zeros([3,3]) #will be calculated

        rz = self.getRz( t1 )
        ryi = self.getRyinv( t2 )
        rz2 = self.getRz( t3 )

        #print 'inside transfer alpha, water: %d at x = %f' %( self.resId, self.o.x)
        #print np.sqrt( (self.alpha[0][0] +self.alpha[0][1] +self.alpha[0][2] )**2 +\
        #               (self.alpha[1][0] +self.alpha[1][1] +self.alpha[1][2] )**2 +\
        #               (self.alpha[2][0] +self.alpha[2][1] +self.alpha[2][2] )**2 )

        for i in range(3):
            for j in range(3):
                for x in range(3):
                    for y in range(3):
                        a_new1[i][j] += rz[i][x] * rz[j][y] * qmalpha[x][y]
        for i in range(3):
            for j in range(3):
                for x in range(3):
                    for y in range(3):
                        a_new2[i][j] += ryi[i][x] * ryi[j][y] * a_new1[x][y]
        for i in range(3):
            for j in range(3):
                for x in range(3):
                    for y in range(3):
                        a_new3[i][j] += rz2[i][x] * rz2[j][y] * a_new2[x][y]
        return a_new3

    def transformBeta( self, qmbeta, t1, t2, t3 ):

        b_new1 = np.zeros([3,3,3]) #will be calculated
        b_new2 = np.zeros([3,3,3]) #will be calculated
        b_new3 = np.zeros([3,3,3]) #will be calculated

        rz =  self.getRz( t1 )
        ryi = self.getRyinv( t2 )
        rz2 = self.getRz( t3 )

        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for x in range(3):
                        for y in range(3):
                            for z in range(3):
                                b_new1[i][j][k] += rz[i][x] * rz[j][y] * rz[k][z] * qmbeta[x][y][z]
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for x in range(3):
                        for y in range(3):
                            for z in range(3):
                                b_new2[i][j][k] += ryi[i][x] * ryi[j][y] * ryi[k][z] * b_new1[x][y][z]
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for x in range(3):
                        for y in range(3):
                            for z in range(3):
                                b_new3[i][j][k] += rz2[i][x] * rz2[j][y] * rz2[k][z] * b_new2[x][y][z]

        return b_new3

    def beta_par(self):
        """
        return beta along molecular dipole z axis as defined by formula:
        b_par = (bxxx + bxyy + bxzz, byxx +byyy + byzz, bzxx + bzyy + bzzz ) x dipole

        """
        beta_par = np.zeros( [3] )
        for i in range(len(self.beta)):
            for j in range(len(self.beta)):
                for k in range(len(self.beta)):
                    if j == k:
                        beta_par[i] += self.beta[i][j][k]
        return np.dot(beta_par, self.dipole) / np.linalg.norm( self.dipole )

    def beta_par2(self):
        """
        return beta paralell as defined by formula:
        b_par = 1/5 * sum (bzii + bizi + biiz) %i = x, y, z

        """
        beta_par = 0.0
        for i in range(len(self.beta)):
            for j in range(len(self.beta)):
                for k in range(len(self.beta)):
                    if j == k:
                        if i == 2:
                            beta_par += self.beta[i][j][k]
                    if i == j:
                        if k == 2:
                            beta_par += self.beta[i][j][k]
                    if i == k:
                        if j == 2:
                            beta_par += self.beta[i][j][k]
        beta_par /= 5.0
        return beta_par

    def getSquareDipole(self):
        return np.sqrt( self.qmDipole[0] **2 + self.qmDipole[1]**2 + self.qmDipole[2]**2 )

    def getSquareBeta(self):
        return np.sqrt( \
                (self.qmBeta[0][0][0] + self.qmBeta[0][1][1] + self.qmBeta[0][2][2] )**2 + \
                (self.qmBeta[1][0][0] + self.qmBeta[1][1][1] + self.qmBeta[1][2][2] )**2 + \
                (self.qmBeta[2][0][0] + self.qmBeta[2][1][1] + self.qmBeta[2][2][2] )**2  )

    def getAlphaTrace(self):
        return  self.qmAlpha[0][0] + self.qmAlpha[1][1] + self.qmAlpha[2][2]

    def toAU(self):
        if self.AA:
            self.h1.toAU()
            self.h2.toAU()
            self.o.toAU()
            self.AA = False
    def toAA(self):
        if not self.AA:
            self.h1.toAA()
            self.h2.toAA()
            self.o.toAA()
            self.AA = True
            
if __name__ == '__main__':

    g = Generator()
    w1 = g.getWater( [0,0,0], 1.0, 104.5 )

    raise SystemExit
# Water bonding parameters:
    r_oh = 0.94 ; theta_hoh = 104.5
# Water rotation parameters
    r = 5.0
    #theta = 45 ; tau = 45
    euler1 = 0.0 ; euler2 = 0.0 ; euler3 = 0.0
    p = m.pi ; c = m.cos ; s = m.sin
    euler1 *= p/180 ;euler2 *= p/180 ;euler3 *= p/180 ;
    theta_hoh *= p/180

#Hardcoded conversions
    for i in np.r_[0 : 180: 8j ]:
        for j in np.r_[0 : 360 : 16j ]:

            theta = i * p/180 
            tau = j * p/180 
            x = r * s(theta) * c(tau )
            y = r * s(theta) * s(tau )
            z = r * c(theta)
            w1 = g.getWater( [ 0, 0, 0], r_oh, theta_hoh )
            w1.resId = 1 ; w1.r = r ; w1.theta = theta ; w1.tau = tau
            w2 = g.getWater( [ x, y, z], r_oh, theta_hoh )
            w2.resId = 2 ; w2.r = r ; w2.theta = theta ; w2.tau = tau
            g.writeMol( [ w1, w2 ])


