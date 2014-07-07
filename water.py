#!/usr/bin/env python
#-*- coding: utf-8 -*-

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt

import numpy as np
import math as m

from generator import *

a0 = 0.52917721092
charge_dic = {"H": 1.0, "C": 6.0, "N": 7.0, "O": 8.0, "S": 16.0}
mass_dict = {"H": 1.008,  "C": 6.0, "N": 7.0, "O": 15.999, "S": 16.0}

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
        self.z = None
        self.resid = None
        self.inWater = False
    def __str__(self):
        return "%s %f %f %f" %(self.element, self.x, self.y, self.z)

    def __sub__(self, other ):
        """return numpy array between this and other atom"""
        return self.getArray() - other.getArray()

    def getArray(self):
        return np.array( [self.x , self.y, self.z ] ).copy()
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

    def transfer_alpha( self, qmalpha, t1, t2 , t3 ):

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
        self.qmAlpha = a_new3

    def transferBeta( self, qmbeta, t1, t2, t3 ):

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

        self.qmBeta = b_new3

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


