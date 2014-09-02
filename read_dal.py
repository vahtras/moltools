#!/usr/bin/env python

import os,sys, re, argparse
import numpy as np
import fractions as fr
import math as m

from particles import *
from calculator import *


a0 = 0.52917721092
lab = [ "X", "Y", "Z"]
charge_dic = {"H": 1.0, "C": 6.0, "N": 7.0, "O": 8.0, "S": 16.0}
mass_dict = {"H": 1.008,  "C": 6.0, "N": 7.0, "O": 15.999, "S": 16.0}

class Water:
    def __init__(self):
        self.atoms = 0
        self.q = 0.0
        self.bNoHydrogens = True
        self.h1 = False
        self.h2 = False
        self.o  = False

#By default put center on oxygen, put center of mass if args.com
        self.center  = False

        self.number = 0
        self.atomlist  = []

    def addAtom(self, atom):
        if self.atoms > 3:
            print "tried to add additional atoms to water, exiting"
            raise SystemExit
        self.atoms += 1
        if atom.atype == "H":
            if self.bNoHydrogens:
                self.h1 = atom
                self.bNoHydrogens = False
            else:
                self.h2 = atom
        if atom.atype == "O":
            self.o = atom
#Keep track of atomlist for easier looping
        self.atomlist.append(atom)
#Define water center, implement center-of mass later
        if (self.h1 and self.h2 and self.o):
            self.center = np.array([self.h1.x + self.h2.x + self.o.x,  \
                    self.h1.y + self.h2.y + self.o.y , \
                    self.h1.z + self.h2.z + self.o.z ]) / 3.0
            hm = mass_dict[ self.h1.atype ]
            om = mass_dict[ self.h1.atype ]
            self.com = np.array([ self.h1.x * hm  + self.h2.x *hm + self.o.x *om,  \
                self.h1.y *hm + self.h2.y *hm + self.o.y *om , \
                self.h1.z *hm + self.h2.z *hm + self.o.z *om ]) /( 2*hm +om)

        if self.number:
            if self.number != atom.resid:
                print "Tried to add %s to %s, exiting" %(atom, self)
                raise SystemExit
        else:
#Initialize water number from atomic resid
            self.number = atom.resid

    def __str__(self):
        return str("WAT") + str(self.number) 

    def distToPoint( self , point ):
        return m.sqrt( (self.center[0] - point[0])**2 +\
                (self.center[1] - point[1])**2  + ( self.center[2] -point[2])**2 )

    def distToWater(self, other):
        xyz1 = self.center
        xyz2 = other.center
        return m.sqrt( (xyz1[0] - xyz2[0])**2 + \
            (xyz1[1] - xyz2[1])**2 + (xyz1[2] - xyz2[2])**2 )

    def get_euler(self):
        def Rzinv( theta ):
            vec = np.array(     [[ m.cos(theta), m.sin(theta), 0],
                                [ -m.sin(theta), m.cos(theta), 0],
                                [ 0,             0,            1]])
            return vec
        def Rz( theta ):
            vec = np.array(    [[ m.cos(theta), -m.sin(theta), 0],
                                [ m.sin(theta),  m.cos(theta), 0],
                                [ 0,         0,                1]])
            return vec
        def Ry( theta ):
            vec = np.array( [[ m.cos(theta),0, m.sin(theta)],
                                [ 0,    1,  0],
                                [ -m.sin(theta), 0, m.cos(theta)]])
            return vec
        def Ryinv( theta ):
            vec = np.array( [[ m.cos(theta),0, -m.sin(theta)],
                                [ 0,    1,  0],
                                [ m.sin(theta), 0, m.cos(theta)]])
            return vec

        oq = -0.5
        hq = 0.25

        H1 = np.array([ self.h1.x, self.h1.y, self.h1.z] )
        H2 = np.array([ self.h2.x, self.h2.y, self.h2.z] )
        O1 = np.array([ self.o.x,  self.o.y,  self.o.z ] )

        dip = oq * O1 + hq *H1 + hq * H2

        origin = O1.copy()
        H1, H2, O1 = H1 - origin, H2 - origin, O1 - origin

        theta1 = m.atan2( dip[1], dip[0])
        H1 =  np.dot( Rzinv( theta1 ) , H1 )
        H2 =  np.dot( Rzinv( theta1 ) , H2 )
        O1 =  np.dot( Rzinv( theta1 ) , O1 )
        dip = np.dot( Rzinv( theta1 ) , dip )
#Rotate by theta around y axis so that the dipole is in the z axis 
        theta2 = m.atan2( -dip[0], dip[2])
        H1 =  np.dot( Ry( theta2 ) , H1 )
        H2 =  np.dot( Ry( theta2 ) , H2 )
        O1 =  np.dot( Ry( theta2 ) , O1 )
        dip = np.dot( Ry( theta2 ) , dip )
#Rotate around Z axis so that hydrogens are in xz plane.
        if H2[1] >0:
            xc = H2[0]
            yc = H2[1]
        else:
            xc = H1[0]
            yc = H1[1]
        theta3 = m.atan2( yc , xc)
        H1 =  np.dot( Rzinv( theta3 ) , H1 )
        H2 =  np.dot( Rzinv( theta3 ) , H2 )
        O1 =  np.dot( Rzinv( theta3 ) , O1 )
        dip = np.dot( Rzinv( theta3 ) , dip )

        self.theta1 = theta1
        self.theta2 = theta2
        self.theta3 = theta3

    def transfer_dipole( self ):
        def Rz( theta ):
            vec = np.array( [[ m.cos(theta), -m.sin(theta), 0],
                                [ m.sin(theta), m.cos(theta), 0],
                                [ 0,    0,  1]])
            return vec
        def Ryinv( theta ):
            vec = np.array( [[ m.cos(theta),0, -m.sin(theta)],
                                [ 0,    1,  0],
                                [ m.sin(theta), 0, m.cos(theta)]])
            return vec
        d_new1 = np.zeros([3]) #will be returned
        d_new2 = np.zeros([3]) #will be returned
        d_new3 = np.zeros([3]) #will be returned

        rz = Rz( self.theta3 )
        ryi = Ryinv( self.theta2 )
        rz2 = Rz( self.theta1 )

        for i in range(3):
            for x in range(3):
                d_new1[i] += rz[i][x] * self.dipole[x]
        for i in range(3):
            for x in range(3):
                d_new2[i] += ryi[i][x] * d_new1[x]
        for i in range(3):
            for x in range(3):
                d_new3[i] += rz2[i][x] * d_new2[x]

        self.dipole = d_new3

    def transfer_alpha( self ):
        def Rz( theta ):
            vec = np.array( [[ m.cos(theta), -m.sin(theta), 0],
                                [ m.sin(theta), m.cos(theta), 0],
                                [ 0,    0,  1]])
            return vec
        def Ryinv( theta ):
            vec = np.array( [[ m.cos(theta),0, -m.sin(theta)],
                                [ 0,    1,  0],
                                [ m.sin(theta), 0, m.cos(theta)]])
            return vec

        a_new1 = np.zeros([3,3]) #will be calculated
        a_new2 = np.zeros([3,3]) #will be calculated
        a_new3 = np.zeros([3,3]) #will be calculated

        rz = Rz( self.theta3 )
        ryi = Ryinv( self.theta2 )
        rz2 = Rz( self.theta1 )

        #print 'inside transfer alpha, water: %d at x = %f' %( self.number, self.o.x)
        #print np.sqrt( (self.alpha[0][0] +self.alpha[0][1] +self.alpha[0][2] )**2 +\
        #               (self.alpha[1][0] +self.alpha[1][1] +self.alpha[1][2] )**2 +\
        #               (self.alpha[2][0] +self.alpha[2][1] +self.alpha[2][2] )**2 )

        for i in range(3):
            for j in range(3):
                for x in range(3):
                    for y in range(3):
                        a_new1[i][j] += rz[i][x] * rz[j][y] * self.alpha[x][y]
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
        self.alpha = a_new3

    def transfer_beta( self ):
        def Rz( theta ):
            vec = np.array( [[ m.cos(theta), -m.sin(theta), 0],
                                [ m.sin(theta), m.cos(theta), 0],
                                [ 0,    0,  1]])
            return vec
        def Ryinv( theta ):
            vec = np.array( [[ m.cos(theta),0, -m.sin(theta)],
                                [ 0,    1,  0],
                                [ m.sin(theta), 0, m.cos(theta)]])
            return vec
        b_new1 = np.zeros([3,3,3]) #will be calculated
        b_new2 = np.zeros([3,3,3]) #will be calculated
        b_new3 = np.zeros([3,3,3]) #will be calculated

        rz = Rz( self.theta3 )
        ryi = Ryinv( self.theta2 )
        rz2 = Rz( self.theta1 )

        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for x in range(3):
                        for y in range(3):
                            for z in range(3):
                                b_new1[i][j][k] += rz[i][x] * rz[j][y] * rz[k][z] * self.beta[x][y][z]
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

        self.beta = b_new3

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

    def alpha_par(self):
        """
        return alpha along molecular dipole z axis as defined by formula:
        a_par = (axx + axy + axz, ayx + ayy + ayz, azx + azy + azz ) x dipole

        """
        alpha_par = np.zeros( [3] )
        for i in range(3):
            for j in range(3):
                    alpha_par[i] += self.alpha[i][j]
        return np.dot( alpha_par, self.dipole) / np.linalg.norm( self.dipole )

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

    def square_dipole(self):
        return np.sqrt( self.dipole[0] **2 + self.dipole[1]**2 + self.dipole[2]**2 )


    def square_beta(self):
        return np.sqrt( \
                (self.beta[0][0][0] + self.beta[0][1][1] + self.beta[0][2][2] )**2 + \
                (self.beta[1][0][0] + self.beta[1][1][1] + self.beta[1][2][2] )**2 + \
                (self.beta[2][0][0] + self.beta[2][1][1] + self.beta[2][2][2] )**2  )

    def alpha_trace(self):
        return  self.alpha[0][0] + self.alpha[1][1] + self.alpha[2][2]
               


class Atom:
    def __init__(self, atype, x, y ,z, resid ):


        self.atype = atype
        self.x = x
        self.y = y
        self.z = z
        self.resid = resid
        self.inWater = False

    def __str__(self):
        return "%s %f %f %f" %(self.atype, self.x, self.y, self.z)

    def toAU(self):
        self.x /= a0
        self.y /= a0
        self.z /= a0

    def toAA(self):
        self.x *= a0
        self.y *= a0
        self.z *= a0

    def dist(self, other):
        return np.sqrt( (self.x - other.x)**2 + (self.y -other.y)**2 + (self.z -other.z)**2 )

class Templates:
    def __init__(self):

        monomer1_hf_cc_pCDZ  =  \
        [
#Dipole
        [ 0.334578 , -0.592384 , 0.430175 ],
#Alpha
        [[ 3.861254 , -1.271224 , -0.711365  ] ,
         [ -1.271224 , 5.014419 , 0.736731  ] ,
         [ -0.711365 , 0.736731 , 6.765428 ]] ,
#Beta
        [[[ -5.346496 , 5.451527 , 1.234515 ],
          [ 5.451527 , -6.627448 , -3.935797 ],
          [ 1.234515 , -3.935796 , -1.538262 ]],

         [[ 5.451527 , -6.627448 ,  -3.935797 ],
          [ -6.627448 , 13.610069 , 4.816182 ],
          [ -3.935796 , 4.816181 , 4.859196 ]],

         [[ 1.234515 , -3.935797 , -1.538263 ],
          [ -3.935797 , 4.816182 , 4.859198 ],
          [ -1.538263 , 4.859198 , -23.416047 ]]]
         ]

        monomer2_hf_cc_pCDZ  =  \
        [
#Dipole
        [ 0.199591 , -0.481684 , 0.616255 ],
#Alpha
        [[ 6.418170 , 1.279270 , 0.579489  ] ,
         [ 1.279270 , 4.552088 , -0.857882  ] ,
         [ 0.579489 , -0.857882 , 4.261423 ]] ,
#Beta
        [[[ -11.012012 , 5.950673 , -11.806065 ],
          [ 5.950674 , 5.599055 , -3.860620 ],
          [ -11.806066 , -3.860620 , -2.227571 ]],

         [[ 5.950673 , 5.599055 ,  -3.860620 ],
          [ 5.599055 , 10.272574 , -4.663896 ],
          [ -3.860620 , -4.663896 , 2.218040 ]],

         [[ -11.806065 , -3.860620 , -2.227570 ],
          [ -3.860620 , -4.663896 , 2.218040 ],
          [ -2.227570 , 2.218040 , -7.124011 ]]]
         ]




        olav_hf_cc_pVDZ =  \
        [
#Dipole
        [ 0.0, 0.0, 0.814458 ],
#Alpha
        [[ 7.204163 , 0.0 , 0.0  ] ,
         [ 0.0 , 3.034600 , 0.0  ] ,
         [ 0.0 , 0.0 , 5.223948 ]] ,
#Beta
        [[[ 0.0 , 0.0, -18.452810 ],
          [ 0.0 , 0.0, 0.0],
          [ -18.452810 , 0.0, 0.0 ]],

         [[ 0.0 , 0.0, 0.0 ],
          [ 0.0 , 0.0, -2.336562 ],
          [ 0.0 , -2.336562 , 0.0 ]],

         [[ -18.452813 , 0.0, 0.0 ],
          [ 0.0 , -2.336562 , 0.0 ],
          [ 0.0 , 0.0, -11.161749 ]]]
         ]

        centered_b3lyp_cc_pVDZ =  {}

#Template properties for CENTERED, r = 0.958019, theta = 104.5 water
        centered_hf_cc_pVDZ =  \
        [
#Dipole
        [ 0.0, 0.0, 0.809400 ],
#Alpha
        [[ 6.922537 , 0.0 , 0.0  ] ,
         [ 0.0 , 3.040036 , 0.0  ] ,
         [ 0.0 , 0.0 , 5.0931372 ]],
#Beta
        [[[ 0.0 , 0.0, -17.217324 ],
          [ 0.0 , 0.0, 0.0],
          [ -17.217324 , 0.0, 0.0 ]],

         [[ 0.0 , 0.0, 0.0 ],
          [ 0.0 , 0.0, -2.339154 ],
          [ 0.0 , -2.339154 , 0.0 ]],

         [[ -17.217324 , 0.0, 0.0 ],
          [ 0.0 , -2.339154 , 0.0 ],
          [ 0.0 , 0.0, -10.671809 ]]]
         ]

#  TIP3P model HF cc-pVDZ
        tip3p_hf_cc_pVDZ =  \
        [
#Dipole
        [ 0.0, 0.0, 0.808971 ],
#Alpha
        [[ 6.906544 , 0.0 , 0.0 ],
         [ 0.0 , 3.040337 , 0.0 ],
         [ 0.0 , 0.0 ,  5.084489 ]],
#Beta
        [[[ 0.0 , 0.0 , -17.144250 ],
          [ 0.0 , 0.0 , 0.0],
          [ -17.144250 , 0.0, 0.0]],

         [[ 0.0 , 0.0, 0.0],
          [ 0.0 , 0.0, -2.338925 ],
          [ 0.0 , -2.338925 , 0.0 ]],

         [[ -17.144250 , 0.0, 0.0],
          [ 0.0 , -2.338925 , 0.0],
          [ 0.0 , 0.0, -10.640297 ]]]
         ]



#Template properties for CENTERED, r = 0.958019, theta = 104.5 water
        centered_b3lyp_middle =  \
        [
#Dipole
        [ 0.0, 0.0, 0.731498 ],
#Alpha
        [[ 9.326228 , 0.0 , 0.0  ] ,
         [ 0.0 , 9.206878 , 0.0  ] ,
         [ 0.0 , 0.0 , 9.037904 ]],
#Beta
        [[[ 0.0 , 0.0, -9.717288 ],
          [ 0.0 , 0.0, 0.0],
          [ -9.717288 , 0.0, 0.0 ]],

         [[ 0.0 , 0.0, 0.0 ],
          [ 0.0 , 0.0, -3.984672 ],
          [ 0.0 , -3.984672 , 0.0 ]],

         [[ -9.717288 , 0.0, 0.0 ],
          [ 0.0 , -3.984672 , 0.0 ],
          [ 0.0 , 0.0, -4.716134 ]]]
         ]


#  TIP3P model
        tip3p_b3lyp_middle =  \
        [
#Dipole
        [ 0.0, 0.0, 0.731575 ],
#Alpha
        [[ 9.333366 , 0.0 , 0.0 ],
         [ 0.0 , 9.208982 , 0.0 ],
         [ 0.0 , 0.0 ,  9.042354 ]],
#Beta
        [[[ 0.0 , 0.0 , -9.744865 ],
          [ 0.0 , 0.0 , 0.0],
          [ -9.744865 , 0.0, 0.0]],

         [[ 0.0 , 0.0, 0.0],
          [ 0.0 , 0.0, -3.975595 ],
          [ 0.0 , -3.975595 , 0.0 ]],

         [[ -9.744865 , 0.0, 0.0],
          [ 0.0 , -3.975595 , 0.0],
          [ 0.0 , 0.0, -4.716610 ]]]
         ]

# SPC water model
# HF
        spc_hf_cc_pVDZ =  \
        [
#Dipole
        [ 0.0, 0.0, 0.792907 ],
#Alpha
        [[ 7.985773 , 0.0 , 0.0  ] ,
         [ 0.0 , 3.019117 , 0.0  ] ,
         [ 0.0 , 0.0 , 5.245267 ]],
#Beta
        [[[ 0.0 , 0.0, -20.972694 ],
          [ 0.0 , 0.0, 0.0],
          [ -20.972694 , 0.0, 0.0 ]],

         [[ 0.0 , 0.0, 0.0 ],
          [ 0.0 , 0.0, -2.232938 ],
          [ 0.0 , -2.232938 , 0.0 ]],

         [[ -20.972694 , 0.0, 0.0],
          [ 0.0 , -2.232938 , 0.0],
          [ 0.0 , 0.0, -11.478545 ]]]
         ]
        spc_b3lyp_middle =  \
        [
#Dipole
        [ 0.0, 0.0, 0.707846 ],
#Alpha
        [[ 10.301482 , 0.0 , 0.0  ] ,
         [ 0.0 , 9.462288 , 0.0  ] ,
         [ 0.0 , 0.0 , 9.493345 ]],
#Beta
        [[[ 0.0 , 0.0, -13.567895 ],
          [ 0.0 , 0.0, 0.0],
          [ -13.567895 , 0.0, 0.0 ]],

         [[ 0.0 , 0.0, 0.0 ],
          [ 0.0 , 0.0, -2.959169 ],
          [ 0.0 , -2.959170 , 0.0 ]],

         [[ -13.567895 , 0.0, 0.0],
          [ 0.0 , -2.959169 , 0.0],
          [ 0.0 , 0.0, -4.996675 ]]]
         ]

    #Dictionaries from basis
        olav_hf_dict = { "PVDZ" : olav_hf_cc_pVDZ }
        centered_b3lyp_dict = { "MIDDLE" : centered_b3lyp_middle,
                "PVDZ" : centered_b3lyp_cc_pVDZ  }
        centered_hf_dict = { "PVDZ" : centered_hf_cc_pVDZ }
        tip3p_b3lyp_dict = { "MIDDLE" : tip3p_b3lyp_middle }
        tip3p_hf_dict = { "PVDZ" : tip3p_hf_cc_pVDZ }
        spc_b3lyp_dict = { "MIDDLE" : spc_b3lyp_middle }
        spc_hf_dict = { "PVDZ" : spc_hf_cc_pVDZ  }

        monomer1_hf_dict = { "PVDZ" : monomer1_hf_cc_pCDZ }
        monomer2_hf_dict = { "PVDZ" : monomer2_hf_cc_pCDZ }

    #Dictionaries from method
        olav_method_dict = { "HF" : olav_hf_dict }
        centered_method_dict = { "B3LYP" : centered_b3lyp_dict,
                "HF" : centered_hf_dict }
        tip3p_method_dict = { "B3LYP" : tip3p_b3lyp_dict,
                "HF" : tip3p_hf_dict }
        spc_method_dict = { "B3LYP": spc_b3lyp_dict, 
                "HF" : spc_hf_dict }
        monomer1_method_dict = { "HF" : monomer1_hf_dict }
        monomer2_method_dict = { "HF" : monomer2_hf_dict }

        self.nameDict = { "OLAV" : olav_method_dict ,
                "MON1" : monomer1_method_dict,
                "MON2" : monomer2_method_dict,
                "CENTERED" : centered_method_dict, \
                "TIP3P": tip3p_method_dict, \
                "SPC" : spc_method_dict }

    def getBeta(self, model, method, basis):

        return self.nameDict[model][method][basis]





def hyperq(vec):
    tmp = []
    vec_new = []
    for i in vec:
        if i[1] not in tmp:
            tmp.append(i[1])
            vec_new.append( i ) 
    return vec_new

def read_waters( args ):

#If the file is plain xyz file

    atoms = []
    if args.x.endswith( ".xyz" ) or args.x.endswith(".mol"):
        pat_xyz = re.compile(r'^\s*(\w+)\s+(-*\d*.+\d+)\s+(-*\d*.+\d+)\s+(-*\d*.+\d+) *$')
        for i in open( args.x ).readlines():
            if pat_xyz.match(i):
                f = pat_xyz.match(i).groups()
                tmpAtom = Atom(f[0][0], float(f[1]), float(f[2]), float(f[3]), 0)

                if args.xAAorAU == "AU":
                    if args.opAAorAU == "AA":
                        tmpAtom.toAA()

                elif args.xAAorAU == "AA":
                    if args.opAAorAU == "AU":
                        tmpAtom.toAU()

                atoms.append( tmpAtom )

    elif args.x.endswith( ".pdb" ):
        pat1 = re.compile(r'^(ATOM|HETATM)')
        for i in open( args.x ).readlines():
            if pat1.search(i):
                #Ignore charge centers for polarizable water models
                if ( i[11:16].strip() == "SW") or (i[11:16] == "DW"):
                    continue
                tmpAtom = Atom(i[11:16].strip()[0], \
                        float(i[30:38].strip()), \
                        float(i[38:46].strip()), \
                        float(i[46:54].strip()), \
                        int(i[22:26].strip()) )

                if args.xAAorAU == "AU":
                    if args.opAAorAU == "AA":
                        tmpAtom.toAA()
                elif args.xAAorAU == "AA":
                    if args.opAAorAU == "AU":
                        tmpAtom.toAU()
                atoms.append( tmpAtom )
    elif args.x.endswith( ".out" ):
        pat_xyz = re.compile(r'^(\w+)\s+(-*\d*.+\d+)\s+(-*\d*.+\d+)\s+(-*\d*.+\d+) *$')
        for i in open( args.x ).readlines():
            if pat_xyz.match(i):
                f = pat_xyz.match(i).groups()
                tmpAtom = Atom(f[0][0], float(f[1]), float(f[2]), float(f[3]), 0)
                if args.xAAorAU == "AU":
                    if args.opAAorAU == "AA":
                        tmpAtom.toAA()
                elif args.xAAorAU == "AA":
                    if args.opAAorAU == "AU":
                        tmpAtom.toAU()
                atoms.append( tmpAtom )
#loop over oxygen and hydrogen and if they are closer than 1 A add them to a water
    waters = []
    cnt = 1

    if args.x.endswith( ".xyz" ) or args.x.endswith(".mol"):
        for i in atoms:
            if i.atype == "H":
                continue
            if i.inWater:
                continue
            tmp = Water(  )
            i.inWater = True
            tmp.addAtom( i )
            for j in atoms:
                if j.atype == "O":
                    continue
                if j.inWater:
                    continue
#If in cartesian:
                if args.opAAorAU == "AA":
                    if i.dist(j) < 1.1:
                        tmp.addAtom ( j )
                        j.inWater = True
                elif args.opAAorAU == "AU":
                    if i.dist(j) < 1.1/a0:
                        tmp.addAtom ( j )
                        j.inWater = True
            tmp.number = cnt
            cnt += 1
            waters.append( tmp )
    elif args.x.endswith( ".pdb" ):
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
            if i.atype != "O":
                continue
            tmp = Water()
            i.inWater= True
#__Water__.addAtom() method will update the waters residue number and center coordinate
#When all atoms are there
#Right now NOT center-of-mass
            tmp.addAtom(i)
            for j in atoms:
                if j.atype != "H":
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
    elif args.x.endswith( ".out" ):
        for i in atoms:
            if i.atype == "H":
                continue
            if i.inWater:
                continue
            tmp = Water(  )
            i.inWater = True
            tmp.addAtom( i )
            for j in atoms:
                if j.atype == "O":
                    continue
                if j.inWater:
                    continue
#If in cartesian:
                if args.opAAorAU == "AA":
                    if i.dist(j) < 1.0:
                        tmp.addAtom ( j )
                        j.inWater = True
                elif args.opAAorAU == "AU":
                    if i.dist(j) < 1.0/a0:
                        tmp.addAtom ( j )
                        j.inWater = True
            tmp.number = cnt
            cnt += 1
            waters.append( tmp )
    return waters

def arg_parser( args ):
    A = argparse.ArgumentParser( description = \
            "This program reads alpha and beta from dalton .out files, obtained for an ideal water molecule centered as oxygen at origo and hydrogens in the xz-plane.\n\n It also read coordinates of arbitrary water molecules and transforms the above read properties to their coordinate reference frame, and writes it to a .pot file." ,add_help= True)
    A.add_argument( "-a", type = str, help="File that contains LINEAR response output with polarizabilities" )
    A.add_argument( "-al", type = str, default = "22", help="Symmetry number alpha [1, 2, 3]" )
    A.add_argument( "-b", type = str,help="File that contains QUADRATIC response output with hyperpolarizabilities" ) 
    A.add_argument( "-bl", type = str, default = "1", help="Symmetry number beta [1, 2, 3]" )
    A.add_argument( "-x", type = str, help = 'Coordinate file with water molecules for the output .pot file. [ xyz , pdb ]')
    A.add_argument( "-xAAorAU", type = str , default = "AA" ,
            help = 'Default coordinate type in AA or AU in -x input water coordinate file, default: "AA"')
    A.add_argument( "-dl", type = str, default = "1", help="Symmetry number dipole [0, 1]" )
    A.add_argument( "-test", action = "store_true", default = False )
    A.add_argument( "-mon", action = "store_true", default = False )
    A.add_argument("-v","--verbose", action='store_true' , default = False)
    A.add_argument("-write", nargs='*', default = [],  help = "Supply any which files to write from a selection: pot, xyz" )
    A.add_argument("-waters", type = int , default = 4, help = "how many waters to take closest to center atom, default: 4")
    A.add_argument( "-op", type = str, default = "conf.pot", help='output name of the .pot file, default: "conf.pot"' )
    A.add_argument( "-opAAorAU", type = str, default = "AU" , help='Default coordinate type AA or AU for -op output potential file, default: "AU"' )

    A.add_argument( "-ox", type = str, default = "conf.xyz", help='output name of the .xyz file, default: "conf.xyz"' )


    A.add_argument( "-com", help="Place point properties on center-of-mass instead of Oxygen",action = 'store_true', default = False)
    A.add_argument( "-template", action = 'store_true', help= "Activate Beta tensor reading from templates, options provided below for choices", default = False)

    A.add_argument( "-tname", type = str, default = "CENTERED", 
            help = "available templates: CENTERED (default), TIP3P, SPC, OLAV")

    A.add_argument( "-tmethod", type = str, default = "HF",
            help = "available methods: HF (default) , B3LYP")

    A.add_argument( "-tbasis", type = str, default = "PVDZ",
            help = "available choices: PVDZ (default, is actually cc-pVDZ), \
                    MIDDLE (ano type) [O: 5s3p2d, H: 3s1p]")


    a = A.parse_args( args[1:] )
    return a

def read_coords_and_dipole( args, **kwargs ):

    if args.a:
        file_ = args.a
    if args.b:
        file_ = args.b

    if "custom_file" in kwargs:
        file_ = kwargs[ "custom_file" ]

    nuc_dip = np.zeros(3)
    el_dip = np.zeros(3)

    pat_xyz = re.compile(r'^\s*(\w+)\s+(-*\d*\.+\d+)\s+(-*\d*\.+\d+)\s+(-*\d*\.+\d+) *$')
    pat_pol = re.compile(r'([XYZ])DIPLEN.*total.*:')
    atoms = []

    for i in open( file_ ).readlines():
        if pat_xyz.match(i):
            f = pat_xyz.match(i).groups()
            tmpAtom = Atom(f[0][0], float(f[1]), float(f[2]), float(f[3]), 1 )
            if args.xAAorAU == "AA":
                tmpAtom.toAU()
            atoms.append( tmpAtom )
        if pat_pol.search(i):
            if pat_pol.search(i).group(1) == "X":
                try:
                    if "D" in i.split()[3]:
                        frac = float(i.split()[3].replace("D","E"))
                    else:
                        frac = float(i.split()[3])
                except IndexError:
                    if "D" in i.split()[2]:
                        frac = float( i.split()[2].strip(":").replace("D","E"))
                    else:
                        frac = float( i.split()[2].strip(":"))
                el_dip[0] += frac
            if pat_pol.search(i).group(1) == "Y":
                try:
                    if "D" in i.split()[3]:
                        frac = float(i.split()[3].replace("D","E"))
                    else:
                        frac = float(i.split()[3])
                except IndexError:
                    if "D" in i.split()[2]:
                        frac = float( i.split()[2].strip(":").replace("D","E"))
                    else:
                        frac = float( i.split()[2].strip(":"))
                el_dip[1] += frac
            if pat_pol.search(i).group(1) == "Z":
                try:
                    if "D" in i.split()[3]:
                        frac = float(i.split()[3].replace("D","E"))
                    else:
                        frac = float(i.split()[3])
                except IndexError:
                    if "D" in i.split()[2]:
                        frac = float( i.split()[2].strip(":").replace("D","E"))
                    else:
                        frac = float( i.split()[2].strip(":"))
                el_dip[2] += frac

    for i in atoms:
        nuc_dip[0] += charge_dic[ i.atype ] * i.x
        nuc_dip[1] += charge_dic[ i.atype ] * i.y
        nuc_dip[2] += charge_dic[ i.atype ] * i.z
    return atoms, nuc_dip - el_dip

def read_alpha( args ):
# Reading in Alpha tensor
    pat_alpha = re.compile(r'@ -<< ([XYZ])DIPLEN.*([XYZ])DIPLEN')
    alpha = np.zeros([3,3])
    lab = ["X", "Y", "Z"]
    for i in open( args.a ).readlines():
        if pat_alpha.match( i ):
            try:
                if "D" in i.split()[-1]:
                    frac = float( i.split()[-1].replace("D","E") )
                else:
                    frac = float( i.split()[-1] )
            except IndexError:
                if "D" in i.split()[-1]:
                    frac = float( i.split()[-1].strip("=").replace("D","E") )
                else:
                    frac = float( i.split()[-1].strip("=") )
            A = pat_alpha.match(i).groups(1)[0]
            B = pat_alpha.match(i).groups(1)[1]
            alpha[ lab.index( A ) , lab.index( B ) ]  = frac
            if A == "X" and B == "Y":
                alpha[ lab.index( B ) , lab.index( A ) ]  = frac
            if A == "X" and B == "Z":
                alpha[ lab.index( B ) , lab.index( A ) ]  = frac
            if A == "Y" and B == "Z":
                alpha[ lab.index( B ) , lab.index( A ) ]  = frac
    return alpha 

def read_beta( args ):

    nuc_dip = np.zeros(3)
    el_dip = np.zeros(3)
    alpha = np.zeros([3,3])
    beta = np.zeros([3,3,3])
    tmp = []
    atoms = []
    missing = {}
    exists = {}
    lab = ["X", "Y", "Z"]

    pat_xyz = re.compile(r'^\s*(\w+)\s+(-*\d*\.+\d+)\s+(-*\d*\.+\d+)\s+(-*\d*\.+\d+) *$')
    pat_pol = re.compile(r'([XYZ])DIPLEN.*total.*:')

# Reading in dipole
    for i in open( args.b ).readlines():
        if pat_xyz.match(i):
            f = pat_xyz.match(i).groups()
            tmpAtom = Atom(f[0][0], float(f[1]), float(f[2]), float(f[3]), 1 )
            if args.xAAorAU == "AA":
                tmpAtom.toAU()
            atoms.append( tmpAtom )
        if pat_pol.search(i):
            if pat_pol.search(i).group(1) == "X":
                try:
                    if "D" in i.split()[3]:
                        frac = float(i.split()[3].replace("D","E"))
                    else:
                        frac = float(i.split()[3])
                except IndexError:
                    if "D" in i.split()[2]:
                        frac = float( i.split()[2].strip(":").replace("D","E"))
                    else:
                        frac = float( i.split()[2].strip(":"))
                el_dip[0] += frac
            if pat_pol.search(i).group(1) == "Y":
                try:
                    if "D" in i.split()[3]:
                        frac = float(i.split()[3].replace("D","E"))
                    else:
                        frac = float(i.split()[3])
                except IndexError:
                    if "D" in i.split()[2]:
                        frac = float( i.split()[2].strip(":").replace("D","E"))
                    else:
                        frac = float( i.split()[2].strip(":"))
                el_dip[1] += frac
            if pat_pol.search(i).group(1) == "Z":
                try:
                    if "D" in i.split()[3]:
                        frac = float(i.split()[3].replace("D","E"))
                    else:
                        frac = float(i.split()[3])
                except IndexError:
                    if "D" in i.split()[2]:
                        frac = float( i.split()[2].strip(":").replace("D","E"))
                    else:
                        frac = float( i.split()[2].strip(":"))
                el_dip[2] += frac

    for i in atoms:
        nuc_dip[0] += charge_dic[ i.atype ] * i.x
        nuc_dip[1] += charge_dic[ i.atype ] * i.y
        nuc_dip[2] += charge_dic[ i.atype ] * i.z





# Reading in Alfa and Beta tensor
    pat_alpha = re.compile(r'@ QRLRVE:.*([XYZ])DIPLEN.*([XYZ])DIPLEN')
    for i in open( args.b ).readlines():
        if pat_alpha.match( i ):
            try:
                if "D" in i.split()[-1]:
                    frac = float( i.split()[-1].replace("D","E") )
                else:
                    frac = float( i.split()[-1] )
            except IndexError:
                if "D" in i.split()[-1]:
                    frac = float( i.split()[-1].strip("=").replace("D","E") )
                else:
                    frac = float( i.split()[-1].strip("=") )
            A = pat_alpha.match(i).groups(1)[0]
            B = pat_alpha.match(i).groups(1)[1]
            alpha[ lab.index( A ) , lab.index( B ) ]  = frac
            if A == "X" and B == "Y":
                alpha[ lab.index( B ) , lab.index( A ) ]  = frac
            if A == "X" and B == "Z":
                alpha[ lab.index( B ) , lab.index( A ) ]  = frac
            if A == "Y" and B == "Z":
                alpha[ lab.index( B ) , lab.index( A ) ]  = frac
    pat_beta = re.compile(r'@ B-freq')
    for i in open( args.b ).readlines():
        if pat_beta.match(i):
            try:
                if i.split()[7].lstrip("beta") in exists:
                    continue
                exists[ i.split()[7].lstrip("beta") ] = float(i.split()[9] )
            except ValueError:
                a, b, c = i.split()[9].lstrip("beta").strip("()").split(",")
                if i.split()[7].lstrip("beta") in missing:
                    continue
                missing[ i.split()[7].lstrip("beta") ] =  "(%s;%s,%s)"%(a,b,c)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                try:
                    beta[i][j][k] = exists[ "(%s;%s,%s)" %(lab[i],lab[j],lab[k])]
                except KeyError:
                    beta[i][j][k] = exists[ missing["(%s;%s,%s)"%(lab[i],lab[j],lab[k]) ] ]

    return atoms, nuc_dip - el_dip, alpha , beta


def main():
    """ 
    Program reads alpha or beta tensor and dipole moment from DALTON output
    """
    args = arg_parser( sys.argv )

    f_alpha = False
    f_beta = False
    f_waters = False

    a_l = args.al
    b_l = args.bl
    d_l = args.dl


#These are used to create string for olavs dipole list class using templates
    dipole = np.zeros( [3] )
    alpha = np.zeros( [3, 3] )
    beta = np.zeros( [3, 3, 3])


#To be read from -b hfqua_file.out

    dipole_qm = np.zeros( [3] )
    alpha_qm = np.zeros( [3, 3] )
    beta_qm = np.zeros( [3, 3, 3])

    waters = np.zeros( [] )

    if args.a:
        f_alpha = args.a

    if args.b:
        f_beta = args.b

    if args.x:
        f_waters = args.x

    pat_bReadDipole = re.compile (r'\.PROPAV')
    pat_bReadAlpha = re.compile (r'\*LINEAR')
    pat_bReadBeta = re.compile (r'\*QUADRATIC')

    bReadDipole = False
    bReadAlpha = False
    bReadBeta = False

#Check what to read if linear calculation is supplied
    if f_alpha:
        bReadDipole = True
        bReadAlpha = True
        for i in open( f_alpha ).readlines():
            if pat_bReadDipole.match(i):
                bReadDipole = True
            if pat_bReadAlpha.match(i):
                bReadAlpha = True
    #if f_alpha:
    #    for i in open( f_alpha ).readlines():
    #        if pat_bReadDipole.match(i):
    #            bReadDipole = True
    #        if pat_bReadAlpha.match(i):
    #            bReadAlpha = True

#Check what to read if quadratic calculation is supplied
    if f_beta:
        bReadDipole = True
        bReadBeta = True
        for i in open( f_beta ).readlines():
            if pat_bReadDipole.match(i):
                bReadDipole = True
            if pat_bReadBeta.match(i):
                bReadBeta = True
    #if f_beta:
    #    for i in open( f_beta ).readlines():
    #        if pat_bReadDipole.match(i):
    #            bReadDipole = True
    #        if pat_bReadBeta.match(i):
    #            bReadBeta = True

#Read dipole, polarizability and hyperpolarizability from out files supplied
    if f_alpha:
        if bReadDipole:
            atoms , dipole_qm = read_coords_and_dipole( args )
        if bReadAlpha:
            alpha_qm = read_alpha( args )

    if f_beta:
        if bReadDipole:
            atoms , dipole_qm = read_coords_and_dipole( args )
        if bReadBeta:
            atoms, dipole_qm , alpha_qm , beta_qm = read_beta( args )


#If no output files were given, use calculated properties at centered water with
#r_oh = 0.958019  theta = 104.5 ,  B3LYP/ 5s3p2d (O), 3s1p (H)

#choose -name [TIP3P, SPC, CENTERED] for which template properties to use, right now only these at B3LYP / 5s3p2d (O), 3s1p (H)

    if args.template:
        tmp = Templates().getBeta( args.tname, args.tmethod, args.tbasis )
        dipole = np.array( tmp[0] )
        alpha = np.array(  tmp[1] )
        beta = np.array(   tmp[2] )

#Read coordinates for water molecules where to put properties to
    if f_waters:
        waters = read_waters( args )

# Read in rotation angles for each water molecule follow by transfer of dipole, alpha and beta to coordinates
    if f_waters:
        for i in waters:
            i.dipole = dipole
            i.alpha = alpha
            i.beta = beta
            i.get_euler()
            i.transfer_dipole()
            i.transfer_alpha()
            i.transfer_beta()


# Read in rotation angles for each water molecule follow by transfer of dipole, alpha and beta to coordinates
    if f_waters and args.mon:
        for i in waters:
            tmp1 = Templates().getBeta( "MON1" , args.tmethod, args.tbasis )
            dipole1 = np.array( tmp1[0] )
            alpha1 = np.array(  tmp1[1] )
            beta1 = np.array(   tmp1[2] )

            tmp2 = Templates().getBeta( "MON2" , args.tmethod, args.tbasis )
            dipole2 = np.array( tmp2[0] )
            alpha2 = np.array(  tmp2[1] )
            beta2 = np.array(   tmp2[2] )

            i.o.toAA()
            if i.o.x < 0.1:
                dipole = dipole2
                alpha = alpha2
                beta = beta2
            else:
                dipole = dipole1
                alpha = alpha1
                beta = beta1
            i.o.toAU()

            i.dipole = dipole
            i.alpha = alpha
            i.beta = beta
            i.get_euler()
            #i.transfer_dipole()
            #i.transfer_alpha()
            #i.transfer_beta()

#Write to file, potential

#Only write output file if input alpha, beta, and coordinates are given
    if ("pot" in args.write) or (f_waters and f_alpha and f_beta):
        write_pot = True
    else:
        write_pot = False

    if "xyz" in args.write:
        write_xyz = True
    else:
        write_xyz = False

    if "mol" in args.write:
        write_mol = True
    else:
        write_mol = False

#Explicit printing to stdout for testing, only the model water from linear / quadratic calc is printed

    if args.verbose:
        for i in atoms:
            print i

    if args.verbose and bReadDipole:
        for i in range(len(dipole)):
            print "Dipole_%s: %f" %( lab[i], dipole_qm[i] )
    if args.verbose and bReadAlpha:
        for i in range(len(alpha)):
            for j in range(len(alpha[i])):
                print "Alpha_%s%s: %f" %( lab[i], lab[j] , alpha_qm[i][j] )

    if args.verbose and bReadBeta:
        for i in range(3):
            for j in range(3):
                print "Alpha_%s%s: %f" %( lab[i], lab[j] , alpha_qm[i][j] )
        for i in range(len(beta)):
            for j in range(len(beta[i])):
                for k in range(len(beta[j])):
                    print "Beta_%s%s%s: %f" %( lab[i], lab[j] , lab[k] , beta_qm[i][j][k] )

#Perform some tests
    if args.test:
        if f_waters:
            print "Performing tests of transfers"
            for i in waters:

                print "Water number %d \nSquare dipole:" % i.number
                print i.square_dipole()

                print "alpha trace: "
                print i.alpha_trace()

                print "Square beta:"
                print i.square_beta()

                print "alpha projected on dipole:"
                print i.alpha_par()


#Inline temporary code to generate static, polarizable and hyperpolarizable string
    

    if f_waters:
#put tensors as centered in center of mass instead of on oxygen
        if args.com:
            for i in waters:
                i.center = i.com
        else:
            for i in waters:
                i.center = [ i.o.x, i.o.y, i.o.z ]
        string_static = "AU\n%d 1 0\n" %len(waters)
        string_polarizable = "AU\n%d 1 2 1\n" %len(waters)
        string_hyperpolarizable = "AU\n%d 1 22 1\n" %len(waters)
        for i in waters:
            string_static += "%d %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n" %((i.number, i.center[0], 
                i.center[1], i.center[2], \
                           i.q, i.dipole[0], i.dipole[1], i.dipole[2], \
                       i.alpha[0][0], i.alpha[0][1], i.alpha[0][2], \
                                      i.alpha[1][1], i.alpha[1][2], \
                                                     i.alpha[2][2], \
                       i.beta[0][0][0], i.beta[0][0][1], i.beta[0][0][2], \
                                        i.beta[0][1][1], i.beta[0][1][2], \
                                                         i.beta[0][2][2], \
                                        i.beta[1][1][1], i.beta[1][1][2], \
                                                         i.beta[1][2][2], \
                                                         i.beta[2][2][2] )) 
            string_polarizable += "%d %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n" %((i.number, i.center[0],
                i.center[1], i.center[2] , \
                           i.q, i.dipole[0], i.dipole[1], i.dipole[2], \
                       i.alpha[0][0], i.alpha[0][1], i.alpha[0][2], \
                                      i.alpha[1][1], i.alpha[1][2], \
                                                     i.alpha[2][2], \
                       i.beta[0][0][0], i.beta[0][0][1], i.beta[0][0][2], \
                                        i.beta[0][1][1], i.beta[0][1][2], \
                                                         i.beta[0][2][2], \
                                        i.beta[1][1][1], i.beta[1][1][2], \
                                                         i.beta[1][2][2], \
                                                         i.beta[2][2][2] )) 
            string_hyperpolarizable += "%d %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n" %((i.number, 
                i.center[0], i.center[1], i.center[2], \
                           i.q, i.dipole[0], i.dipole[1], i.dipole[2], \
                       i.alpha[0][0], i.alpha[0][1], i.alpha[0][2], \
                                      i.alpha[1][1], i.alpha[1][2], \
                                                     i.alpha[2][2], \
                       i.beta[0][0][0], i.beta[0][0][1], i.beta[0][0][2], \
                                        i.beta[0][1][1], i.beta[0][1][2], \
                                                         i.beta[0][2][2], \
                                        i.beta[1][1][1], i.beta[1][1][2], \
                                                         i.beta[1][2][2], \
                                                         i.beta[2][2][2] )) 


# Here start writing the .pot file

    if write_pot:
        f_ = open( args.op , "w" )
        if f_waters:
# Write different output depending on d_l, a_l, b_l
# First way is default, which is upper triangular for alpha and beta
            if d_l == "1" and a_l == "22" and b_l == "1":
                f_.write( "%s\n" %args.opAAorAU + str(len(waters)) + " %s %s %s\n" \
                        %( d_l, a_l, b_l) )
                for i in waters:
                    f_.write( "%d %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n" %((i.number, i.o.x, i.o.y, i.o.z, \
                           i.q, i.dipole[0], i.dipole[1], i.dipole[2], \
                       i.alpha[0][0], i.alpha[0][1], i.alpha[0][2], \
                                      i.alpha[1][1], i.alpha[1][2], \
                                                     i.alpha[2][2], \
                       i.beta[0][0][0], i.beta[0][0][1], i.beta[0][0][2], \
                                        i.beta[0][1][1], i.beta[0][1][2], \
                                                         i.beta[0][2][2], \
                                        i.beta[1][1][1], i.beta[1][1][2], \
                                                         i.beta[1][2][2], \
                                                         i.beta[2][2][2] )) )
# If dipole, full alpha, and full tensor beta, (1, 3, 3) OLD ROUTINE 
            if d_l == "1" and a_l == "3" and b_l == "3":
                f_.write( "AA\n" + str(len(waters)) + " %s %s %s\n" \
                        %( d_l, a_l, b_l) )
                for i in waters:
                    f_.write( "%d %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n" %((i.number, i.o.x, i.o.y, i.o.z, \

                       i.q, i.dipole[0], i.dipole[1], i.dipole[2], \

                       i.alpha[0][0], i.alpha[0][1], i.alpha[0][2], \
                       i.alpha[1][0], i.alpha[1][1], i.alpha[1][2], \
                       i.alpha[2][0], i.alpha[2][1], i.alpha[2][2], \

                       i.beta[0][0][0], i.beta[0][0][1], i.beta[0][0][2], \
                       i.beta[0][1][0], i.beta[0][1][1], i.beta[0][1][2], \
                       i.beta[0][2][0], i.beta[0][2][1], i.beta[0][2][2], \
                       i.beta[1][0][0], i.beta[1][0][1], i.beta[1][0][2], \
                       i.beta[1][1][0], i.beta[1][1][1], i.beta[1][1][2], \
                       i.beta[1][2][0], i.beta[1][2][1], i.beta[1][2][2], \
                       i.beta[2][0][0], i.beta[2][0][1], i.beta[2][0][2], \
                       i.beta[2][1][0], i.beta[2][1][1], i.beta[2][1][2], \
                       i.beta[2][2][0], i.beta[2][2][1], i.beta[2][2][2] )) )
        f_.close()

#Write the mol file for target cluster, if the corresponding qua_ file
#already exist( i.e. calculation has been run, perform analysis on those)

    if write_mol:
        if args.x.endswith(".pdb"):
            name = args.x.split(".")[0] + "_" + str(args.waters) + ".mol"
        elif args.x.endswith( ".xyz" ):
            name = args.x.split(".")[0] + ".mol"

        f_ = open( name , "w" )
        f_.write( "ATOMBASIS\n\nComment\nAtomtypes=2 Charge=0 Nosymm Angstrom\n")
        if not f_waters:
            "Can't write to .mol file, didn't read water molecules"
            raise SystemExit

        hCnt = len(waters) * 2
        oCnt = len(waters)
        f_.write( "Charge=1.0 Atoms=%d Basis=cc-pVDZ\n" % hCnt)

        for i in waters:
            for j in i.atomlist:
                if j.atype == "H":
                    j.toAA()
                    f_.write( "%s   %.5f   %.5f   %.5f\n" %( j.atype, j.x, j.y, j.z ))
                    j.toAU()

        f_.write( "Charge=8.0 Atoms=%d Basis=cc-pVDZ\n" % oCnt)
        for i in waters:
            for j in i.atomlist:
                if j.atype == "O":
                    j.toAA()
                    f_.write( "%s   %.5f   %.5f   %.5f\n" %( j.atype, j.x, j.y, j.z ))
                    j.toAU()
        raise SystemExit

#Write the xyz file for target cluster
    if write_xyz:
        f_ = open( args.ox , "w" )
        f_.write( str( 3* len(waters)) + "\n\n" )
        if f_waters:
            for i in waters:
                for j in i.atomlist:
                    f_.write( "%s %.5f %.5f %.5f\n" %( j.atype, j.x, j.y, j.z ))
# Read in QM dipole moment from dalton .out files if they exist for supplied .xyz/.pdb file
# lin_tip3p3_10.out is the corresponding out file for tip3p md configuration 3 using 10
# water molecules obtained as linear response with dipole moment
#
    #if args.x:
    #    if args.x.endswith( ".pdb" ):
    #        lin_outfile = "lin_" + args.x.split('.')[0]+"_" + str( args.waters ) + ".out"
    #        qua_outfile = "qua_" + args.x.split('.')[0]+"_" + str( args.waters ) + ".out"
    #    elif args.x.endswith( ".xyz" ):
    #        lin_outfile = "lin_" + args.x.split('.')[0] + ".out"
    #        qua_outfile = "qua_" + args.x.split('.')[0] + ".out"
    #    if os.path.isfile( lin_outfile ):
    #        atoms, qm_dipole = read_coords_and_dipole( args, custom_file = lin_outfile )
    #        print "Found QM dipole moment in %s" %lin_outfile
    #        print qm_dipole

# Do olav calculations for the generated .pot file for the supplied .xyz/.pdb file:
    if not f_waters:
        raise SystemExit

    #if args.op:
    #    string  =  open( args.op ).read()

    print string_hyperpolarizable

    if f_waters:
        static = PointDipoleList.from_string( string_static )
        polarizable = PointDipoleList.from_string( string_polarizable )
        hyperpolarizable = PointDipoleList.from_string( string_hyperpolarizable )

        static.solve_scf()
        polarizable.solve_scf()
        hyperpolarizable.solve_scf()

#Dipole section

    select = [ (0, 0, 2), (1, 1, 2), (2, 2, 2)]
    ref = [ dipole_qm, alpha_qm, beta_qm ]

    c =  Calculator()
    c.writeLog()


    reference = dipole_qm

    print '\n\nDipole: p_x, p_y, p_z\n'
    print "Quantum Mech: ", reference
    print "static dipole", static.total_dipole_moment()
    print "polarizable dipole", polarizable.total_dipole_moment()
    print "hyperpolarizable dipole", hyperpolarizable.total_dipole_moment()

    print "\n"
    print "Relative error Dipole static:" , \
            [(this-ref)/ref for this, ref in zip(
                    static.total_dipole_moment(),
                            reference )]
    print "Relative error Dipole polarizable:" , \
         [(this-ref)/ref for this, ref in zip(
                    polarizable.total_dipole_moment(),
                            reference )]
    print "Relative error Dipole hyperpolarizable:" , \
        [(this-ref)/ref for this, ref in zip(
                    hyperpolarizable.total_dipole_moment(),
                            reference )]
#Alpha section
    reference = alpha_qm.diagonal()
    print '\n\nAlfa: a_xx, a_yy, a_zz\n'
    print "Quantum Mech: ", reference
    print "static alpha", static.alpha().diagonal()
    print "polarizable alpha", polarizable.alpha().diagonal()
    print "hyperpolarizable alpha", hyperpolarizable.alpha().diagonal()

    print "\n"
    print "Relative error Alpha polarizable:" ,\
            [(this-ref)/ref for this, ref in zip(
                polarizable.alpha().diagonal(), 
                        reference )]
    print "Relative error Alpha hyperpolarizable:",\
            [(this-ref)/ref for this, ref  \
            in zip( hyperpolarizable.alpha().diagonal(), 
                        reference )]



#Beta section
#Relative error for xxz, yyz, zzz

    select = [ (0, 0, 2), (1, 1, 2), (2, 2, 2)]
    reference = [beta_qm[i, j, k] for i, j, k in select]
    print "Relative error for xxz, yyz, zzz"
    print '\n\nBeta: b_xxz, b_yyz, b_zzz\n'
    print "Quantum Mech: ", reference
    print "Static:", [static.beta()[i, j, k] for i, j, k in select]
    print "Polarizable:" , [polarizable.beta()[i, j, k] for i, j, k in select]
    print "Hyperpolarizable:" ,[hyperpolarizable.beta()[i, j, k] for i, j, k in select]

    print "\n"
    print "Relative error Beta:" ,[(this-ref)/ref for this, ref in zip([
                             hyperpolarizable.beta()[i, j, k] for i, j, k in select],
                        reference)]


if __name__ == '__main__':
    main()
