#!/usr/bin/env python
#-*- coding: utf-8 -*-

import os, re, subprocess
import numpy as np
import  math as m
from particles import *

a0 = 0.52917721092
charge_dic = {"H": 1.0, "C": 6.0, "N": 7.0, "O": 8.0, "S": 16.0}
mass_dict = {"H": 1.008,  "C": 6.0, "N": 7.0, "O": 15.999, "S": 16.0}

class Calc:
    def __init__(self):
        self.files = [f.rstrip(".mol") for f in os.listdir( os.getcwd() ) \
            if f.endswith( ".mol") ]
        self.molFiles = [f for f in os.listdir( os.getcwd() ) \
            if f.endswith( ".mol")  ]

    def getDipoleError(self):
        rel_static_dip =  \
                 [(this-ref)/ref for this, ref in zip(
                        self.static.total_dipole_moment(),
                                self.ref[0] )] 
        rel_polar_dip =  \
                 [(this-ref)/ref for this, ref in zip(
                        self.polarizable.total_dipole_moment(),
                                self.ref[0] )] 
        rel_hyp_dip =  \
                 [(this-ref)/ref for this, ref in zip(
                        self.hyperpolarizable.total_dipole_moment(),
                                self.ref[0] )] 
        return  rel_static_dip, rel_polar_dip, rel_hyp_dip
#Alpha section
    def getAlphaError(self):
        reference = self.ref[1].diagonal()
        rel_polar_alpha =  [(this-ref)/ref for this, ref in zip(  self.polarizable.alpha().diagonal() , reference  ) ]   
        rel_hyp_alpha =  [(this-ref)/ref for this, ref in zip(  self.hyperpolarizable.alpha().diagonal() , reference  ) ]  

        return rel_polar_alpha , rel_hyp_alpha
    def getBetaError(self):
        select = [ (0, 0, 2), (1, 1, 2), (2, 2, 2)]
        reference = [ self.ref[2][i, j, k] for i, j, k in select ]
        rel_hyp_beta =  [(this-ref)/ref for this, ref in zip( [ self.hyperpolarizable.beta()[i, j, k] for i, j, k in select], reference  ) ] 

        return rel_hyp_beta
    def writeLog(self):
        p = subprocess.Popen( 'rm *.log', shell=True )
        for i in self.files:
            theta , tau = i.split('-')[1] , i.split('-')[2]
            f = open( "rel_%s.log" % tau , 'a' )

            atoms, dip_qm , alpha_qm , beta_qm =  self.getQM( "hfqua_" + i + ".out" )
            self.ref = range(3)
            self.ref[0] = dip_qm
            self.ref[1] = alpha_qm
            self.ref[2] = beta_qm

            tmp = Templates().getBeta( "CENTERED", "HF", "PVDZ" )
            dipole = np.array( tmp[0] )
            alpha = np.array(  tmp[1] )
            beta = np.array(   tmp[2] )

            waters =  self.getWaters( i + ".mol" )
            for i in waters:

                i.dipole = dipole
                i.alpha = alpha
                i.beta = beta
                i.getEuler()

                i.transfer_dipole()
                i.transfer_alpha()
                i.transfer_beta()

            strings = self.getStrings( waters )

            self.static = PointDipoleList.from_string( strings[0] )
            self.polarizable = PointDipoleList.from_string( strings[1] )
            self.hyperpolarizable = PointDipoleList.from_string( strings[2] )

            self.static.solve_scf()
            self.polarizable.solve_scf()
            self.hyperpolarizable.solve_scf()

            rel_static_dip , rel_polar_dip , rel_hyp_dip  = self.getDipoleError()
            rel_polar_alpha, rel_hyp_alpha = self.getAlphaError()
            rel_hyp_beta  = self.getBetaError()

            #self.f_.write( base + "  " + rel_static_dip[2] + "\n")
            #self.f_.write( base + "  " + rel_polar_dip[2][2] + "\n")
            f.write( theta + "  " + str(rel_hyp_dip[2]) +"\n" )# \
            f.write( theta + "  " + str(rel_hyp_alpha[2]) +"\n" )# \
                    #+  " " + \
                    #str(rel_hyp_alpha[2]) + " " + \
                    #str(rel_hyp_beta[2]) + " " + \
                    #"\n") 
            f.close()

    def getQM(self, f_ ):
        """f is a output from dalton file"""
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
        for i in open( f_ ).readlines():
            if pat_xyz.match(i):
                f = pat_xyz.match(i).groups()
                tmpAtom = Atom(f[0][0], float(f[1]), float(f[2]), float(f[3]), 1 )
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
            nuc_dip[0] += charge_dic[ i.element ] * i.x
            nuc_dip[1] += charge_dic[ i.element ] * i.y
            nuc_dip[2] += charge_dic[ i.element ] * i.z





# Reading in Alfa and Beta tensor
        pat_alpha = re.compile(r'@ QRLRVE:.*([XYZ])DIPLEN.*([XYZ])DIPLEN')
        for i in open( f_ ).readlines():
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
        for i in open( f_ ).readlines():
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
    def getStrings( self, waters ):
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
        return string_static, string_polarizable, string_hyperpolarizable
    def getWaters(self, f):
        """ f is a .mol file, will return array of water molecules """

        atoms = []
        pat_xyz = re.compile(r'^\s*(\w+)\s+(-*\d*.+\d+)\s+(-*\d*.+\d+)\s+(-*\d*.+\d+) *$')
        for i in open( f ).readlines():
            if pat_xyz.match(i):
                f = pat_xyz.match(i).groups()
                tmpAtom = Atom(f[0][0], float(f[1]), float(f[2]), float(f[3]), 0)
                tmpAtom.toAU()
                atoms.append( tmpAtom )
#loop over oxygen and hydrogen and if they are closer than 1 A add them to a water
        waters = []
        cnt = 1
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
                if i.dist(j) < 1.89:
                    tmp.addAtom ( j )
                    j.inWater = True
            tmp.number = cnt
            cnt += 1
            waters.append( tmp )
        return waters
class Water:
    def __init__(self):
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
        """Return euler angles"""
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

        return theta3, theta2, theta1

    def rotate(self, t1, t2, t3):
        """Rotate self by t1, t2 and t3
        first Rz with theta1, then Ry^-1 by theta2, then Rz with theta 3"""
        d1, d2, d3 = self.getEuler()
        
# Place water molecule in origo
        H1 = self.h1.getArray() ; H2 = self.h2.getArray() ; O = self.o.getArray()
        TMP = O.copy()
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

# Rotate to defined angles t1 t2 t3

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
        self.o.x = O[0] ;self.o.y = O[1] ;self.o.z = O[2] 

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

        #print 'inside transfer alpha, water: %d at x = %f' %( self.resId, self.o.x)
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

    def __str__(self):
        return "%s %f %f %f" %(self.element, self.x, self.y, self.z)

    def __sub__(self, other ):
        """return numpy array between this and other atom"""
        return self.getArray() - other.getArray()

    def getArray(self):
        return np.array( [self.x , self.y, self.z ] ).copy()

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

class Atom:
    def __init__(self, element, x, y ,z, resid ):
        self.element = element
        self.x = x
        self.y = y
        self.z = z
        self.resid = resid
        self.inWater = False

    def __str__(self):
        return "%s %f %f %f" %(self.element, self.x, self.y, self.z)

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

if __name__ == '__main__':
    c = Calc()
    c.writeLog()
