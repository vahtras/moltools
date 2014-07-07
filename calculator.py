#!/usr/bin/env python
#-*- coding: utf-8 -*-

import os, re, subprocess
import numpy as np
import  math as m
from particles import *

from generator import *
from water import *
from template import *
from R import *

a0 = 0.52917721092
charge_dic = {"H": 1.0, "C": 6.0, "N": 7.0, "O": 8.0, "S": 16.0}
mass_dict = {"H": 1.008,  "C": 6.0, "N": 7.0, "O": 15.999, "S": 16.0}

class Calculator:
    """ Class container for retrieving and calculating data
    retrieving water classes from DALTON .mol / .pdb / .xyz
    retrieving water classes  from .pdb / .xyz
    retrieves properties from DALTON .out
    """

    def __init__(self):
        pass
    def getMolFiles(self):
        molFiles = [f for f in os.listdir( os.getcwd() ) \
            if f.endswith( ".mol")  ]
        return molFiles
    def getXyzFiles(self):
        xyzFiles = [f for f in os.listdir( os.getcwd() ) \
            if f.endswith( ".xyz")  ]
        return xyzFiles
    def getOutFiles(self):
        outFiles = [f for f in os.listdir( os.getcwd() ) \
            if f.endswith( ".out")  ]
        return outFiles
    def getQmDipole( self, fname ):
        """fname is an output file from DALTON calculation with *PROPAV"""
        nuc_dip = np.zeros(3)
        el_dip = np.zeros(3)
        atoms = []
        pat_xyz = re.compile(r'^\s*(\w+)\s+(-*\d*\.+\d+)\s+(-*\d*\.+\d+)\s+(-*\d*\.+\d+) *$')
        pat_pol = re.compile(r'([XYZ])DIPLEN.*total.*:')

# Reading in dipole
#
        for i in open( fname ).readlines():
            if pat_xyz.match(i):
                f = pat_xyz.match(i).groups()
                tmpAtom = Atom()
                tmpAtom.AA = True
                tmpAtom.element = f[0][0]
                tmpAtom.x = float(f[1]); tmpAtom.y = float(f[2]); tmpAtom.z = float(f[3])
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
        return nuc_dip - el_dip

    def getQmAlpha( self, fname):
        alpha = np.zeros([3,3])
        lab = ["X", "Y", "Z"]
# Reading in Alfa and Beta tensor
        pat_alpha = re.compile(r'@ QRLRVE:.*([XYZ])DIPLEN.*([XYZ])DIPLEN')
        for i in open( fname ).readlines():
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
    def getQmBeta(self, fname ):
        tmp = []
        missing = {}
        exists = {}
        lab = ["X", "Y", "Z"]
        beta = np.zeros([3,3,3])
        pat_beta = re.compile(r'@ B-freq')
        for i in open( fname ).readlines():
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

        return beta
    def getDipoleError( self ):
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
    def getAlphaError(self, reference):

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

    def getStrings( self, waters ):
        """ Converts list of waters into Olav string for .pot"""
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

if __name__ == '__main__':

    c = Calculator()
    g = Generator()
    r = R()
    w1 = g.getWater( [0, 0, 0] , 1.0 , 104.5/180 * m.pi )


    refDipole, refAlpha, refBeta = Template().getReference("TIP3P", "HF","PVDZ")

    for i in c.getOutFiles():
        for j in c.getMolFiles():
            if i.rstrip(".out").lstrip("hfqua_") == j.rstrip(".mol"):
                base = j.rstrip(".mol")
                for k in g.readWaters( j ):
                    print r.obj.FloatVector( xrange(4) )


