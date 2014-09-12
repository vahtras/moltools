#!/usr/bin/env python
#-*- coding: utf-8 -*-

import os, re, subprocess
import numpy as np
import  math as m
from particles import *
from quadrupole import *
from gaussian import *

from generator import *
from water import *
from template import *
from R import *
from owndict import *

a0 = 0.52917721092
charge_dic = {"H": 1.0, "C": 6.0, "N": 7.0, "O": 8.0, "S": 16.0}
mass_dict = {"H": 1.008,  "C": 6.0, "N": 7.0, "O": 15.999, "S": 16.0}
indexDict = { 
        "static" : {"dipole":{"X": [0,0] , "Y": [0,1], "Z": [0,2]}},

        "polar" : {"dipole":{"X": [1,0] , "Y": [1,1], "Z": [1,2]},
                   "alpha":{"X": [3,0] , "Y": [3,1], "Z": [3,2]} },
        "hyper" : { "dipole":{"X": [2,0], "Y": [2,1], "Z": [2,2]},
                   "alpha":{"X": [4,0] , "Y": [4,1], "Z": [4,2] } ,
                   "beta":{"X": [5,0] , "Y": [5,1], "Z": [5,2]} }}
lineThickDict = { 
        "static" : {"dipole":{"X": 0, "Y": 0, "Z": 0}},

        "polar" : {"dipole":{"X": 2 , "Y": 2, "Z": 2},
                   "alpha":{"X":  2 , "Y": 2, "Z": 2} },
        "hyper" : { "dipole":{"X": 3 , "Y": 3, "Z": 3 },
                   "alpha":{  "X": 3 , "Y": 3, "Z": 3 },
                   "beta":{   "X": 3 , "Y": 3, "Z": 3 }}}
lineStyleDict = { 
        "static" : {"dipole":{"X": 1, "Y": 1, "Z": 1}},

        "polar" : {"dipole":{"X": 2 , "Y": 2, "Z": 2},
                   "alpha":{"X":  2 , "Y": 2, "Z": 2} },
        "hyper" : { "dipole":{"X": 3, "Y": 3, "Z": 3 },
                   "alpha":{  "X": 3, "Y": 3, "Z": 3 },
                   "beta":{   "X": 3, "Y": 3, "Z": 3 }}}


class Calculator:
    """ Class container for retrieving and calculating data
    retrieving water classes from DALTON .mol / .pdb / .xyz
    retrieving water classes  from .pdb / .xyz
    retrieves properties from DALTON .out
    """

    def __init__(self):

        self.opts = { 
             "r"      : {"constant" : "5.00"  },
             "theta"  : {"constant" : "0.00" },
             "tau"    : {"constant" : "0.00" },
             "rho1"   : {"constant" : "0.00" },
             "rho2"   : {"constant" : "0.00" },
             "rho3"   : {"constant" : "0.00" }
             }

        self.Dict = Dic()

    def get_x_and_y( self ):
        x = []
        y = []
        for i in self.get_matching_out_and_mol():
            r, theta, tau, rho1, rho2, rho3 = i.split('-')
            if self.opts["r"].has_key( "constant" ):
                if r != self.opts["r"]["constant"]:
                    continue
            if self.opts["theta"].has_key( "constant" ):
                if theta != self.opts["theta"]["constant"]:
                    continue
            if self.opts["tau"].has_key( "constant" ):
                if tau != self.opts["tau"]["constant"]:
                    continue
            if self.opts["rho1"].has_key( "constant" ):
                if rho1 != self.opts["rho1"]["constant"]:
                    continue
            if self.opts["rho2"].has_key( "constant" ):
                if rho2 != self.opts["rho2"]["constant"]:
                    continue
            if self.opts["rho3"].has_key( "constant" ):
                if rho3 != self.opts["rho3"]["constant"]:
                    continue
            if self.opts["r"].has_key( "vary" ):
                x.append( r )
                y.append( self.Dict.getVal( r, theta, tau, rho1, rho2, rho3 ) )
            if self.opts["theta"].has_key( "vary" ):
                x.append( theta )
                y.append( self.Dict.getVal( r, theta, tau, rho1, rho2, rho3 ) )
            if self.opts["tau"].has_key( "vary" ):
                x.append( tau )
                y.append( self.Dict.getVal( r, theta, tau, rho1, rho2, rho3 ) )
            if self.opts["rho1"].has_key( "vary" ):
                x.append( rho1 )
                y.append( self.Dict.getVal( r, theta, tau, rho1, rho2, rho3 ) )
            if self.opts["rho2"].has_key( "vary" ):
                x.append( rho2 )
                y.append( self.Dict.getVal( r, theta, tau, rho1, rho2, rho3 ) )
            if self.opts["rho3"].has_key( "vary" ):
                x.append( rho3 )
                y.append( self.Dict.getVal( r, theta, tau, rho1, rho2, rho3 ) )

        return  x , y 

    def get_rel_error( self, args ):
        for i in self.get_matching_out_and_mol():
            r, theta, tau, rho1, rho2, rho3 = i.split('-')
            if self.opts["r"].has_key( "constant" ):
                if r != self.opts["r"]["constant"]:
                    continue
            if self.opts["theta"].has_key( "constant" ):
                if theta != self.opts["theta"]["constant"]:
                    continue
            if self.opts["tau"].has_key( "constant" ):
                if tau != self.opts["tau"]["constant"]:
                    continue
            if self.opts["rho1"].has_key( "constant" ):
                if rho1 != self.opts["rho1"]["constant"]:
                    continue
            if self.opts["rho2"].has_key( "constant" ):
                if rho2 != self.opts["rho2"]["constant"]:
                    continue
            if self.opts["rho3"].has_key( "constant" ):
                if rho3 != self.opts["rho3"]["constant"]:
                    continue

            qmDipole = self.get_qm_dipole( "hfqua_" + i + ".out" )
            qmAlpha = self.get_qm_alpha(  "hfqua_" + i + ".out" )
            qmBeta = self.get_qm_beta(   "hfqua_" + i + ".out" )
            tmp_waters = []

            for j in self.read_waters( i + ".mol" ):
                t1, t2, t3 =  j.get_euler()

                kwargs = Template().get_dist_data( "OLAV", "HF", "PVDZ" )

                if args.model: # == "gaussian":

                    p = Property.from_dist_template( **kwargs )
                    j.Property = p

                else:
                    j.dipole = j.transform_dist_dipole( t_dipole, t1, t2, t3 )
                    j.quadrupole = j.transform_dist_quadrupole( t_quad, t1, t2, t3 )
                    j.alpha = j.transform_dist_alpha( t_alpha, t1, t2, t3 )
                    j.beta = j.transform_dist_beta( t_beta, t1, t2, t3 )

                tmp_waters.append( j )
#
            if args.model == "pointdipole":
                pass
                #static= PointDipoleList.from_string( self.get_static_string( tmp_waters ),
                #        max_l = 1, pol = 0, hyper = 0, dist = True ))

                #polar = PointDipoleList.from_string( self.get_polar_string(  tmp_waters ),
                #        max_l = 1, pol = 2, hyper = 1, dist = True ))

                #hyper = PointDipoleList.from_string( self.get_hyper_string(  tmp_waters,
                #        max_l = 1, pol = 22, hyper = 1, dist = True ))

            if args.model == "gaussian":

                tmp_Rq = float(args.Rq)
                tmp_Rp = float(args.Rp)

                if "static" in args.l:
                    static = GaussianQuadrupoleList.from_string( self.get_string(  tmp_waters,
                        max_l =1, pol = 0 , hyper = 0, dist = True ))
                    for j in static:
                        j._R_q = tmp_Rq
                        j._R_p = tmp_Rp
                if "polar" in args.l:
                    polar = GaussianQuadrupoleList.from_string( self.get_string(  tmp_waters,
                        max_l =2, pol = 2 , hyper = 1, dist = True ))
                    for j in polar:
                        j._R_q = tmp_Rq
                        j._R_p = tmp_Rp
                if "hyper" in args.l:
                    hyper = GaussianQuadrupoleList.from_string( self.get_string(  tmp_waters,
                        max_l =2, pol = 22, hyper = 1, dist = True ))
                    for j in hyper:
                        j._R_q = tmp_Rq
                        j._R_p = tmp_Rp
                print hyper
                raise SystemExit


            if args.model == "quadrupole":
                static= QuadrupoleList.from_string( self.get_static_string( tmp_waters ))
                polar = QuadrupoleList.from_string( self.get_polar_string(  tmp_waters ))
                hyper = QuadrupoleList.from_string( self.get_hyper_string(  tmp_waters ))

            try:
                static.solve_scf()
                polar.solve_scf()
                hyper.solve_scf()
            except:
                print i
                #print self.get_static_string( tmp_waters )

            d_static = []
            d_polar = []
            d_hyper = []

            a_polar = []
            a_hyper = []

            b_hyper = []

            if "static" in args.l:
                d_static =  \
                         [(this-ref)/ref for this, ref in zip(
                                static.total_dipole_moment(),
                                        qmDipole )] 
            if "polar" in args.l:
                d_polar =  \
                         [(this-ref)/ref for this, ref in zip(
                                polar.total_dipole_moment(),
                                        qmDipole )] 
                a_polar = \
                         [(this-ref)/ref for this, ref in zip(
                                polar.alpha().diagonal(),
                                        qmAlpha.diagonal() )] 
            if "hyper" in args.l:
                d_hyper = \
                         [(this-ref)/ref for this, ref in zip(
                                hyper.total_dipole_moment(),
                                        qmDipole )] 
                a_hyper = \
                         [(this-ref)/ref for this, ref in zip(
                                hyper.alpha().diagonal(),
                                        qmAlpha.diagonal() )] 

                select = [ (0, 0, 2), (1, 1, 2), (2, 2, 2)]
                reference = [ qmBeta[ii, jj, kk] for ii, jj, kk in select ]

                b_hyper =  [(this-ref)/ref for this, ref in zip( [ hyper.beta()[ii, jj, kk] for ii, jj, kk in select], reference  ) ] 

            val = [ d_static, d_polar, d_hyper, a_polar, a_hyper, b_hyper ]

            if args.param:
                r, theta, tau, rho1, rho2, rho3 = i.split('-')
                self.Dict.setVal( r, theta, tau, rho1, rho2, rho3, val)

    def read_waters(self, fname):
        """From file with name fname, return a list of all waters encountered"""
#If the file is plain xyz file

        atoms = []
        if fname.endswith( ".xyz" ) or fname.endswith(".mol"):
            pat_xyz = re.compile(r'^\s*(\w+)\s+(-*\d*.+\d+)\s+(-*\d*.+\d+)\s+(-*\d*.+\d+) *$')
            for i in open( fname ).readlines():
                if pat_xyz.match(i):
                    f = pat_xyz.match(i).groups()
                    tmpAtom = Atom()
                    tmpAtom.x = float(f[1])
                    tmpAtom.y = float(f[2])
                    tmpAtom.z = float(f[3])
                    tmpAtom.element = f[0][0]
                    tmpAtom.isAA = True
                    #tmpAtom.to_au()
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
                            tmpAtom.to_au()
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
                            tmpAtom.to_au()
                    atoms.append( tmpAtom )
#loop over oxygen and hydrogen and if they are closer than 1 A add them to a water
        waters = []
        cnt = 1

        if fname.endswith( ".xyz" ) or fname.endswith(".mol"):
            for i in atoms:
                if i.element == "H":
                    continue
                if i.in_water:
                    continue
                tmp = Water(  )
                i.in_water = True
                tmp.append( i )
                for j in atoms:
                    if j.element == "O":
                        continue
                    if j.in_water:
                        continue
#If in cartesian:
                    if j.AA:
                        if i.dist_to_atom(j) < 1.1:
                            tmp.append ( j )
                            j.in_water = True
                    else:
                        if i.dist_to_atom(j) < 1.1/a0:
                            tmp.append ( j )
                            j.in_water = True
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
                i.in_water= True
#__Water__.append() method will update the waters residue number and center coordinate
#When all atoms are there
#Right now NOT center-of-mass
                tmp.append(i)
                for j in atoms:
                    if j.element != "H":
                        continue
                    if j.in_water:
                        continue
#1.05 because sometimes spc water lengths can be over 1.01
                        
                    if args.opAAorAU == "AA":
                        if i.dist(j) <= 1.05:
                            j.in_water = True
                            tmp.append( j )
                            if len(tmp.atomlist) == 3:
                                break
                    elif args.opAAorAU == "AU":
                        if i.dist(j) <= 1.05/a0:
                            j.in_water = True
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
                if i.in_water:
                    continue
                tmp = Water(  )
                i.in_water = True
                tmp.append( i )
                for j in atoms:
                    if j.element == "O":
                        continue
                    if j.in_water:
                        continue
#If in cartesian:
                    if i.AA:
                        if i.dist(j) < 1.0:
                            tmp.append ( j )
                            j.in_water = True
                    else:
                        if i.dist(j) < 1.0/a0:
                            tmp.append ( j )
                            j.in_water = True
                tmp.number = cnt
                cnt += 1
                waters.append( tmp )
        for wat in waters:
            for atom in wat:
                atom.res_id = wat.number
        return waters
    def get_matching_out_and_mol(self):
        tmp = []
        tmpOut = self.get_out_files()
        tmpFound = {}
        cnt = 0
        for j in self.get_mol_files():
            if os.path.isfile( os.path.join( os.getcwd(),  "hfqua_" + j.rstrip(".mol") + ".out" ) ):
                tmp.append( j.rstrip( ".mol" ))
                continue
        return tmp
    def get_mol_files(self):
        molFiles = [f for f in os.listdir( os.getcwd() ) \
            if f.endswith( ".mol")  ]
        return molFiles

    def get_out_files(self):
        outFiles = [f for f in os.listdir( os.getcwd() ) \
            if f.endswith( ".out")  ]
        return outFiles
    def get_qm_dipole( self, fname ):
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
                tmpAtom.to_au()
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
    def get_qm_alpha( self, fname):
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
    def get_qm_beta(self, fname ):
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
#Alpha section

    def get_xvg_string( self, args ):
        """ kwargs is a dictionary where user specifies which components, dipole models and properties
        will be returned and printed for a finished formatted .xvg xmgrace plot"""
        x, y = self.get_x_and_y()
        localCounter = 0
        string = ""


        for level in args.l:
            for prop in args.p:
                for component in args.c:

                    kwargs = { "level" : level, "prop" : prop,  "component" : component,"variable" : args.var }
                    level = kwargs.get( "level" , "hyper" )
                    prop = kwargs.get( "prop", "dipole" )
                    component = kwargs.get( "component", "X" )
                    var = kwargs.get( "variable", "r" )

                    try:
                        in1, in2 =  indexDict[ level ][ prop ][ component ]
                        width = lineThickDict[ level] [ prop ] [ component ]
                        style = lineStyleDict[ level] [ prop ] [ component ]
                    except KeyError:
                        print "Skipping (%s, %s, %s, )" %( level, prop, component )
                        continue

                    string += '@TITLE "Relative errors as a function of %s"\n' \
                        % Xms(var).makeGreek() 

                    string += '@SUBTITLE "Using Rp, Rq: (%s, %s) ;' \
                        %( args.Rq, args.Rp )

                    if args.vary_r:
                        #string += '@SUBTITLE "Constant %s: %s, %s: %s, %s: %s, %s: %s, %s: %s" \n'\
                        string += ' Constant %s: %s, %s: %s, %s: %s, %s: %s, %s: %s" \n'\
                                %(
                                  Xms("theta").makeGreek(), self.opts["theta"]["constant"],
                                  Xms("tau").makeGreek(),   self.opts["tau"]["constant"],
                                  Xms("rho1").makeGreek(),  self.opts["rho1"]["constant"],
                                  Xms("rho2").makeGreek(),  self.opts["rho2"]["constant"],
                                  Xms("rho3").makeGreek(),  self.opts["rho3"]["constant"])
                    if args.vary_theta:
                        #string += '@SUBTITLE "Constant %s: %s, %s: %s, %s: %s, %s: %s, %s: %s" \n'\
                        string += ' Constant %s: %s, %s: %s, %s: %s, %s: %s, %s: %s" \n'\
                                %(
                                  Xms("r").makeGreek(), self.opts["r"]["constant"],
                                  Xms("tau").makeGreek(), self.opts["tau"]["constant"],
                                  Xms("rho1").makeGreek(), self.opts["rho1"]["constant"],
                                  Xms("rho2").makeGreek(), self.opts["rho2"]["constant"],
                                  Xms("rho3").makeGreek(), self.opts["rho3"]["constant"])
                    if args.vary_tau:
                        #string += '@SUBTITLE "Constant %s: %s, %s: %s, %s: %s, %s: %s, %s: %s" \n'\
                        string += ' Constant %s: %s, %s: %s, %s: %s, %s: %s, %s: %s" \n'\
                                %(
                                  Xms("r").makeGreek(), self.opts["r"]["constant"],
                                  Xms("theta").makeGreek(), self.opts["theta"]["constant"],
                                  Xms("rho1").makeGreek(), self.opts["rho1"]["constant"],
                                  Xms("rho2").makeGreek(), self.opts["rho2"]["constant"],
                                  Xms("rho3").makeGreek(), self.opts["rho3"]["constant"])
                    if args.vary_rho1:
                        #string += '@SUBTITLE "Constant %s: %s, %s: %s, %s: %s, %s: %s, %s: %s" \n'\
                        string += ' Constant %s: %s, %s: %s, %s: %s, %s: %s, %s: %s" \n'\
                                %(
                                  Xms("r").makeGreek(), self.opts["r"]["constant"],
                                  Xms("theta").makeGreek(), self.opts["tau"]["constant"],
                                  Xms("tau").makeGreek(), self.opts["tau"]["constant"],
                                  Xms("rho2").makeGreek(), self.opts["rho2"]["constant"],
                                  Xms("rho3").makeGreek(), self.opts["rho3"]["constant"])
                    if args.vary_rho2:
                        #string += '@SUBTITLE "Constant %s: %s, %s: %s, %s: %s, %s: %s, %s: %s" \n'\
                        string += ' Constant %s: %s, %s: %s, %s: %s, %s: %s, %s: %s" \n'\
                                %(
                                  Xms("r").makeGreek(), self.opts["r"]["constant"],
                                  Xms("theta").makeGreek(), self.opts["tau"]["constant"],
                                  Xms("tau").makeGreek(), self.opts["tau"]["constant"],
                                  Xms("rho1").makeGreek(), self.opts["rho1"]["constant"],
                                  Xms("rho3").makeGreek(), self.opts["rho3"]["constant"])
                    if args.vary_rho3:
                        #string += '@SUBTITLE "Constant %s: %s, %s: %s, %s: %s, %s: %s, %s: %s" \n'\
                        string += ' Constant %s: %s, %s: %s, %s: %s, %s: %s, %s: %s" \n'\
                                %(Xms("r").makeGreek(), self.opts["r"]["constant"],
                                  Xms("theta").makeGreek(), self.opts["theta"]["constant"],
                                  Xms("tau").makeGreek(), self.opts["tau"]["constant"],
                                  Xms("rho1").makeGreek(), self.opts["rho1"]["constant"],
                                  Xms("rho2").makeGreek(), self.opts["rho2"]["constant"])

                    string +=  '@VIEW 0.15, 0.08, 1.15, 0.85\n'
                    string +=  '@LEGEND ON\n'
                    string +=  '@LEGEND BOX ON\n'
                    string +=  '@LEGEND BOX FILL OFF\n'
                    string +=  '@LEGEND LOCTYPE VIEW\n'
                    string +=  '@LEGEND 0.50, 0.84\n' 
                    string +=  '@ s%d LEGEND "%s, %s, %s"\n' %( localCounter, level, prop, component) 
                    string +=  '@ XAXIS LABEL "%s"\n' % Xms( var ).makeGreek()
                    string +=  '@ YAXIS LABEL "Relative error"\n' 

#add the actual data x and y

                    for i in range(len( x )):
                        string += "%s %.4f\n" %( x[i], y[i][in1][in2] )
                    string += '@ SORT s%d X ASCENDING\n' % localCounter
                    string += '@ s%d LINEWIDTH %d\n' % (localCounter, width )
                    string += '@ s%d LINESTYLE %d\n' % (localCounter, style )
                    localCounter += 1

        return string

    def write_log(self, **kwargs):
        #p = subprocess.Popen( 'rm *.log', shell=True )
        x, y = self.get_x_and_y()

        if kwargs is None:
            f = open( "rel_error.log" , 'w' )
            for i in range(len( x )):
                f.write( "%s %.4f\n" %( x[i], y[i][0][2] ) )
        else:
            level = kwargs.get( "level" , "hyper" )
            prop = kwargs.get( "prop", "dipole" )
            component = kwargs.get( "component", "X" )
            var = kwargs.get( "variable", "r" )
            in1, in2 =  indexDict[ level ][ prop ][ component ]
            f = open( "%s_%s_%s_%s.log" %(level, prop, component, var) , 'w')

            if prop == "dipole":
                f.write('@TITLE "Dipole moment error as a function of %s"\n' %var )
            if prop == "alpha":
                f.write('@TITLE "Alpha error as a function of %s"\n' %var )
            if prop == "beta":
                f.write('@TITLE "Beta error as a function of %s"\n' %var )

            f.write( '@LEGEND ON\n@LEGEND BOX ON\n@LEGEND LOCTYPE VIEW\n@LEGEND 0.80, 0.80\n' )
            f.write( '@ s0 LEGEND "%s, %s"\n' %(level, prop) )
            f.write( '@ XAXIS LABEL "%s"\n' %var )
            f.write( '@ YAXIS LABEL "Relative error"\n' )


#Write out the actual data x and y
            for i in range(len( x )):
                f.write( "%s %.4f\n" %( x[i], y[i][in1][in2] ))

            f.write( '@ SORT s0 X ASCENDING\n' )

            print "wrote %s_%s_%s_%s.log" %(level, prop, component, var) 

        f.close()

    def get_static_string( self, waters, to_au = True, dist = False ):
        """ Converts list of waters into Olav string for .pot"""

        if dist:
            pass
        else:
            for i in waters:
                if to_au:
                    i.to_au()
                i.center = [ i.o.x, i.o.y, i.o.z ]
            string_static = "AU\n%d 1 0\n" %len(waters)
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
        return string_static
    def get_polar_string( self, waters ):
        """ Converts list of waters into Olav string for .pot"""
        string_polarizable = "AU\n%d 1 2 1\n" %len(waters)
        for i in waters:
            i.center = [ i.o.x, i.o.y, i.o.z ]
        for i in waters:
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

        return string_polarizable

    def get_string( self, waters, max_l = 1, pol =22 , hyper = 2, dist = False ):
        """ Converts list of waters into Olav string for hyperpolarizable .pot"""
# If the properties are in distributed form, I. E. starts from Oxygen, then H in +x and H -x
        if dist:
            string_hyperpolarizable = "AA\n%d %d %d %d\n" % ( len(waters)*3,
                    max_l, pol, hyper )
            for i in waters:

                i.set_property_on_each_atom()
                string_hyperpolarizable +=  i.o.potline()
                string_hyperpolarizable +=  i.h1.potline()
                string_hyperpolarizable +=  i.h2.potline()
            return string_hyperpolarizable

        else:
            string_hyperpolarizable = "AU\n%d 1 22 2\n" %len(waters)
            for i in waters:
                i.center = [ i.o.x, i.o.y, i.o.z ]
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
        return string_hyperpolarizable

    def get_waters(self, f):
        """ f is a .mol file, will return array of water molecules """
        atoms = []
        pat_xyz = re.compile(r'^\s*(\w+)\s+(-*\d*.+\d+)\s+(-*\d*.+\d+)\s+(-*\d*.+\d+) *$')
        for i in open( f ).readlines():
            if pat_xyz.match(i):
                f = pat_xyz.match(i).groups()
                tmpAtom = Atom(f[0][0], float(f[1]), float(f[2]), float(f[3]), 0)
                tmpAtom.AA = True
                #tmpAtom.to_au()
                atoms.append( tmpAtom )
#loop over oxygen and hydrogen and if they are closer than 1 A add them to a water
        waters = []
        cnt = 1
        for i in atoms:
            if i.element == "H":
                continue
            if i.in_water:
                continue
            tmp = Water(  )
            i.in_water = True
            i.res_id = cnt
            tmp.append( i )
            for j in atoms:
                if j.element == "O":
                    continue
                if j.in_water:
                    continue
#If in cartesian:
                if i.dist(j) < 1.89:
                    tmp.append ( j )
                    j.in_water = True
            tmp.number = cnt
            cnt += 1
            waters.append( tmp )
        return waters

class Xms:

    def __init__(self, char):
        self.char = char
    def makeGreek(self):
        if self.char not in ["r", "theta", "tau", "rho1", "rho2", "rho3"]:
            print "wrong char in XmgraceStyle class, exiting"; raise SystemExit
        if self.char == "r"    :return r"r"            
        if self.char == "theta":return r"\f{Symbol}q"  
        if self.char == "tau"  :return r"\f{Symbol}t" 
        if self.char == "rho1" :return r"\f{Symbol}r\s1\N" 
        if self.char == "rho2" :return r"\f{Symbol}r\s2\N" 
        if self.char == "rho3" :return r"\f{Symbol}r\s3\N" 

if __name__ == '__main__':

    c = Calculator()
    g = Generator()
    r = R()
    w1 = g.getWater( [0, 0, 0] , 1.0 , 104.5/180 * m.pi )

