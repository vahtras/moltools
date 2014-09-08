#!/usr/bin/env python
#-*- coding: utf-8 -*-

import os, re, subprocess
import numpy as np
import math as m
import argparse

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt

from particles import *
from quadrupole import *
from gaussian import *

a0 = 0.52917721092
charge_dic = {"H": 1.0, "C": 6.0, "N": 7.0, "O": 8.0, "S": 16.0}
mass_dict = {"H": 1.008,  "C": 6.0, "N": 7.0, "O": 15.999, "S": 16.0}
index_dict = { 
        "static" : {"dipole":{"X": [0,0] , "Y": [0,1], "Z": [0,2]}},

        "polar" : {"dipole":{"X": [1,0] , "Y": [1,1], "Z": [1,2]},
                   "alpha":{"X": [3,0] , "Y": [3,1], "Z": [3,2]} },
        "hyper" : { "dipole":{"X": [2,0], "Y": [2,1], "Z": [2,2]},
                   "alpha":{"X": [4,0] , "Y": [4,1], "Z": [4,2] } ,
                   "beta":{"X": [5,0] , "Y": [5,1], "Z": [5,2]} }}
line_thick_dict = { 
        "static" : {"dipole":{"X": 0, "Y": 0, "Z": 0}},

        "polar" : {"dipole":{"X": 2 , "Y": 2, "Z": 2},
                   "alpha":{"X":  2 , "Y": 2, "Z": 2} },
        "hyper" : { "dipole":{"X": 3 , "Y": 3, "Z": 3 },
                   "alpha":{  "X": 3 , "Y": 3, "Z": 3 },
                   "beta":{   "X": 3 , "Y": 3, "Z": 3 }}}
line_style_dict = { 
        "static" : {"dipole":{"X": 1, "Y": 1, "Z": 1}},

        "polar" : {"dipole":{"X": 2 , "Y": 2, "Z": 2},
                   "alpha":{"X":  2 , "Y": 2, "Z": 2} },
        "hyper" : { "dipole":{"X": 3, "Y": 3, "Z": 3 },
                   "alpha":{  "X": 3, "Y": 3, "Z": 3 },
                   "beta":{   "X": 3, "Y": 3, "Z": 3 }}}


class Template:
    def __init__(self):
        monomer1_hf_cc_p_cdz  =  \
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

        monomer2_hf_cc_p_cdz  =  \
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




        """ 
        
        Coordinates for this model: 
        theta = 104.5
        r = 0.972

        """
        olav_hf_cc_p_vdz =  \
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

        centered_b3lyp_cc_p_vdz =  {}

#Template properties for CENTERED, r = 0.958019, theta = 104.5 water
        centered_hf_cc_p_vdz =  \
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

#  TIP3P model HF cc-p_vdz
        tip3p_hf_cc_p_vdz =  \
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
        spc_hf_cc_p_vdz =  \
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
        olav_hf_dict = { "PVDZ" : olav_hf_cc_p_vdz }
        centered_b3lyp_dict = { "MIDDLE" : centered_b3lyp_middle,
                "PVDZ" : centered_b3lyp_cc_p_vdz  }
        centered_hf_dict = { "PVDZ" : centered_hf_cc_p_vdz }
        tip3p_b3lyp_dict = { "MIDDLE" : tip3p_b3lyp_middle }
        tip3p_hf_dict = { "PVDZ" : tip3p_hf_cc_p_vdz }
        spc_b3lyp_dict = { "MIDDLE" : spc_b3lyp_middle }
        spc_hf_dict = { "PVDZ" : spc_hf_cc_p_vdz  }

        monomer1_hf_dict = { "PVDZ" : monomer1_hf_cc_p_cdz }
        monomer2_hf_dict = { "PVDZ" : monomer2_hf_cc_p_cdz }

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

        self.name_dict = { "OLAV" : olav_method_dict ,
                "MON1" : monomer1_method_dict,
                "MON2" : monomer2_method_dict,
                "CENTERED" : centered_method_dict, \
                "TIP3P": tip3p_method_dict, \
                "SPC" : spc_method_dict }
    def get_data(self, model, method, basis):
        return self.name_dict[model][method][basis]

class Dic:
    def __init__(self):
        rho3 = { 0.00 : [] }
        rho2 = { 0.00 : rho3 }
        rho1 = { 0.00 : rho2 }
        tau = { 0.00 : rho1}
        theta = { 0.00 : tau }
        r = { 0.00 : theta }
        self.dic = r
    def get_val(self, r, theta, tau, rho1, rho2, rho3):
        if self.dic.has_key( r ):
            if self.dic[r].has_key( theta ):
                if self.dic[r][theta].has_key(tau):
                    if self.dic[r][theta][tau].has_key(rho1):
                        if self.dic[r][theta][tau][rho1].has_key(rho2):
                            if self.dic[r][theta][tau][rho1][rho2].has_key(rho3):
                                return self.dic[r][theta][tau][rho1][rho2][rho3]
    def set_val(self, r, theta, tau, rho1, rho2, rho3, val ): #, rho1, rho2, rho3, val):
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
            self.options = { "is_aa" : True  }

        self.dic = Dic()

        self.r_oh = 0.97167
        self.theta_hoh = np.pi * 104.5/ 180.0

        self.vary_r =     False
        self.vary_theta = False
        self.vary_tau =   False
        self.vary_rho1 =  False
        self.vary_rho2 =  False
        self.vary_rho3 =  False

        self.options_r =     {  "min": 5.00, "max" : 5.00, "points": 1 }
        self.options_theta = {  "min": 0.00, "max" : 0.00, "points": 1 }
        self.options_tau =   {  "min": 0.00, "max" : 0.00, "points": 1 }
        self.options_rho1 =  {  "min": 0.00, "max" : 0.00, "points": 1 }
        self.options_rho2 =  {  "min": 0.00, "max" : 0.00, "points": 1 }
        self.options_rho3 =  {  "min": 0.00, "max" : 0.00, "points": 1 }

        opts = { "r" : {"min":2.30, "max":5.00,  "points":100} ,
             "theta" : {"max": np.pi , "min": 0.00 , "points":10},
             "tau"  : {"max": np.pi/2 , "min": 0.00 , "points":10},
             "rho1" : {"max": np.pi , "min": 0.00 , "points":1},
             "rho2" : {"max": np.pi , "min": 0.00 , "points":1},
             "rho3" : {"max": np.pi , "min": 0.00 , "points":1},
             }
        self.vary_parameters()

    def gen_mols(self):
        r = np.r_[ self.options_r[ "min" ] : self.options_r[ "max" ] : \
                complex( "%sj"%self.options_r[ "points" ] ) ]
        theta = np.r_[ self.options_theta[ "min" ] : self.options_theta[ "max" ] : \
                complex( "%sj"%self.options_theta[ "points" ] ) ]
        tau = np.r_[ self.options_tau[ "min" ] : self.options_tau[ "max" ] : \
                complex( "%sj"%self.options_tau[ "points" ] ) ]
        rho1 = np.r_[ self.options_rho1[ "min" ] : self.options_rho1[ "max" ] : \
                complex( "%sj"%self.options_rho1[ "points" ] ) ]
        rho2 = np.r_[ self.options_rho2[ "min" ] : self.options_rho2[ "max" ] : \
                complex( "%sj"%self.options_rho2[ "points" ] ) ]
        rho3 = np.r_[ self.options_rho3[ "min" ] : self.options_rho3[ "max" ] : \
                complex( "%sj"%self.options_rho3[ "points" ] ) ]
        for i in r:
            for j in theta:
                for k in tau:
                    for l in rho1:
                        for m in rho2:
                            for n in rho3:
                                w1 = self.get_water( [0, 0, 0], self.r_oh, self.theta_hoh)
                                x, y, z = self.get_cartesian_from_degree( i, j, k )
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
                        self.vary_r = True
                        self.options_r = j[i]
                    if i == "theta":
                        self.vary_theta = True
                        self.options_theta = j[i]
                    if i == "tau":
                        self.vary_tau = True
                        self.options_tau = j[i]
                    if i == "rho1":
                        self.vary_rho1 = True
                        self.options_rho1 = j[i]
                    if i == "rho2":
                        self.vary_rho2 = True
                        self.options_rho2 = j[i]
                    if i == "rho3":
                        self.vary_rho3 = True
                        self.options_rho3 = j[i]

    def get_water(self, origin, r, theta, AA = True ):
        h1 = Atom() ; h2 = Atom() ; o = Atom()
        d = (m.pi/2 - theta/2)
        o.element = "O" ; h1.element = "H" ; h2.element = "H"
        o.x = origin[0] ; o.y = origin[1] ; o.z = origin[2] 
        h1.x = (origin[0] + r * m.cos(d)) ; h1.y = origin[1] ; h1.z = (origin[2] + r*m.sin(d))
        h2.x = (origin[0] - r * m.cos(d)) ; h2.y = origin[1] ; h2.z = (origin[2] + r*m.sin(d))
        w = Water(); w.add_atom( o) ;w.add_atom( h2 ) ;w.add_atom( h1 ) 
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
    def read_waters(self, fname):
        """From file with name fname, return a list of all waters encountered"""
#If the file is plain xyz file

        atoms = []
        if fname.endswith( ".xyz" ) or fname.endswith(".mol"):
            pat_xyz = re.compile(r'^\s*(\w+)\s+(-*\d*.+\d+)\s+(-*\d*.+\d+)\s+(-*\d*.+\d+) *$')
            for i in open( fname ).readlines():
                if pat_xyz.match(i):
                    f = pat_xyz.match(i).groups()
                    tmp_atom = Atom()
                    tmp_atom.AA = True
                    tmp_atom.x = float(f[1])
                    tmp_atom.y = float(f[2])
                    tmp_atom.z = float(f[3])
                    tmp_atom.element = f[0][0]
                    atoms.append( tmp_atom )

        elif fname.endswith( ".pdb" ):
            pat1 = re.compile(r'^(ATOM|HETATM)')
            for i in open( fname ).readlines():
                if pat1.search(i):
                    #Ignore charge centers for polarizable water models
                    if ( i[11:16].strip() == "SW") or (i[11:16] == "DW"):
                        continue
                    tmp_atom = Atom(i[11:16].strip()[0], \
                            float(i[30:38].strip()), \
                            float(i[38:46].strip()), \
                            float(i[46:54].strip()), \
                            int(i[22:26].strip()) )

                    if fname_aaor_au == "AU":
                        if args.op_aaor_au == "AA":
                            tmp_atom.to_aa()
                    elif fname_aaor_au == "AA":
                        if args.op_aaor_au == "AU":
                            tmp_atom.to_au()
                    atoms.append( tmp_atom )
        elif fname.endswith( ".out" ):
            pat_xyz = re.compile(r'^(\w+)\s+(-*\d*.+\d+)\s+(-*\d*.+\d+)\s+(-*\d*.+\d+) *$')
            for i in open( fname ).readlines():
                if pat_xyz.match(i):
                    f = pat_xyz.match(i).groups()
                    tmp_atom = Atom(f[0][0], float(f[1]), float(f[2]), float(f[3]), 0)
                    if fname_aaor_au == "AU":
                        if args.op_aaor_au == "AA":
                            tmp_atom.to_aa()
                    elif fname_aaor_au == "AA":
                        if args.op_aaor_au == "AU":
                            tmp_atom.to_au()
                    atoms.append( tmp_atom )
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
                tmp.add_atom( i )
                for j in atoms:
                    if j.element == "O":
                        continue
                    if j.in_water:
                        continue
#If in cartesian:
                    if j.AA:
                        if i.dist_to_atom(j) < 1.1:
                            tmp.add_atom ( j )
                            j.in_water = True
                    else:
                        if i.dist_to_atom(j) < 1.1/a0:
                            tmp.add_atom ( j )
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
#__Water__.add_atom() method will update the waters residue number and center coordinate
#When all atoms are there
#Right now NOT center-of-mass
                tmp.add_atom(i)
                for j in atoms:
                    if j.element != "H":
                        continue
                    if j.in_water:
                        continue
#1.05 because sometimes spc water lengths can be over 1.01
                        
                    if args.op_aaor_au == "AA":
                        if i.dist(j) <= 1.05:
                            j.in_water = True
                            tmp.add_atom( j )
                            if len(tmp.atomlist) == 3:
                                break
                    elif args.op_aaor_au == "AU":
                        if i.dist(j) <= 1.05/a0:
                            j.in_water = True
                            tmp.add_atom( j )
                            if len(tmp.atomlist) == 3:
                                break
                wlist.append( tmp )
            wlist.sort( key = lambda x: x.dist_to_point( center ))
            center_water = wlist[0]
            cent_wlist = wlist[1:]
            cent_wlist.sort( key= lambda x: x.dist_to_water( center_water) )
            waters = [center_water] + cent_wlist[ 0:args.waters - 1 ]
        elif fname.endswith( ".out" ):
            for i in atoms:
                if i.element == "H":
                    continue
                if i.in_water:
                    continue
                tmp = Water(  )
                i.in_water = True
                tmp.add_atom( i )
                for j in atoms:
                    if j.element == "O":
                        continue
                    if j.in_water:
                        continue
#If in cartesian:
                    if i.AA:
                        if i.dist(j) < 1.0:
                            tmp.add_atom ( j )
                            j.in_water = True
                    else:
                        if i.dist(j) < 1.0/a0:
                            tmp.add_atom ( j )
                            j.in_water = True
                tmp.number = cnt
                cnt += 1
                waters.append( tmp )
        return waters
    def write_mol(self, wlist, name = "tmp.mol" ):
        f_ = open (name, 'w')
        f_.write("ATOMBASIS\n\n\n_atomtypes=2 Charge=0 Angstrom Nosymm\n")
        f_.write("Charge=1.0 Atoms=4 Basis=cc-p_vdz\n")
        for i in wlist:
            for j in i.atomlist:
                if j.element != "H":
                    continue
                f_.write( "%s %.5f %.5f %.5f\n" %(j.element, j.x, j.y, j.z ) )
        f_.write("Charge=8.0 Atoms=2 Basis=cc-p_vdz\n")
        for i in wlist:
            for j in i.atomlist:
                if j.element != "O":
                    continue
                f_.write( "%s %.5f %.5f %.5f\n" %(j.element, j.x, j.y, j.z ) )
        f_.close()

    def get_cartesian_from_degree(self, r, theta, tau):
        return r* m.sin( m.pi*theta/180.0 )*m.cos( m.pi*tau/180.0) \
               , r* m.sin( m.pi*theta/180.0 )*m.sin( m.pi*tau/180.0)  \
               , r* m.cos( m.pi*theta/180.0 )

class Atom(object):
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
        self.in_water = False

        self.AA = True

    def __str__(self):
        return "%s %f %f %f" %(self.element, self.x, self.y, self.z)
    def __sub__(self, other ):
        """return numpy array between this and other atom"""
        return self.get_array() - other.get_array()
    def get_array(self):
        return np.array( [self.x , self.y, self.z ] ).copy()
    def dist_to_atom(self, other):
        return np.sqrt( (self.x - other.x)**2 + (self.y -other.y)**2 + (self.z -other.z)**2 )
    def to_au(self):
        if self.AA:
            self.x /= a0
            self.y /= a0
            self.z /= a0
            self.AA = False
    def to_aa(self):
        if not self.AA:
            self.x *= a0
            self.y *= a0
            self.z *= a0
            self.AA = True

class Water(list):
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

        self.b_no_hydrogens = True
        self.h1 = False
        self.h2 = False
        self.o  = False
#By default put center on oxygen, put center of mass if args.com
        self.center  = False
        self.res_id = 0
        self.atomlist  = []

        self.AA = True
    def add_atom(self, atom):
        if self.atoms > 3:
            print "tried to add additional atoms to water, exiting"
            raise System_exit
        self.atoms += 1
        if atom.element == "H":
            if self.b_no_hydrogens:
                self.h1 = atom
                self.b_no_hydrogens = False
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

        if self.res_id:
            if self.res_id != atom.resid:
                print "Tried to add %s to %s, exiting" %(atom, self)
                raise System_exit
        else:
#Initialize water res_id from atomic resid
            self.res_id = atom.resid

    def __str__(self):
        return str("WAT") + str(self.res_id) 

    def get_angle_rho(self, other):
        d1 = self.get_dipole()
        d2 = other.get_dipole()
        return np.arccos( np.dot( d1, d2)/ ( np.linalg.norm(d1) * np.linalg.norm(d2) ) )

    def get_angle_tau(self, other):
        r1= self.get_norm()
        r2= other.get_norm()
        return np.arccos( np.dot( r1, r2 ) / (np.linalg.norm( r1 ) * np.linalg.norm( r2 )))

    def get_dipole(self):
        hq = 0.25
        oq = -0.5
        return self.h1.get_array() * hq + self.h2.get_array() * hq + self.o.get_array() * oq

    def get_norm(self):
        r1 = self.h1 - self.o
        r2 = self.h2 - self.o
        return np.cross( r1, r2 )

    def dist_to_point( self , point ):
        return m.sqrt( (self.center[0] - point[0])**2 +\
                (self.center[1] - point[1])**2  + ( self.center[2] -point[2])**2 )

    def dist_to_water(self, other):
        xyz1 = self.center
        xyz2 = other.center
        return m.sqrt( (xyz1[0] - xyz2[0])**2 + \
            (xyz1[1] - xyz2[1])**2 + (xyz1[2] - xyz2[2])**2 )

    def get_euler(self):
        """Return euler angles required to rotate water in oxygen at origo to current"""

        H1 = self.h1.get_array()
        H2 = self.h2.get_array()
        O1 = self.o.get_array()

        dip = self.get_dipole()

        origin = O1.copy()
        H1, H2, O1 = H1 - origin, H2 - origin, O1 - origin

        theta1 = m.atan2( dip[1], dip[0])
        H1 =  np.dot( self.get_rzinv( theta1 ) , H1 )
        H2 =  np.dot( self.get_rzinv( theta1 ) , H2 )
        O1 =  np.dot( self.get_rzinv( theta1 ) , O1 )
        dip = np.dot( self.get_rzinv( theta1 ) , dip )
#Rotate by theta around y axis so that the dipole is in the z axis 
        theta2 = m.atan2( -dip[0], dip[2])
        H1 =  np.dot( self.get_ry( theta2 ) , H1 )
        H2 =  np.dot( self.get_ry( theta2 ) , H2 )
        O1 =  np.dot( self.get_ry( theta2 ) , O1 )
        dip = np.dot( self.get_ry( theta2 ) , dip )
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
        first Rz with theta1, then Ry^-1 by theta2, then Rz with theta 3

        R all in radians

        """
        d1, d2, d3 = self.get_euler()
        
# Place water molecule in origo, and rotate it so hydrogens in xz plane
        H1 = self.h1.get_array() ; H2 = self.h2.get_array() ; O = self.o.get_array()
        TMP = self.o.get_array()
        H1 -= TMP ; H2 -= TMP; O -= TMP

        H1 = np.dot( self.get_rzinv(d3) , H1 )
        H1 = np.dot( self.get_ry(d2) , H1 )
        H1 = np.dot( self.get_rzinv(d1) , H1 )

        H2 = np.dot( self.get_rzinv(d3) , H2 )
        H2 = np.dot( self.get_ry(d2) , H2 )
        H2 = np.dot( self.get_rzinv(d1) , H2 )

        O = np.dot( self.get_rzinv(d3) , O )
        O = np.dot( self.get_ry(d2) , O )
        O = np.dot( self.get_rzinv(d1) , O )

# Rotate with angles t1, t2, t3

        H1 = np.dot( self.get_rz(t1) , H1 )
        H1 = np.dot( self.get_ryinv(t2) , H1 )
        H1 = np.dot( self.get_rz(t3) , H1 )

        H2 = np.dot( self.get_rz(t1) , H2 )
        H2 = np.dot( self.get_ryinv(t2) , H2 )
        H2 = np.dot( self.get_rz(t3) , H2 )

        O = np.dot( self.get_rz(t1) , O )
        O = np.dot( self.get_ryinv(t2) , O )
        O = np.dot( self.get_rz(t3) , O )

#Put back in oxygen original point
        H1 += TMP ; H2 += TMP; O += TMP

        self.h1.x = H1[0] ;self.h1.y = H1[1] ;self.h1.z = H1[2] 
        self.h2.x = H2[0] ;self.h2.y = H2[1] ;self.h2.z = H2[2] 
        self.o.x  =  O[0] ;  self.o.y = O[1] ;  self.o.z = O[2] 
    def get_rz( self, theta ):
        vec = np.array( [[ m.cos(theta), -m.sin(theta), 0],
                            [ m.sin(theta), m.cos(theta), 0],
                            [ 0,    0,  1]])
        return vec
    def get_rzinv( self, theta ):
        vec = np.array(     [[ m.cos(theta), m.sin(theta), 0],
                            [ -m.sin(theta), m.cos(theta), 0],
                            [ 0,             0,            1]])
        return vec
    def get_ry( self, theta ):
        vec = np.array( [[ m.cos(theta),0, m.sin(theta)],
                            [ 0,    1,  0],
                            [ -m.sin(theta), 0, m.cos(theta)]])
        return vec
    def get_ryinv( self, theta ):
        vec = np.array( [[ m.cos(theta),0, -m.sin(theta)],
                            [ 0,    1,  0],
                            [ m.sin(theta), 0, m.cos(theta)]])
        return vec

    def plot_water(self ):
#Plot water molecule in green and  nice xyz axis
        O1, H1, H2 = self.o, self.h1, self.h2
        fig = plt.figure()
        dip = self.get_dipole()
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

    def transform_dipole( self, qmdipole, t1, t2, t3 ):

        d_new1 = np.zeros([3]) #will be returned
        d_new2 = np.zeros([3]) #will be returned
        d_new3 = np.zeros([3]) #will be returned

        rz  = self.get_rz( t1 )
        ryi = self.get_ryinv( t2 )
        rz2 = self.get_rz( t3 )

        for i in range(3):
            for x in range(3):
                d_new1[i] += rz[i][x] * qmdipole[x]
        for i in range(3):
            for x in range(3):
                d_new2[i] += ryi[i][x] * d_new1[x]
        for i in range(3):
            for x in range(3):
                d_new3[i] += rz2[i][x] * d_new2[x]
        self.qm_dipole = d_new3
        return d_new3

    def transform_alpha( self, qmalpha, t1, t2 , t3 ):

        a_new1 = np.zeros([3,3]) #will be calculated
        a_new2 = np.zeros([3,3]) #will be calculated
        a_new3 = np.zeros([3,3]) #will be calculated

        rz = self.get_rz( t1 )
        ryi = self.get_ryinv( t2 )
        rz2 = self.get_rz( t3 )

        #print 'inside transfer alpha, water: %d at x = %f' %( self.res_id, self.o.x)
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

    def transform_beta( self, qmbeta, t1, t2, t3 ):

        b_new1 = np.zeros([3,3,3]) #will be calculated
        b_new2 = np.zeros([3,3,3]) #will be calculated
        b_new3 = np.zeros([3,3,3]) #will be calculated

        rz =  self.get_rz( t1 )
        ryi = self.get_ryinv( t2 )
        rz2 = self.get_rz( t3 )

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

    def get_square_dipole(self):
        return np.sqrt( self.qm_dipole[0] **2 + self.qm_dipole[1]**2 + self.qm_dipole[2]**2 )

    def get_square_beta(self):
        return np.sqrt( \
                (self.qm_beta[0][0][0] + self.qm_beta[0][1][1] + self.qm_beta[0][2][2] )**2 + \
                (self.qm_beta[1][0][0] + self.qm_beta[1][1][1] + self.qm_beta[1][2][2] )**2 + \
                (self.qm_beta[2][0][0] + self.qm_beta[2][1][1] + self.qm_beta[2][2][2] )**2  )

    def get_alpha_trace(self):
        return  self.qm_alpha[0][0] + self.qm_alpha[1][1] + self.qm_alpha[2][2]

    def to_au(self):
        if self.AA:
            self.h1.to_au()
            self.h2.to_au()
            self.o.to_au()
            self.AA = False
    def to_aa(self):
        if not self.AA:
            self.h1.to_aa()
            self.h2.to_aa()
            self.o.to_aa()
            self.AA = True
            
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
                y.append( self.Dict.get_val( r, theta, tau, rho1, rho2, rho3 ) )
            if self.opts["theta"].has_key( "vary" ):
                x.append( theta )
                y.append( self.Dict.get_val( r, theta, tau, rho1, rho2, rho3 ) )
            if self.opts["tau"].has_key( "vary" ):
                x.append( tau )
                y.append( self.Dict.get_val( r, theta, tau, rho1, rho2, rho3 ) )
            if self.opts["rho1"].has_key( "vary" ):
                x.append( rho1 )
                y.append( self.Dict.get_val( r, theta, tau, rho1, rho2, rho3 ) )
            if self.opts["rho2"].has_key( "vary" ):
                x.append( rho2 )
                y.append( self.Dict.get_val( r, theta, tau, rho1, rho2, rho3 ) )
            if self.opts["rho3"].has_key( "vary" ):
                x.append( rho3 )
                y.append( self.Dict.get_val( r, theta, tau, rho1, rho2, rho3 ) )

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

            qm_dipole = self.get_qm_dipole( "hfqua_" + i + ".out" )
            qm_alpha = self.get_qm_alpha(  "hfqua_" + i + ".out" )
            qm_beta = self.get_qm_beta(   "hfqua_" + i + ".out" )
            tmp_waters = []

            for j in self.read_waters( i + ".mol" ):
                t1, t2, t3 =  j.get_euler()
                tmp_dipole, tmp_alpha, tmp_beta = Template().get_data( "OLAV", "HF", "PVDZ" )
                j.dipole = j.transform_dipole( tmp_dipole, t1, t2, t3 )
                j.alpha = j.transform_alpha( tmp_alpha, t1, t2, t3 )
                j.beta = j.transform_beta( tmp_beta, t1, t2, t3 )
                tmp_waters.append( j )

#
            if args.model == "dipole":
                static= PointDipoleList.from_string( self.get_static_string( tmp_waters ))
                polar = PointDipoleList.from_string( self.get_polar_string(  tmp_waters ))
                hyper = PointDipoleList.from_string( self.get_hyper_string(  tmp_waters ))

            if args.model == "gaussian":
                static= GaussianQuadrupoleList.from_string( self.get_static_string( tmp_waters ))
                polar = GaussianQuadrupoleList.from_string( self.get_polar_string(  tmp_waters ))
                hyper = GaussianQuadrupoleList.from_string( self.get_hyper_string(  tmp_waters ))

                tmp = float(args.stdev)
                for j in static:
                    j._R_p = tmp
                for j in polar:
                    j._R_p = tmp
                for j in hyper:
                    j._R_p = tmp

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
                print self.get_static_string( tmp_waters )

            d_static =  \
                     [(this-ref)/ref for this, ref in zip(
                            static.total_dipole_moment(),
                                    qm_dipole )] 
            d_polar =  \
                     [(this-ref)/ref for this, ref in zip(
                            polar.total_dipole_moment(),
                                    qm_dipole )] 
            d_hyper = \
                     [(this-ref)/ref for this, ref in zip(
                            hyper.total_dipole_moment(),
                                    qm_dipole )] 
            a_polar = \
                     [(this-ref)/ref for this, ref in zip(
                            polar.alpha().diagonal(),
                                    qm_alpha.diagonal() )] 
            a_hyper = \
                     [(this-ref)/ref for this, ref in zip(
                            hyper.alpha().diagonal(),
                                    qm_alpha.diagonal() )] 

            select = [ (0, 0, 2), (1, 1, 2), (2, 2, 2)]
            reference = [ qm_beta[ii, jj, kk] for ii, jj, kk in select ]

            b_hyper =  [(this-ref)/ref for this, ref in zip( [ hyper.beta()[ii, jj, kk] for ii, jj, kk in select], reference  ) ] 

            val = [ d_static, d_polar, d_hyper, a_polar, a_hyper, b_hyper ]
            if args.param:
                r, theta, tau, rho1, rho2, rho3 = i.split('-')
                self.Dict.set_val( r, theta, tau, rho1, rho2, rho3, val)

    def read_waters(self, fname):
        """From file with name fname, return a list of all waters encountered"""
#If the file is plain xyz file

        atoms = []
        if fname.endswith( ".xyz" ) or fname.endswith(".mol"):
            pat_xyz = re.compile(r'^\s*(\w+)\s+(-*\d*.+\d+)\s+(-*\d*.+\d+)\s+(-*\d*.+\d+) *$')
            for i in open( fname ).readlines():
                if pat_xyz.match(i):
                    f = pat_xyz.match(i).groups()
                    tmp_atom = Atom()
                    tmp_atom.AA = True
                    tmp_atom.x = float(f[1])
                    tmp_atom.y = float(f[2])
                    tmp_atom.z = float(f[3])
                    tmp_atom.element = f[0][0]
                    atoms.append( tmp_atom )

        elif fname.endswith( ".pdb" ):
            pat1 = re.compile(r'^(ATOM|HETATM)')
            for i in open( fname ).readlines():
                if pat1.search(i):
                    #Ignore charge centers for polarizable water models
                    if ( i[11:16].strip() == "SW") or (i[11:16] == "DW"):
                        continue
                    tmp_atom = Atom(i[11:16].strip()[0], \
                            float(i[30:38].strip()), \
                            float(i[38:46].strip()), \
                            float(i[46:54].strip()), \
                            int(i[22:26].strip()) )

                    if fname_aaor_au == "AU":
                        if args.op_aaor_au == "AA":
                            tmp_atom.to_aa()
                    elif fname_aaor_au == "AA":
                        if args.op_aaor_au == "AU":
                            tmp_atom.to_au()
                    atoms.append( tmp_atom )
        elif fname.endswith( ".out" ):
            pat_xyz = re.compile(r'^(\w+)\s+(-*\d*.+\d+)\s+(-*\d*.+\d+)\s+(-*\d*.+\d+) *$')
            for i in open( fname ).readlines():
                if pat_xyz.match(i):
                    f = pat_xyz.match(i).groups()
                    tmp_atom = Atom(f[0][0], float(f[1]), float(f[2]), float(f[3]), 0)
                    if fname_aaor_au == "AU":
                        if args.op_aaor_au == "AA":
                            tmp_atom.to_aa()
                    elif fname_aaor_au == "AA":
                        if args.op_aaor_au == "AU":
                            tmp_atom.to_au()
                    atoms.append( tmp_atom )
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
                tmp.add_atom( i )
                for j in atoms:
                    if j.element == "O":
                        continue
                    if j.in_water:
                        continue
#If in cartesian:
                    if j.AA:
                        if i.dist_to_atom(j) < 1.1:
                            tmp.add_atom ( j )
                            j.in_water = True
                    else:
                        if i.dist_to_atom(j) < 1.1/a0:
                            tmp.add_atom ( j )
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
#__Water__.add_atom() method will update the waters residue number and center coordinate
#When all atoms are there
#Right now NOT center-of-mass
                tmp.add_atom(i)
                for j in atoms:
                    if j.element != "H":
                        continue
                    if j.in_water:
                        continue
#1.05 because sometimes spc water lengths can be over 1.01
                        
                    if args.op_aaor_au == "AA":
                        if i.dist(j) <= 1.05:
                            j.in_water = True
                            tmp.add_atom( j )
                            if len(tmp.atomlist) == 3:
                                break
                    elif args.op_aaor_au == "AU":
                        if i.dist(j) <= 1.05/a0:
                            j.in_water = True
                            tmp.add_atom( j )
                            if len(tmp.atomlist) == 3:
                                break
                wlist.append( tmp )
            wlist.sort( key = lambda x: x.dist_to_point( center ))
            center_water = wlist[0]
            cent_wlist = wlist[1:]
            cent_wlist.sort( key= lambda x: x.dist_to_water( center_water) )
            waters = [center_water] + cent_wlist[ 0:args.waters - 1 ]
        elif fname.endswith( ".out" ):
            for i in atoms:
                if i.element == "H":
                    continue
                if i.in_water:
                    continue
                tmp = Water(  )
                i.in_water = True
                tmp.add_atom( i )
                for j in atoms:
                    if j.element == "O":
                        continue
                    if j.in_water:
                        continue
#If in cartesian:
                    if i.AA:
                        if i.dist(j) < 1.0:
                            tmp.add_atom ( j )
                            j.in_water = True
                    else:
                        if i.dist(j) < 1.0/a0:
                            tmp.add_atom ( j )
                            j.in_water = True
                tmp.number = cnt
                cnt += 1
                waters.append( tmp )
        return waters
    def get_matching_out_and_mol(self):
        tmp = []
        tmp_out = self.get_out_files()
        tmp_found = {}
        cnt = 0
        for j in self.get_mol_files():
            if os.path.isfile( os.path.join( os.getcwd(),  "hfqua_" + j.rstrip(".mol") + ".out" ) ):
                tmp.append( j.rstrip( ".mol" ))
                continue
        return tmp
    def get_mol_files(self):
        mol_files = [f for f in os.listdir( os.getcwd() ) \
            if f.endswith( ".mol")  ]
        return mol_files

    def get_out_files(self):
        out_files = [f for f in os.listdir( os.getcwd() ) \
            if f.endswith( ".out")  ]
        return out_files
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
                tmp_atom = Atom()
                tmp_atom.AA = True
                tmp_atom.element = f[0][0]
                tmp_atom.x = float(f[1]); tmp_atom.y = float(f[2]); tmp_atom.z = float(f[3])
                tmp_atom.to_au()
                atoms.append( tmp_atom )
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
        local_counter = 0
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
                        in1, in2 =  index_dict[ level ][ prop ][ component ]
                        width = line_thick_dict[ level] [ prop ] [ component ]
                        style = line_style_dict[ level] [ prop ] [ component ]
                    except KeyError:
                        print "Skipping (%s, %s, %s, )" %( level, prop, component )
                        continue

                    string += '@TITLE "Relative errors as a function of %s"\n' \
                        % Xms(var).make_greek() 

                    string += '@SUBTITLE "Using: %s (%s) ;' \
                        %( args.model, args.stdev )

                    if args.vary_r:
                        #string += '@SUBTITLE "Constant %s: %s, %s: %s, %s: %s, %s: %s, %s: %s" \n'\
                        string += ' Constant %s: %s, %s: %s, %s: %s, %s: %s, %s: %s" \n'\
                                %(
                                  Xms("theta").make_greek(), self.opts["theta"]["constant"],
                                  Xms("tau").make_greek(),   self.opts["tau"]["constant"],
                                  Xms("rho1").make_greek(),  self.opts["rho1"]["constant"],
                                  Xms("rho2").make_greek(),  self.opts["rho2"]["constant"],
                                  Xms("rho3").make_greek(),  self.opts["rho3"]["constant"])
                    if args.vary_theta:
                        #string += '@SUBTITLE "Constant %s: %s, %s: %s, %s: %s, %s: %s, %s: %s" \n'\
                        string += ' Constant %s: %s, %s: %s, %s: %s, %s: %s, %s: %s" \n'\
                                %(
                                  Xms("r").make_greek(), self.opts["r"]["constant"],
                                  Xms("tau").make_greek(), self.opts["tau"]["constant"],
                                  Xms("rho1").make_greek(), self.opts["rho1"]["constant"],
                                  Xms("rho2").make_greek(), self.opts["rho2"]["constant"],
                                  Xms("rho3").make_greek(), self.opts["rho3"]["constant"])
                    if args.vary_tau:
                        #string += '@SUBTITLE "Constant %s: %s, %s: %s, %s: %s, %s: %s, %s: %s" \n'\
                        string += ' Constant %s: %s, %s: %s, %s: %s, %s: %s, %s: %s" \n'\
                                %(
                                  Xms("r").make_greek(), self.opts["r"]["constant"],
                                  Xms("theta").make_greek(), self.opts["theta"]["constant"],
                                  Xms("rho1").make_greek(), self.opts["rho1"]["constant"],
                                  Xms("rho2").make_greek(), self.opts["rho2"]["constant"],
                                  Xms("rho3").make_greek(), self.opts["rho3"]["constant"])
                    if args.vary_rho1:
                        #string += '@SUBTITLE "Constant %s: %s, %s: %s, %s: %s, %s: %s, %s: %s" \n'\
                        string += ' Constant %s: %s, %s: %s, %s: %s, %s: %s, %s: %s" \n'\
                                %(
                                  Xms("r").make_greek(), self.opts["r"]["constant"],
                                  Xms("theta").make_greek(), self.opts["tau"]["constant"],
                                  Xms("tau").make_greek(), self.opts["tau"]["constant"],
                                  Xms("rho2").make_greek(), self.opts["rho2"]["constant"],
                                  Xms("rho3").make_greek(), self.opts["rho3"]["constant"])
                    if args.vary_rho2:
                        #string += '@SUBTITLE "Constant %s: %s, %s: %s, %s: %s, %s: %s, %s: %s" \n'\
                        string += ' Constant %s: %s, %s: %s, %s: %s, %s: %s, %s: %s" \n'\
                                %(
                                  Xms("r").make_greek(), self.opts["r"]["constant"],
                                  Xms("theta").make_greek(), self.opts["tau"]["constant"],
                                  Xms("tau").make_greek(), self.opts["tau"]["constant"],
                                  Xms("rho1").make_greek(), self.opts["rho1"]["constant"],
                                  Xms("rho3").make_greek(), self.opts["rho3"]["constant"])
                    if args.vary_rho3:
                        #string += '@SUBTITLE "Constant %s: %s, %s: %s, %s: %s, %s: %s, %s: %s" \n'\
                        string += ' Constant %s: %s, %s: %s, %s: %s, %s: %s, %s: %s" \n'\
                                %(Xms("r").make_greek(), self.opts["r"]["constant"],
                                  Xms("theta").make_greek(), self.opts["theta"]["constant"],
                                  Xms("tau").make_greek(), self.opts["tau"]["constant"],
                                  Xms("rho1").make_greek(), self.opts["rho1"]["constant"],
                                  Xms("rho2").make_greek(), self.opts["rho2"]["constant"])

                    string +=  '@VIEW 0.15, 0.08, 1.15, 0.85\n'
                    string +=  '@LEGEND ON\n'
                    string +=  '@LEGEND BOX ON\n'
                    string +=  '@LEGEND BOX FILL OFF\n'
                    string +=  '@LEGEND LOCTYPE VIEW\n'
                    string +=  '@LEGEND 0.50, 0.84\n' 
                    string +=  '@ s%d LEGEND "%s, %s, %s"\n' %( local_counter, level, prop, component) 
                    string +=  '@ XAXIS LABEL "%s"\n' % Xms( var ).make_greek()
                    string +=  '@ YAXIS LABEL "Relative error"\n' 

#add the actual data x and y

                    for i in range(len( x )):
                        string += "%s %.4f\n" %( x[i], y[i][in1][in2] )
                    string += '@ SORT s%d X ASCENDING\n' % local_counter
                    string += '@ s%d LINEWIDTH %d\n' % (local_counter, width )
                    string += '@ s%d LINESTYLE %d\n' % (local_counter, style )
                    local_counter += 1

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
            in1, in2 =  index_dict[ level ][ prop ][ component ]
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

    def get_static_string( self, waters, to_au = True ):
        """ Converts list of waters into Olav string for .pot"""
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
    def get_hyper_string( self, waters ):
        """ Converts list of waters into Olav string for .pot"""
        string_hyperpolarizable = "AU\n%d 1 22 1\n" %len(waters)
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
                tmp_atom = Atom(f[0][0], float(f[1]), float(f[2]), float(f[3]), 0)
                tmp_atom.to_au()
                atoms.append( tmp_atom )
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
            tmp.add_atom( i )
            for j in atoms:
                if j.element == "O":
                    continue
                if j.in_water:
                    continue
#If in cartesian:
                if i.dist(j) < 1.89:
                    tmp.add_atom ( j )
                    j.in_water = True
            tmp.number = cnt
            cnt += 1
            waters.append( tmp )
        return waters

class Xms:

    def __init__(self, char):
        self.char = char
    def make_greek(self):
        if self.char not in ["r", "theta", "tau", "rho1", "rho2", "rho3"]:
            print "wrong char in Xmgrace_style class, exiting"; raise System_exit
        if self.char == "r"    :return r"r"            
        if self.char == "theta":return r"\f{Symbol}q"  
        if self.char == "tau"  :return r"\f{Symbol}t" 
        if self.char == "rho1" :return r"\f{Symbol}r\s1\N" 
        if self.char == "rho2" :return r"\f{Symbol}r\s2\N" 
        if self.char == "rho3" :return r"\f{Symbol}r\s3\N" 





if __name__ == '__main__':

    A = argparse.ArgumentParser( add_help= True)

    A.add_argument( "-model"  , type = str, default = "dipole", choices = [ "dipole", "quadrupole", "gaussian"] )
    A.add_argument( "-stdev"  , type = str, default = "0.00001" )
    A.add_argument( "-param"  , default = True )


    A.add_argument( "-l", nargs = "*", default = ["static"] )
    A.add_argument( "-p", nargs = "*", default = ["dipole"] )
    A.add_argument( "-c", nargs = "*", default = ["Z"] )

    A.add_argument( "-vary_r"   , default = False , action = 'store_true' ) 
    A.add_argument( "-vary_theta" , default = False , action = 'store_true' )
    A.add_argument( "-vary_tau"   , default = False , action = 'store_true' ) 
    A.add_argument( "-vary_rho1", default = False , action = 'store_true' ) 
    A.add_argument( "-vary_rho2", default = False , action = 'store_true' ) 
    A.add_argument( "-vary_rho3", default = False , action = 'store_true' ) 

    A.add_argument( "-r"   , type = str ) 
    A.add_argument( "-theta" , type = str ) 
    A.add_argument( "-tau"   , type = str ) 
    A.add_argument( "-rho1", type = str ) 
    A.add_argument( "-rho2", type = str ) 
    A.add_argument( "-rho3", type = str ) 

    args = A.parse_args()

    c = Calculator()
    if args.r:
        c.opts[ "r" ] = { "constant" : args.r }
    if args.theta:
        c.opts[ "theta" ] = { "constant" : args.theta }
    if args.tau:
        c.opts[ "tau" ] = { "constant" : args.tau }
    if args.rho1:
        c.opts[ "rho1" ] = { "constant" : args.rho1 }
    if args.rho2:
        c.opts[ "rho2" ] = { "constant" : args.rho2 }
    if args.rho3:
        c.opts[ "rho3" ] = { "constant" : args.rho3 }

    if args.vary_r:
        c.opts[ "r" ] = { "vary" : True }
        args.var = "r"
    if args.vary_theta:
        c.opts[ "theta" ] = { "vary" : True }
        args.var = "theta"
    if args.vary_tau:
        c.opts[ "tau" ] = { "vary" : True }
        args.var = "tau"
    if args.vary_rho1:
        c.opts[ "rho1" ] = { "vary" : True }
        args.var = "rho1"
    if args.vary_rho2:
        c.opts[ "rho2" ] = { "vary" : True }
        args.var = "rho2"
    if args.vary_rho3:
        c.opts[ "rho3" ] = { "vary" : True }
        args.var = "rho3"
    c.get_rel_error( args )

    string = c.get_xvg_string( args )
    open('tmp.xvg', 'w').write( string )

    #c.write_log( level = args.l , prop = args.p,  component = args.c, variable = var )


