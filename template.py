#!/usr/bin/env python
#-*- coding: utf-8 -*-

class Template(dict):

    def __init__(self):

        self[ ("OLAV", "HF", "PVDZ", False, 0.0 ) ] = \
        {
#Charge
"charge" : [0.000], 
#Dipx     Dipy       Dipz
"dipole" :  [0.000, 0.000, 0.81458], 
#Quadxx    Quadxy    Quadxz    Quadyy     Quadyz    Quadzz
"quadrupole" : [ -3.443 , 0.000, 0.000, -5.254, 0.00 , -4.254],
#Alphaxx   Alphaxy   Alphaxz   Alphayy   Alphayz   Alphazz
"alpha" :   [7.21103, 0.00000, -0.00000, 3.03447,  -0.00000, 5.22711],
          #Betaxxx    Betaxxy    Betaxxz    Betaxyy     Betaxyz   Betaxzz   Betayyy     Betayyz    Betayzz    Betazzz
"beta" : [0.00000, 0.00000, -18.48240, -0.00000, -0.00000, -0.00000, 0.00000, -2.33726, 0.00000, -11.17619]
        }

        self[ ("OLAV", "HF", "PVDZ", True, 0.0 ) ] = \
        {
#Charge
"charge" : [-0.651, 0.325, 0.325 ],

#Dipx     Dipy       Dipz
"dipole" :  [[-0.000,     0.000,     0.342],
             [-0.167,     0.000,    -0.129],
             [ 0.167,     0.000,    -0.129]],

#Quadxx    Quadxy    Quadxz    Quadyy     Quadyz    Quadzz
"quadrupole" : [[-3.213,   -0.000,    -0.000,    -4.362,     0.000,    -3.822],
               [-0.115,    -0.000,     0.291,    -0.446,    -0.000,    -0.216],
               [-0.115,     0.000,    -0.291,    -0.446,    -0.000,    -0.216]],

#Alphaxx   Alphaxy   Alphaxz   Alphayy   Alphayz   Alphazz
"alpha" :   [[3.268,     0.000,      0.000,   1.531,     -0.000,     2.611],
             [1.971,    -0.000,      1.131,   0.752,     -0.000,     1.308],
             [1.971,     0.000,     -1.131,   0.752,     -0.000,     1.308]],

          #Betaxxx    Betaxxy    Betaxxz    Betaxyy     Betaxyz   Betaxzz   Betayyy     Betayyz    Betayzz    Betazzz
"beta" : [[-0.000,     0.000,    -8.567,     0.000,    -0.000,     0.000,     0.000,    -0.373,     0.000,    -4.010],
          [-7.465,     0.000,    -4.959,    -1.692,     0.000,    -4.101,     0.000,    -0.982,     0.000,    -3.583],
          [ 7.465,     0.000,    -4.959,     1.692,     0.000,     4.101,     0.000,    -0.982,     0.000,    -3.583]]

        }

        self.olav_hf_cc_pVDZ_dist = \
        {
#Distributed Dipole for O, H1, H2
#Charge
"charge" : [-0.651, 0.325, 0.325 ],

#Dipx     Dipy       Dipz
"dipole" :  [[-0.000,     0.000,     0.342],
             [-0.167,     0.000,    -0.129],
             [ 0.167,     0.000,    -0.129]],

#Quadxx    Quadxy    Quadxz    Quadyy     Quadyz    Quadzz
"quadrupole" : [[-3.213,   -0.000,    -0.000,    -4.362,     0.000,    -3.822],
               [-0.115,    -0.000,     0.291,    -0.446,    -0.000,    -0.216],
               [-0.115,     0.000,    -0.291,    -0.446,    -0.000,    -0.216]],

#Alphaxx   Alphaxy   Alphaxz   Alphayy   Alphayz   Alphazz
"alpha" :   [[3.268,     0.000,      0.000,   1.531,     -0.000,     2.611],
             [1.971,    -0.000,      1.131,   0.752,     -0.000,     1.308],
             [1.971,     0.000,     -1.131,   0.752,     -0.000,     1.308]],

          #Betaxxx    Betaxxy    Betaxxz    Betaxyy     Betaxyz   Betaxzz   Betayyy     Betayyz    Betayzz    Betazzz
"beta" : [[-0.000,     0.000,    -8.567,     0.000,    -0.000,     0.000,     0.000,    -0.373,     0.000,    -4.010],
          [-7.465,     0.000,    -4.959,    -1.692,     0.000,    -4.101,     0.000,    -0.982,     0.000,    -3.583],
          [ 7.465,     0.000,    -4.959,     1.692,     0.000,     4.101,     0.000,    -0.982,     0.000,    -3.583]]

        }

        self.olav_hf_cc_pVDZ = \
        {
#Charge
"charge" : 0.000, 

#Dipx     Dipy       Dipz
"dipole" :  [0.000, 0.000, 0.81458], 

#Quadxx    Quadxy    Quadxz    Quadyy     Quadyz    Quadzz
"quadrupole" : [ -3.443 , 0.000, 0.000, -5.254, 0.00 , -4.254],

#Alphaxx   Alphaxy   Alphaxz   Alphayy   Alphayz   Alphazz
"alpha" :   [7.21103, 0.00000, -0.00000, 3.03447,  -0.00000, 5.22711],

          #Betaxxx    Betaxxy    Betaxxz    Betaxyy     Betaxyz   Betaxzz   Betayyy     Betayyz    Betayzz    Betazzz
"beta" : [0.00000, 0.00000, -18.48240, -0.00000, -0.00000, -0.00000, 0.00000, -2.33726, 0.00000, -11.17619]
        }

        self[ ("OLAV", "HF", "PVDZ", True , "0.0" ) ] = \
        {
#Distributed Dipole for O, H1, H2
#Charge
("O1","charge") : [-0.651],
("H2","charge") : [ 0.325],
("H3","charge") : [ 0.325],
             
#Dipx     Dipy       Dipz
("O1","dipole") :[-0.000,     0.000,     0.342],
("H2", "dipole"):[-0.167,     0.000,    -0.129],
("H3", "dipole"):[ 0.167,     0.000,    -0.129],
               
#Quadxx    Quadxy    Quadxz    Quadyy     Quadyz    Quadzz
("O1","quadrupole") :[-3.213,   -0.000,    -0.000,    -4.362,     0.000,    -3.822],
("H2", "quadrupole"):[-0.115,    -0.000,     0.291,    -0.446,    -0.000,    -0.216],
("H3", "quadrupole"):[-0.115,     0.000,    -0.291,    -0.446,    -0.000,    -0.216],

#Alphaxx   Alphaxy   Alphaxz   Alphayy   Alphayz   Alphazz
("O1","alpha") :[3.268,     0.000,      0.000,   1.531,     -0.000,     2.611],
("H2", "alpha"):[1.971,    -0.000,      1.131,   0.752,     -0.000,     1.308],
("H3", "alpha"):[1.971,     0.000,     -1.131,   0.752,     -0.000,     1.308],
#xxx    xxy    xxz    axyy     taxyz   xzz   yyy     Betayyz    Betayzz    Betazzz
("O1", "beta"):[-0.000,     0.000,    -8.567,     0.000,    -0.000,     0.000,     0.000,    -0.373,     0.000,    -4.010],
("H2", "beta"):[-7.465,     0.000,    -4.959,    -1.692,     0.000,    -4.101,     0.000,    -0.982,     0.000,    -3.583],
("H3", "beta"):[ 7.465,     0.000,    -4.959,     1.692,     0.000,     4.101,     0.000,    -0.982,     0.000,    -3.583]
        }

# -------------------------------------------
#
#      TIP3P water model static, 1907, 1064 and 589 nm, both full and LoProp
#
# -------------------------------------------

        self[("TIP3P", "HF", "PVDZ", False, "0.0")] = \
        {
#Charge
("O1","charge"):[0.000],
("H2","charge") :[ 0.],
("H3","charge") :[ 0.],

#Dipx     Dipy       Dipz
("O1","dipole"):[0.000, 0.000, 0.80897], 
("H2", "dipole"):[0., 0., 0.],
("H3", "dipole"):[0., 0., 0.],

#Quadxx    Quadxy    Quadxz    Quadyy     Quadyz    Quadzz
("O1","quadrupole"):[ -3.443 , 0.000, 0.000, -5.254, 0.00 , -4.254],
("H2", "quadrupole"):[0., 0., 0., 0., 0. ,0.],
("H3", "quadrupole"):[0., 0., 0., 0., 0. ,0.],

#Alphaxx   Alphaxy   Alphaxz   Alphayy   Alphayz   Alphazz
("O1","alpha"):[6.90654, 0.00000, -0.00000, 3.04034,  -0.00000, 5.08449],
("H2", "alpha"):[0., 0., 0., 0., 0. ,0.],
("H3", "alpha"):[0., 0., 0., 0., 0. ,0.],
          #Betaxxx    Betaxxy    Betaxxz    Betaxyy     Betaxyz   Betaxzz   Betayyy     Betayyz    Betayzz    Betazzz
("O1","beta"):[0.00000, 0.00000, -17.14425, -0.00000, -0.00000, -0.00000, 0.00000, -2.33892, 0.00000, -10.64030],
("H2", "beta"):[0., 0., 0., 0., 0. ,0., 0., 0., 0., 0.],
("H3", "beta"):[0., 0., 0., 0., 0. ,0., 0., 0., 0., 0.],
        }

        self[("TIP3P", "HF", "PVDZ", True, "0.0")] = \
        {
("O1","charge") : [-0.659],
("H2","charge") : [ 0.329],
("H3","charge") : [ 0.329],

("O1", "dipole") :  [ -0.000,    -0.000,     0.339],
("H2", "dipole") :  [ -0.170,    -0.000,    -0.130],
("H3", "dipole") :  [ 0.170,      0.000,    -0.130 ],

#Quadxx    Quadxy    Quadxz    Quadyy     Quadyz    Quadzz
("O1", "quadrupole"):[ -3.209, 0.000, 0.000, -4.352,-0.000,-3.820],
("H2", "quadrupole"):[ -0.114, 0.000, 0.284, -0.438, 0.000,-0.216],
("H3", "quadrupole"):[ -0.114, 0.000,-0.284, -0.438,-0.000,-0.216],

#Alphaxx   Alphaxy   Alphaxz   Alphayy   Alphayz   Alphazz
("O1", "alpha"):[ 3.134, 0.000, 0.000, 1.537, -0.000, 2.543,],
("H2", "alpha"):[ 1.886, 0.000, 1.080, 0.752,  0.000, 1.271,],
("H3", "alpha"):[ 1.886, 0.000,-1.080, 0.752, -0.000, 1.271,],

          #Betaxxx    Betaxxy    Betaxxz    Betaxyy     Betaxyz   Betaxzz   Betayyy     Betayyz    Betayzz    Betazzz
("O1", "beta"):[-0.000, 0.000, -7.893,  -0.000,  -0.000,  -0.000,  -0.000,  -0.387,   0.000, -3.752,],
("H2", "beta"):[-6.968, -0.000, -4.627,  -1.676,  -0.000,  -3.876,  -0.000,  -0.976,  -0.000, -3.446,],
("H3", "beta"):[ 6.968,  0.000, -4.627,   1.676,  -0.000,   3.876,  -0.000,  -0.976,   0.000, -3.446]
        }
#---------------------------------
#
# Frequency dependent properties
#
#---------------------------------

# 1907 nm

        self[("TIP3P", "HF", "PVDZ", False, "0.0238927")] = \
        {
("O1","charge") : [0.000],
("H2","charge") : [0.000],
("H3","charge") : [0.000],

("O1","dipole"):    [  0.000,     0.000,     0.80897 ], 
("H2", "dipole") :  [ -0.000,    -0.000,    -0.000   ],
("H3", "dipole") :  [  0.000,     0.000,    -0.000   ],

("O1", "alpha"):[ 6.91536,  -0.000,    0.000,  3.04475,  -0.000, 5.09129,],
("H2", "alpha"):[ 0.000,    -0.000,    0.000,  0.000,    -0.000,     0.000,],
("H3", "alpha"):[ 0.000,     0.000,   -0.000,  0.000,     0.000,     0.000,],
("O1", "beta"):[0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
("H2", "beta"):[0., 0., 0., 0., 0. ,0., 0., 0., 0., 0.],
("H3", "beta"):[0., 0., 0., 0., 0. ,0., 0., 0., 0., 0.],
        }

        self[("TIP3P", "HF", "PVDZ", True, "0.0238927")] = \
        {
("O1","charge") : [-0.659],
("H2","charge") : [ 0.329],
("H3","charge") : [ 0.329],

("O1", "dipole") :  [ -0.000,    -0.000,     0.339],
("H2", "dipole") :  [ -0.170,    -0.000,    -0.130],
("H3", "dipole") :  [  0.170,     0.000,    -0.130 ],

("O1", "alpha"):[ 3.124,    -0.000,    0.000,  1.483,    -0.000,     2.504,],
("H2", "alpha"):[ 1.896,    -0.000,    1.104,  0.781,    -0.000,     1.293,],
("H3", "alpha"):[ 1.896,     0.000,   -1.104,  0.781,     0.000,     1.293,],
("O1", "beta"):[0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
("H2", "beta"):[0., 0., 0., 0., 0. ,0., 0., 0., 0., 0.],
("H3", "beta"):[0., 0., 0., 0., 0. ,0., 0., 0., 0., 0.],
        }

        self[("TIP3P", "HF", "PVDZ", False, "0.0428227")] = \
        {
("O1","charge") : [ 0.000],
("H2","charge") : [ 0.000],
("H3","charge") : [ 0.000],

("O1", "dipole") :  [ -0.000,    -0.000,     0.80897],
("H2", "dipole") :  [ -0.000,    -0.000,    -0.000],
("H3", "dipole") :  [  0.000,     0.000,    -0.000 ],

("O1", "alpha"):[ 6.93497,  -0.000,    0.000,  3.05464, -0.000,     5.10645,],
("H2", "alpha"):[ 0.000,    -0.000,    0.000,  0.000,    -0.000,     0.000,],
("H3", "alpha"):[ 0.000,     0.000,   -0.000,  0.000,     0.000,     0.000,],
("O1", "beta"):[0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
("H2", "beta"):[0., 0., 0., 0., 0. ,0., 0., 0., 0., 0.],
("H3", "beta"):[0., 0., 0., 0., 0. ,0., 0., 0., 0., 0.],

        }

        self[("TIP3P", "HF", "PVDZ", True, "0.0428227")] = \
        {
("O1","charge") : [-0.659],
("H2","charge") : [ 0.329],
("H3","charge") : [ 0.329],

("O1", "dipole") :  [ -0.000,    -0.000,     0.339],
("H2", "dipole") :  [ -0.170,    -0.000,    -0.130],
("H3", "dipole") :  [  0.170,     0.000,    -0.130 ],

("O1", "alpha"):[ 3.120,    -0.000,    0.000,  1.439,    -0.000,     2.477,],
("H2", "alpha"):[ 1.907,    -0.000,    1.104,  0.808,    -0.000,     1.315,],
("H3", "alpha"):[ 1.907,     0.000,   -1.104,  0.808,     0.000,     1.315,],
("O1", "beta"):[0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
("H2", "beta"):[0., 0., 0., 0., 0. ,0., 0., 0., 0., 0.],
("H3", "beta"):[0., 0., 0., 0., 0. ,0., 0., 0., 0., 0.],

        }
        self[("TIP3P", "HF", "PVDZ", False, "0.0773571")] = \
        {
("O1","charge") : [0.0],
("H2","charge") : [0.0],
("H3","charge") : [0.0],

("O1", "dipole") :  [ -0.000,  -0.000,   0.80897],
("H2", "dipole") :  [ 0.000, 0.000, 0.000, ],
("H3", "dipole") :  [ 0.000, 0.000, 0.000, ],

("O1", "alpha"): [7.00055, -0.00000,  0.00000,  3.08864, -0.00000,  5.15748,],
("H2", "alpha"):[ 0., 0., 0., 0., 0., 0., ], 
("H3", "alpha"):[ 0., 0., 0., 0., 0., 0., ],
("O1", "beta"):[0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
("H2", "beta"):[0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
("H3", "beta"):[0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
        }

        self[("TIP3P", "HF", "PVDZ", True, "0.0773571")] = \
        {
("O1","charge") : [-0.659],
("H2","charge") : [ 0.329],
("H3","charge") : [ 0.329],

("O1", "dipole") :  [ -0.000,    -0.000,     0.339],
("H2", "dipole") :  [ -0.170,    -0.000,    -0.130],
("H3", "dipole") :  [ 0.170,      0.000,    -0.130 ],

("O1", "alpha"):[ 3.125,    -0.000,      0.000,    1.357,    -0.000,     2.435, ],
("H2", "alpha"):[ 1.938,    -0.000,      1.134,    0.866,    -0.000,     1.361, ],
("H3", "alpha"):[ 1.938,     0.000,     -1.134,    0.866,     0.000,     1.361, ],
("O1", "beta"):[0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
("H2", "beta"):[0., 0., 0., 0., 0. ,0., 0., 0., 0., 0.],
("H3", "beta"):[0., 0., 0., 0., 0. ,0., 0., 0., 0., 0.],
        }

###########
#       Ano-1 basis set 2s1p contrction on H / 3s2p1d contr on Oxygen
##########
        self[("TIP3P", "HF", "ANOPVDZ", True, "0.0")] = \
        {
("O1","charge"): [-0.678,],
("H2","charge"): [0.339, ],
("H3","charge"): [0.339, ],

("O1", "dipole"): [  0.000,    -0.000,     0.298, ],
("H2", "dipole"): [ -0.154,     0.000,    -0.131, ],
("H3", "dipole"): [  0.154,     0.000,    -0.131, ],

("O1", "alpha"): [ 3.873, -0.000,  -0.000,   3.343,  -0.000,  3.773, ],
("H2", "alpha"): [ 1.819, -0.000,   1.069,   1.198,  -0.000,  1.544, ],
("H3", "alpha"): [ 1.819, -0.000,  -1.069,   1.198,   0.000,  1.544, ],

("O1", "beta"):[0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
("H2", "beta"):[0., 0., 0., 0., 0. ,0., 0., 0., 0., 0.],
("H3", "beta"):[0., 0., 0., 0., 0. ,0., 0., 0., 0., 0.],
}

        self[("TIP3P", "HF", "ANOPVDZ", False, "0.0")] = \
                {
("O1","charge"): [0.],
("H2","charge"): [0.],
("H3","charge"): [0.],

("O1", "dipole"): [ 0.00000,  0.00000,  0.78706,   ],
("H2", "dipole"): [ 0.0,   0.0,   0.0,    ],
("H3", "dipole"): [ 0.0,   0.0,   0.0,    ],

("O1", "alpha"): [ 7.51089, -0.00000,  -0.00000, 5.74008, -0.00000,  6.86089,  ],
("H2", "alpha"): [ 0.00,   0.0, 0.0,  0.0,  0.0,  0.0,  ],
("H3", "alpha"): [ 0.00,   0.0, 0.0,  0.0,  0.0,  0.0,  ],
("O1", "beta"):[0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
("H2", "beta"):[0., 0., 0., 0., 0. ,0., 0., 0., 0., 0.],
("H3", "beta"):[0., 0., 0., 0., 0. ,0., 0., 0., 0., 0.],
                        }
# omega = 1907 nm distributed
#
        self[("TIP3P", "HF", "ANOPVDZ", True, "0.0238927")] = \
{
("O1","charge"):[ -0.678,], 
("H2","charge"):[  0.339, ], 
("H3","charge"):[  0.339, ], 

("O1", "dipole"): [  0.000,    -0.000,     0.298, ],
("H2", "dipole"): [ -0.154,     0.000,    -0.131, ],
("H3", "dipole"): [  0.154,     0.000,    -0.131, ],

("O1", "alpha"): [ 3.843, -0.000, -0.000,  3.245, -0.000,  3.710,],
("H2", "alpha"): [ 1.838, -0.000,  1.091,  1.254, -0.000,  1.580,],
("H3", "alpha"): [ 1.838, -0.000, -1.091,  1.254,  0.000,  1.580,],
("O1", "beta"):[0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
("H2", "beta"):[0., 0., 0., 0., 0. ,0., 0., 0., 0., 0.],
("H3", "beta"):[0., 0., 0., 0., 0. ,0., 0., 0., 0., 0.],
}

        self[("TIP3P", "HF", "ANOPVDZ", False, "0.0238927")] = \
                {
("O1","charge"): [0.],
("H2","charge"): [0.],
("H3","charge"): [0.],

("O1", "dipole"): [ 0.00000,  0.00000,  0.78706,   ],
("H2", "dipole"): [ 0.0,   0.0,   0.0,    ],
("H3", "dipole"): [ 0.0,   0.0,   0.0,    ],

("O1", "alpha"): [7.51884, -0.00000,  -0.00000, 5.75300,  -0.00000,  6.87101,   ],
("H2", "alpha"): [ 0.00,   0.0, 0.0,  0.0,  0.0,  0.0,  ],
("H3", "alpha"): [ 0.00,   0.0, 0.0,  0.0,  0.0,  0.0,  ],

("O1", "beta") : [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
("H2", "beta") : [0., 0., 0., 0., 0. ,0., 0., 0., 0., 0.],
("H3", "beta") : [0., 0., 0., 0., 0. ,0., 0., 0., 0., 0.],
                        }

# At 1064 nm

        self[("TIP3P", "HF", "ANOPVDZ", True, "0.0428227")] = \
{
("O1","charge"):  [-0.678,],
("H2","charge"):  [0.339, ],
("H3","charge"):  [0.339, ],

("O1", "dipole"): [  0.000,    -0.000,     0.298, ], 
("H2", "dipole"): [ -0.154,     0.000,    -0.131, ],
("H3", "dipole"): [  0.154,     0.000,    -0.131, ],

("O1", "alpha"): [ 3.825, -0.000,  -0.000, 3.172, -0.000,  3.667,],
("H2", "alpha"): [ 1.856, -0.000,   1.112, 1.305, -0.000,  1.613,],
("H3", "alpha"): [ 1.856, -0.000,  -1.112, 1.305,  0.000,  1.613,],
("O1", "beta"):[0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
("H2", "beta"):[0., 0., 0., 0., 0. ,0., 0., 0., 0., 0.],
("H3", "beta"):[0., 0., 0., 0., 0. ,0., 0., 0., 0., 0.],

        }

        self[("TIP3P", "HF", "ANOPVDZ", False, "0.0428227")] = \
                {

("O1","charge"): [0.],
("H2","charge"): [0.],
("H3","charge"): [0.],

("O1", "dipole"): [ 0.00000,  0.00000,  0.78706,   ],
("H2", "dipole"): [ 0.0,   0.0,   0.0,    ],
("H3", "dipole"): [ 0.0,   0.0,   0.0,    ],

("O1", "alpha"): [ 7.53651, -0.00000,  -0.00000, 5.78202,  -0.00000,  6.89360,   ],
("H2", "alpha"): [ 0.00,   0.0, 0.0,  0.0,  0.0,  0.0,  ],
("H3", "alpha"): [ 0.00,   0.0, 0.0,  0.0,  0.0,  0.0,  ],
("O1", "beta"):[0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
("H2", "beta"):[0., 0., 0., 0., 0. ,0., 0., 0., 0., 0.],
("H3", "beta"):[0., 0., 0., 0., 0. ,0., 0., 0., 0., 0.],
                }

# At 589 nm
        self[("TIP3P", "HF", "ANOPVDZ", True, "0.0773571")] = \
                {
("O1","charge"): [-0.678,],
("H2","charge"): [0.339, ],
("H3","charge"): [0.339, ],

("O1", "dipole"): [  0.000,    -0.000,     0.298, ], 
("H2", "dipole"): [ -0.154,     0.000,    -0.131, ],
("H3", "dipole"): [  0.154,     0.000,    -0.131, ],

("O1", "alpha"): [ 3.801, -0.000,  0.000,  3.045, -0.000,  3.599,],
("H2", "alpha"): [ 1.897, -0.000,  1.157,  1.419, -0.000,  1.686,],
("H3", "alpha"): [ 1.897, -0.000, -1.157,  1.419,  0.000,  1.686,],
("O1", "beta"):[0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
("H2", "beta"):[0., 0., 0., 0., 0. ,0., 0., 0., 0., 0.],
("H3", "beta"):[0., 0., 0., 0., 0. ,0., 0., 0., 0., 0.],
        }

        self[("TIP3P", "HF", "ANOPVDZ", False, "0.0773571")] = \
                {

("O1","charge"): [0.],
("H2","charge"): [0.],
("H3","charge"): [0.],

("O1", "dipole"): [ 0.00000,  0.00000,  0.78706,   ],
("H2", "dipole"): [ 0.0,   0.0,   0.0,    ],
("H3", "dipole"): [ 0.0,   0.0,   0.0,    ],

("O1", "alpha"): [ 7.59551, -0.00000,  0.00000, 5.88201, -0.00000,  6.96982,],
("H2", "alpha"): [ 0.00,   0.0, 0.0,  0.0,  0.0,  0.0,  ],
("H3", "alpha"): [ 0.00,   0.0, 0.0,  0.0,  0.0,  0.0,  ],
("O1", "beta"):[0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
("H2", "beta"):[0., 0., 0., 0., 0. ,0., 0., 0., 0., 0.],
("H3", "beta"):[0., 0., 0., 0., 0. ,0., 0., 0., 0., 0.],

}

###########
#       Ano-1 basis set 3s2p1d contr on H / 4s3p2d1f contr on O
##########


        self[("TIP3P", "HF", "ANOPVTZ", True, "0.0")] = \
        {
("O1","charge"): [-0.689,],
("H2","charge"): [0.345, ],
("H3","charge"): [0.345, ],

("O1", "dipole"): [ -0.000,0.000,0.269,  ], 
("H2", "dipole"): [ -0.149,-0.000,-0.128,], 
("H3", "dipole"): [ 0.149,0.000,-0.128,  ], 

("O1", "alpha"): [ 5.000,0.000,-0.000,5.425,-0.000,5.214],
("H2", "alpha"): [2.067,-0.000,1.110, 1.188,-0.000,1.614],
("H3", "alpha"): [2.067,-0.000,-1.110,1.188,-0.000,1.614],
("O1", "beta"):[0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
("H2", "beta"):[0., 0., 0., 0., 0. ,0., 0., 0., 0., 0.],
("H3", "beta"):[0., 0., 0., 0., 0. ,0., 0., 0., 0., 0.],
}

        self[("TIP3P", "HF", "ANOPVTZ", False, "0.0")] = \
                {
("O1","charge"): [0.],
("H2","charge"): [0.],
("H3","charge"): [0.],

("O1", "dipole"): [ -0.00000,  0.00000,  0.77468 ],
("H2", "dipole"): [ 0.0,   0.0,   0.0,    ],
("H3", "dipole"): [ 0.0,   0.0,   0.0,    ],

("O1", "alpha"): [ 9.13307,  0.00000,  -0.00000,  7.80083,  -0.00000,  8.44263, ],
("H2", "alpha"): [ 0.00,   0.0, 0.0,  0.0,  0.0,  0.0,  ],
("H3", "alpha"): [ 0.00,   0.0, 0.0,  0.0,  0.0,  0.0,  ],
("O1", "beta"):[0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
("H2", "beta"):[0., 0., 0., 0., 0. ,0., 0., 0., 0., 0.],
("H3", "beta"):[0., 0., 0., 0., 0. ,0., 0., 0., 0., 0.],
                        }

        self[("TIP3P", "HF", "ANOPVTZ", True, "0.0238927")] = \
{
("O1","charge"): [-0.689],
("H2","charge"): [0.345], 
("H3","charge"): [0.345], 

("O1", "dipole"):   [-0.000,0.000,0.269], 
("H2", "dipole"):  [-0.149,-0.000,-0.128],
("H3", "dipole"):  [0.149,0.000,-0.128],  

("O1", "alpha"):[4.973,0.000,  -0.000, 5.350,  -0.000,5.163],
("H2", "alpha"):[2.086,-0.000, 1.129,  1.233,  -0.000,1.646],
("H3", "alpha"):[2.086,-0.000, -1.129, 1.233,  -0.000,1.646],
("O1", "beta"):[0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
("H2", "beta"):[0., 0., 0., 0., 0. ,0., 0., 0., 0., 0.],
("H3", "beta"):[0., 0., 0., 0., 0. ,0., 0., 0., 0., 0.],
}

        self[("TIP3P", "HF", "ANOPVTZ", False, "0.0238927")] = \
                {
("O1","charge"): [0.],
("H2","charge"): [0.],
("H3","charge"): [0.],

("O1", "dipole"): [ -0.00000,  0.00000,  0.77468 ],
("H2", "dipole"): [ 0.0,   0.0,   0.0,    ],
("H3", "dipole"): [ 0.0,   0.0,   0.0,    ],

("O1", "alpha"): [ 9.14435,  0.00000,  -0.00000, 7.81558, -0.00000, 8.45521],
("H2", "alpha"): [ 0.00,   0.0, 0.0,  0.0,  0.0,  0.0,  ],
("H3", "alpha"): [ 0.00,   0.0, 0.0,  0.0,  0.0,  0.0,  ],
("O1", "beta"):[0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
("H2", "beta"):[0., 0., 0., 0., 0. ,0., 0., 0., 0., 0.],
("H3", "beta"):[0., 0., 0., 0., 0. ,0., 0., 0., 0., 0.],
                        }


        self[("TIP3P", "HF", "ANOPVTZ", True, "0.0428227")] = \
{
("O1","charge"):[-0.689],
("H2","charge"):[0.345], 
("H3","charge"):[0.345], 

("O1", "dipole"):  [-0.000,0.000,0.269], 
("H2", "dipole"): [-0.149,-0.000,-0.128],
("H3", "dipole"): [0.149,0.000,-0.128],  

("O1", "alpha"): [4.960,  0.000, -0.000, 5.301, -0.000, 5.131],
("H2", "alpha"): [2.105, -0.000,  1.146, 1.274, -0.000, 1.676],
("H3", "alpha"): [2.105, -0.000, -1.146, 1.274, -0.000, 1.676],
("O1", "beta"):[0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
("H2", "beta"):[0., 0., 0., 0., 0. ,0., 0., 0., 0., 0.],
("H3", "beta"):[0., 0., 0., 0., 0. ,0., 0., 0., 0., 0.],

        }

        self[("TIP3P", "HF", "ANOPVTZ", False, "0.0428227")] = \
                {

("O1","charge"): [0.],
("H2","charge"): [0.],
("H3","charge"): [0.],

("O1", "dipole"): [ -0.00000,  0.00000,  0.77468 ],
("H2", "dipole"): [ 0.0,   0.0,   0.0,    ],
("H3", "dipole"): [ 0.0,   0.0,   0.0,    ],

("O1", "alpha"): [ 9.16940,  0.00000,  -0.00000,7.84863,  -0.00000,  8.48327, ],
("H2", "alpha"): [ 0.00,   0.0, 0.0,  0.0,  0.0,  0.0,  ],
("H3", "alpha"): [ 0.00,   0.0, 0.0,  0.0,  0.0,  0.0,  ],

("O1", "beta"):[0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
("H2", "beta"):[0., 0., 0., 0., 0. ,0., 0., 0., 0., 0.],
("H3", "beta"):[0., 0., 0., 0., 0. ,0., 0., 0., 0., 0.],
                }
        self[("TIP3P", "HF", "ANOPVTZ", True, "0.0773571")] = \
                {
("O1","charge"): [-0.689],
("H2","charge"): [0.345], 
("H3","charge"): [0.345], 

("O1", "dipole"):   [-0.000,0.000,0.269], 
("H2", "dipole"):  [-0.149,-0.000,-0.128],
("H3", "dipole"):  [0.149,0.000,-0.128],  

("O1", "alpha"):[4.951,0.000,  -0.000,5.229, -0.000,5.090],
("H2", "alpha"):[2.151,-0.000,1.187, 1.366, -0.000,1.744],
("H3", "alpha"):[2.151,-0.000,  -1.187,1.366, -0.000,1.744],

("O1", "beta"):[0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
("H2", "beta"):[0., 0., 0., 0., 0. ,0., 0., 0., 0., 0.],
("H3", "beta"):[0., 0., 0., 0., 0. ,0., 0., 0., 0., 0.],
        }

        self[("TIP3P", "HF", "ANOPVTZ", False, "0.0773571")] = \
                {

("O1","charge"): [0.],
("H2","charge"): [0.],
("H3","charge"): [0.],

("O1", "dipole"): [ -0.00000,  0.00000,  0.77468 ],
("H2", "dipole"): [ 0.0,   0.0,   0.0,    ],
("H3", "dipole"): [ 0.0,   0.0,   0.0,    ],

("O1", "alpha"): [ 9.25300,  0.00000,  -0.00000,7.96160,  -0.00000, 8.57770 ],
("H2", "alpha"): [ 0.00,   0.0, 0.0,  0.0,  0.0,  0.0,  ],
("H3", "alpha"): [ 0.00,   0.0, 0.0,  0.0,  0.0,  0.0,  ],

("O1", "beta"):[0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
("H2", "beta"):[0., 0., 0., 0., 0. ,0., 0., 0., 0., 0.],
("H3", "beta"):[0., 0., 0., 0., 0. ,0., 0., 0., 0., 0.],

}



    def get(self, model = "OLAV", method = "HF", basis ="PVDZ",
            dist = False, freq = "0.0" ):
        return self[ ( model, method, basis, dist, freq) ]

if __name__ == '__main__':
    #Perform tests om Tempaltes
    pass
   
