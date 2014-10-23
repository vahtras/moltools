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

        self[("TIP3P", "HF", "PVDZ", False, "0.0773571")] = \
        {
("O1","charge") : [0.0],
("H2","charge") : [0.0],
("H3","charge") : [0.0],

("O1", "dipole") :  [ -0.000,    -0.000,   0.80897],
("H2", "dipole") :  [ 0.000, 0.000, 0.000, ],
("H3", "dipole") :  [ 0.000, 0.000, 0.000, ],

##Quadxx    Quadxy    Quadxz    Quadyy     Quadyz    Quadzz
#("O1", "quadrupole"):[ -3.209, 0.000, 0.000, -4.352,-0.000,-3.820],
#("H2", "quadrupole"):[ -0.114, 0.000, 0.284, -0.438, 0.000,-0.216],
#("H3", "quadrupole"):[ -0.114, 0.000,-0.284, -0.438,-0.000,-0.216],

#Alphaxx   Alphaxy   Alphaxz   Alphayy   Alphayz   Alphazz
("O1", "alpha"): [7.00055, -0.00000,  3.08864,  0.00000, -0.00000,  5.15748,],
("H2", "alpha"):[ 0., 0., 0., 0., 0., 0., ], 
("H3", "alpha"):[ 0., 0., 0., 0., 0., 0., ]

          #Betaxxx    Betaxxy    Betaxxz    Betaxyy     Betaxyz   Betaxzz   Betayyy     Betayyz    Betayzz    Betazzz
#("O1", "beta"):[-0.000, 0.000, -7.893,  -0.000,  -0.000,  -0.000,  -0.000,  -0.387,   0.000, -3.752,],
#("H2", "beta"):[-6.968, -0.000, -4.627,  -1.676,  -0.000,  -3.876,  -0.000,  -0.976,  -0.000, -3.446,],
#("H3", "beta"):[ 6.968,  0.000, -4.627,   1.676,  -0.000,   3.876,  -0.000,  -0.976,   0.000, -3.446]
        }

        self[("TIP3P", "HF", "PVDZ", True, "0.0773571")] = \
        {
("O1","charge") : [-0.659],
("H2","charge") : [ 0.329],
("H3","charge") : [ 0.329],

("O1", "dipole") :  [ -0.000,    -0.000,     0.339],
("H2", "dipole") :  [ -0.170,    -0.000,    -0.130],
("H3", "dipole") :  [ 0.170,      0.000,    -0.130 ],

##Quadxx    Quadxy    Quadxz    Quadyy     Quadyz    Quadzz
#("O1", "quadrupole"):[ -3.209, 0.000, 0.000, -4.352,-0.000,-3.820],
#("H2", "quadrupole"):[ -0.114, 0.000, 0.284, -0.438, 0.000,-0.216],
#("H3", "quadrupole"):[ -0.114, 0.000,-0.284, -0.438,-0.000,-0.216],

#Alphaxx   Alphaxy   Alphaxz   Alphayy   Alphayz   Alphazz
("O1", "alpha"):[ 3.125,    -0.000,      0.000,    1.357,    -0.000,     2.435, ],
("H2", "alpha"):[ 1.938,    -0.000,      1.134,    0.866,    -0.000,     1.361, ],
("H3", "alpha"):[ 1.938,     0.000,     -1.134,    0.866,     0.000,     1.361, ],

          #Betaxxx    Betaxxy    Betaxxz    Betaxyy     Betaxyz   Betaxzz   Betayyy     Betayyz    Betayzz    Betazzz
#("O1", "beta"):[-0.000, 0.000, -7.893,  -0.000,  -0.000,  -0.000,  -0.000,  -0.387,   0.000, -3.752,],
#("H2", "beta"):[-6.968, -0.000, -4.627,  -1.676,  -0.000,  -3.876,  -0.000,  -0.976,  -0.000, -3.446,],
#("H3", "beta"):[ 6.968,  0.000, -4.627,   1.676,  -0.000,   3.876,  -0.000,  -0.976,   0.000, -3.446]
        }


    def get(self, model = "OLAV", method = "HF", basis ="PVDZ",
            dist = False, freq = "0.0" ):
        return self[ ( model, method, basis, dist, freq) ]

if __name__ == '__main__':
    #Perform tests om Tempaltes
    pass
   
