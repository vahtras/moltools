#!/usr/bin/env python
#-*- coding: utf-8 -*-

class Template:
    def __init__(self):


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


    #Dictionaries from basis
        olav_hf_dict = { "PVDZ" : olav_hf_cc_pVDZ }
    #Dictionaries from method
        olav_method_dict = { "HF" : olav_hf_dict }

        self.nameDict = { "OLAV" : olav_method_dict }


    def getData(self, model, method, basis):
        return self.nameDict[model][method][basis]

    def get_dist_data(self, model = "OLAV", method = "HF", basis ="PVDZ"):
        return self.olav_hf_cc_pVDZ_dist

    def get_data(self, model = "OLAV", method = "HF", basis ="PVDZ"):
        return self.olav_hf_cc_pVDZ

if __name__ == '__main__':
    #Perform tests om Tempaltes
    pass
   
