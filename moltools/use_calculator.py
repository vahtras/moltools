#!/usr/bin/env python
#-*- coding: utf-8 -*-

__all__ = [ 'Calculator' ]

"""
Purpose
------------------

Read in output from dalton .out files, produced by the .mol files created with :ref:`use_generator`.


Usage
------------------

Run the script in the same directory where all the .out and .mol files exists from the dalton \*\*QUADRATIC response calculation.

Typical dalton file will look like:

::

    **DALTON INPUT
    .RUN RESPONSE
    .DIRECT
    .PARALLELL
    **WAVE FUNCTION
    .HF
    **RESPONSE
    .PROPAV
    XDIPLEN
    .PROPAV
    YDIPLEN
    .PROPAV
    ZDIPLEN
    *QUADRATIC
    .DIPLEN
    **END OF DALTON INPUT
    
.. note::

    Right now the dalton input must be named hfqua.dal, which will give a hfqua\_ prefix for all the .out files.

-----------------


.. code:: bash

    $ use_calculator.py -vary_r 

Will generate tmp.xvg, which can be plotted as

.. code:: bash

    $ xmgrace tmp.xvg

"""

import os, re, subprocess, argparse

import numpy as np
import math as m


from .template import Template
from .molecules import *

from .read_dal import *

try:
    from applequist.gaussian import *
except ImportError:
    pass

a0 = 0.52917721092
charge_dic = {"H": 1.0, "C": 6.0, "N": 7.0, "O": 8.0, "S": 16.0}
mass_dict = {"H": 1.008,  "C": 6.0, "N": 7.0, "O": 15.999, "S": 16.0}
index_dict = { 
        ('zeroth', 'dipole', 'X'): (0, 0,),
        ('zeroth', 'dipole', 'Y'): (0, 1,),
        ('zeroth', 'dipole', 'Z'): (0, 2,),
        ('linear', 'dipole', 'X'): (1, 0,),
        ('linear', 'dipole', 'Y'): (1, 1,),
        ('linear', 'dipole', 'Z'): (1, 2,),
        ('linear', 'alpha', 'X'): (3, 0,),
        ('linear', 'alpha', 'Y'): (3, 1,),
        ('linear', 'alpha', 'Z'): (3, 2,),
        ('quadratic', 'dipole', 'X'): (2, 0,),
        ('quadratic', 'dipole', 'Y'): (2, 1,),
        ('quadratic', 'dipole', 'Z'): (2, 2,),
        ('quadratic', 'alpha', 'X'): (4, 0,),
        ('quadratic', 'alpha', 'Y'): (4, 1,),
        ('quadratic', 'alpha', 'Z'): (4, 2,),
        ('quadratic', 'beta', 'X'): (5, 0,),
        ('quadratic', 'beta', 'Y'): (5, 1,),
        ('quadratic', 'beta', 'Z'): (5, 2,),
        }

line_thick_dict = { 
        "zeroth" : {"dipole":{"X": 0, "Y": 0, "Z": 0}},

        "linear" : {"dipole":{"X": 2 , "Y": 2, "Z": 2},
                   "alpha":{"X":  2 , "Y": 2, "Z": 2} },
        "quadratic" : { "dipole":{"X": 3 , "Y": 3, "Z": 3 },
                   "alpha":{  "X": 3 , "Y": 3, "Z": 3 },
                   "beta":{   "X": 3 , "Y": 3, "Z": 3 }}}
line_style_dict = { 
        "zeroth" : {"dipole":{"X": 1, "Y": 1, "Z": 1}},

        "linear" : {"dipole":{"X": 2 , "Y": 2, "Z": 2},
                   "alpha":{"X":  2 , "Y": 2, "Z": 2} },
        "quadratic" : { "dipole":{"X": 3, "Y": 3, "Z": 3 },
                   "alpha":{  "X": 3, "Y": 3, "Z": 3 },
                   "beta":{   "X": 3, "Y": 3, "Z": 3 }}}

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

class Calculator( dict ):
    """ Class container for retrieving and calculating data
    retrieving water classes from DALTON .mol / .pdb / .xyz
    retrieving water classes  from .pdb / .xyz
    retrieves properties from DALTON .out
    """

    def __init__(self):

        self.opts = { 
             "r"      : {"constant" : "5.00"  },
             "tau"  : {"constant" : "0.00" },
             "theta"    : {"constant" : "0.00" },
             "rho1"   : {"constant" : "0.00" },
             "rho2"   : {"constant" : "0.00" },
             "rho3"   : {"constant" : "0.00" }
             }
        self.zeroth = {}
        self.linear = {}
        self.quadratic = {}

    def get_x_and_y( self, 
            level = 'linear',
            var = 'r',
            in_AA = False,
            out_AA = False,
            max_l = 1,
            dist = 0, 
            Rq = "0.00",
            Rp = "0.00",
            out_degree = True,
            ):

        x = []
        y = []
        for i in self.get_matching_out_and_mol():
            r, tau, theta, rho1, rho2, rho3 = i.split('-')
            if self.continue_if_not_relevant_var( r, tau, theta,
                    rho1, rho2, rho3):
                continue

            if var == 'r':
#To be able to plot angstrom if water was calculated in AU
                r_x = "%.2f" % float(r)
                if out_AA and not in_AA:
                    r_x = "%.2f" % (float(r)*a0)
                x.append( r_x )
                y.append( self[ ( r, tau, theta, rho1, rho2, rho3, max_l, dist, Rq, Rp ) ] )
                continue

            if var == 'tau':
#To be able to plot angstrom if water was calculated in AU
                r_x = float(tau)
                if out_degree:
                    r_x *= 180/np.pi
                x.append( r_x )
                y.append( self[ ( r, tau, theta, rho1, rho2, rho3, max_l, dist, Rq, Rp ) ] )

            if var == 'theta':
#To be able to plot angstrom if water was calculated in AU
                r_x = float(theta)
                if out_degree:
                    r_x *= 180/np.pi
                x.append( r_x )
                y.append( self[ ( r, tau, theta, rho1, rho2, rho3, max_l, dist, Rq, Rp ) ] )
            if var == 'rho1':
#To be able to plot angstrom if water was calculated in AU
                r_x = float(rho1)
                if out_degree:
                    r_x *= 180/np.pi
                x.append( r_x )
                y.append( self[ ( r, tau, theta, rho1, rho2, rho3, max_l, dist, Rq, Rp ) ] )
            if var == 'rho2':
#To be able to plot angstrom if water was calculated in AU
                r_x = float(rho2)
                if out_degree:
                    r_x *= 180/np.pi
                x.append( r_x )
                y.append( self[ ( r, tau, theta, rho1, rho2, rho3, max_l, dist, Rq, Rp ) ] )
            if var == 'rho3':
#To be able to plot angstrom if water was calculated in AU
                r_x = float(rho3)
                if out_degree:
                    r_x *= 180/np.pi
                x.append( r_x )
                y.append( self[ ( r, tau, theta, rho1, rho2, rho3, max_l, dist, Rq, Rp ) ] )
        return  x , y 

    def choose_relevant_var( self,r, tau, theta, rho1, rho2, rho3):
        if self.opts["r"].has_key( "vary" ):
            return r
        if self.opts["tau"].has_key( "vary" ):
            return tau
        if self.opts["theta"].has_key( "vary" ):
            return theta
        if self.opts["rho1"].has_key( "vary" ):
            return rho1
        if self.opts["rho2"].has_key( "vary" ):
            return rho2
        if self.opts["rho3"].has_key( "vary" ):
            return rho3
    def two_mols( self, 
            levels = ['zeroth'],
            props = ['dipole'],
            comps = ['Z'],
            mol_model = 'tip3p',
            basis = "anopvdz",
            rel = True,
            noqm = False,
            h5 = False,
            model = "gaussian",
            qm_method = "hfqua",
            max_ls = [1],
            dists = [0],
            in_AA = False,
            out_AA = False,
            Rps = [0.000001],
            Rqs = [0.000001],
            worst = False,
            out = 'data.h5',
            ):
        """combine get_rel_error and get_abs_value"""
#
# XXZ (2), YYZ(7) and ZZZ(9) in ut form
        select = [ (0, 0, 2), (1, 1, 2), (2, 2, 2)]
        for i in self.get_matching_out_and_mol():
            r, tau, theta, rho1, rho2, rho3 = i.split('-')
            for max_l in max_ls:
                for dist in dists:
                    for Rp in Rps:
                        for Rq in Rqs:
# Gather the water molecules from the .mol file
                            Rp_ind = "%.2f" %float(Rp)
                            Rq_ind = "%.2f" %float(Rq)
                            tmp_waters = []
                            for j in Water.read_waters( i + ".mol", in_AA = in_AA,
                                    out_AA = out_AA):

                                templ = Template().get(*( mol_model.upper(), "HF",basis.upper(), dist==1,"0.0"))

#If we want to analyze the properties from the worst-case obained from MD:
#We create new water since h1 and h2 might now be propertly defined for this read
#water due to C1 symmetry
                                if j.is_worst() and worst:
                                    w1 = Water.get_standard( AA = in_AA )
                                    w1.populate_bonds()
                                    w1.populate_angles()
                                    w1.h1.scale_angle( 0.988 )
                                    w1.h1.scale_bond( 0.985 )
                                    w1.h2.scale_bond( 1.015 )
                                    templ = Template().get(*( "%s_WORST"
                                        %mol_model.upper(),
                                        "HF",basis.upper(),
                                        dist == 1,"0.0"))
                                    for at in w1:
                                        Property.add_prop_from_template( at, templ )
                                    tmp_waters.append( w1 )
                                    continue
#End of worst case tip3p
                                for at in j:
                                    Property.add_prop_from_template( at, templ )
                                t1, t2, t3 =  j.get_euler()
                                for at in j:
                                    Property.transform_ut_properties( at.Property, t1, t2, t3 )
                                tmp_waters.append( j )

#Gather the QM result and save it in the dictionary for the relative error
                            if rel:
                                file_name = qm_method + "_" + i + ".out"
                                if self.continue_if_not_relevant_var( r, tau, 
                                        theta, rho1, rho2, rho3):
                                    continue

                                if qm_method == "ccsd":
                                    qm_dipole, qm_alpha, qm_beta = self.get_props_ccsd( file_name )
                                else:
                                    qm_dipole = self.get_qm_dipole( file_name )
                                    qm_alpha = self.get_qm_alpha(  file_name )
                                    qm_beta = self.get_qm_beta(   file_name )
#Defaults to this model

                                if model == "pointdipole":
                                    if "zeroth" in levels:
                                        zeroth= PointDipoleList.from_string( Water.get_string_from_waters( tmp_waters , max_l = max_l , pol = 0, hyper = 0, dist = dist ) )
                                        zeroth.solve_scf()
                                    if "linear" in levels:
                                        linear = PointDipoleList.from_string( Water.get_string_from_waters(  tmp_waters , max_l = max_l, pol = 2, hyper = 1, dist = dist ))
                                        linear.solve_scf()
                                    if "quadratic" in levels:
                                        quadratic = PointDipoleList.from_string( Water.get_string_from_waters(  tmp_waters, max_l = max_l, pol = 22, hyper = 1, dist = dist ))
                                        quadratic.solve_scf()
                                if model == "gaussian":
                                    tmp_Rq = float(Rq)
                                    tmp_Rp = float(Rp)
                                    if "zeroth" in levels:
                                        zeroth = GaussianQuadrupoleList.from_string( Water.get_string_from_waters(  tmp_waters, max_l = max_l , pol = 0 , hyper = 0 ))
                                        zeroth.set_damping( tmp_Rq, tmp_Rp )
                                        zeroth.solve_scf()

                                    if "linear" in levels:
                                        linear = GaussianQuadrupoleList.from_string( Water.get_string_from_waters(  tmp_waters,
                                            max_l = max_l , pol = 2 , hyper = 0 ))
                                        linear.set_damping( tmp_Rq, tmp_Rp )
                                        linear.solve_scf()

                                    if "quadratic" in levels:
                                        quadratic = GaussianQuadrupoleList.from_string( Water.get_string_from_waters(  tmp_waters,
                                            max_l = max_l , pol = 22,  hyper = 1 ))
                                        quadratic.set_damping( tmp_Rq, tmp_Rp )
                                        quadratic.solve_scf()
                                d_zeroth = []
                                d_linear = []
                                d_quadratic = []

                                a_linear = []
                                a_quadratic = []

                                b_quadratic = []


                                if "zeroth" in levels:
                                    d_zeroth =  \
                                             [(this-ref)/ref for this, ref in zip(
                                                    zeroth.total_dipole_moment(),
                                                            qm_dipole )] 

                                if "linear" in levels:
                                    d_linear =  \
                                             [(this-ref)/ref for this, ref in zip(
                                                    linear.total_dipole_moment(),
                                                            qm_dipole )] 

                                    a_linear = \
                                             [(this-ref)/ref for this, ref in zip(
                                                    linear.alpha().diagonal(),
                                                            qm_alpha.diagonal() )] 
                                if "quadratic" in levels:
                                    d_quadratic = \
                                             [(this-ref)/ref for this, ref in zip(
                                                    quadratic.total_dipole_moment(),
                                                            qm_dipole )] 
                                    a_quadratic = \
                                             [(this-ref)/ref for this, ref in zip(
                                                    quadratic.alpha().diagonal(),
                                                            qm_alpha.diagonal() )] 
                                    reference = [ qm_beta[ii, jj, kk] for ii, jj, kk in select ]
                                    b_quadratic =  [(this-ref)/ref for this, ref in zip( [ quadratic.beta()[ii, jj, kk] for ii, jj, kk in select], reference  ) ] 

                                val = [d_zeroth, d_linear, d_quadratic, a_linear, a_quadratic, b_quadratic ]

# Dictionary key for damping coeff is 2 decimal strings
                                Rp_ind = "%.2f" %float(Rp)
                                Rq_ind = "%.2f" %float(Rq)
                                self[ (r, tau, theta, rho1, rho2, rho3, max_l, dist, Rq_ind, Rp_ind )] = val
                            else:
# Get absolute values using the quadratic model
                                if noqm:
                                    if model == "pointdipole":
                                        zeroth= PointDipoleList.from_string( self.get_string( tmp_waters ,
                                                max_l = max_l, pol = 0, quadratic = 0, dist = dist ) )
                                        linear = PointDipoleList.from_string( self.get_string(  tmp_waters ,
                                                max_l = max_l, pol = 2, quadratic = 1, dist = dist ))
                                        quadratic = PointDipoleList.from_string( self.get_string(  tmp_waters,
                                                max_l = max_l, pol = 22, quadratic = 1, dist = dist ))

                                    if model == "gaussian":
                                        tmp_Rq = float(Rq)
                                        tmp_Rp = float(Rp)
                                        if "zeroth" in levels:
                                            zeroth = GaussianQuadrupoleList.from_string( Water.get_string_from_waters(  tmp_waters,
                                                max_l = max_l , pol = 0 , hyper = 0 ,  dist = dist ))
                                            zeroth.set_damping( tmp_Rq, tmp_Rp )
                                        if "linear" in levels:
                                            linear = GaussianQuadrupoleList.from_string( Water.get_string_from_waters(  tmp_waters,
                                                max_l = max_l , pol = 2 , hyper = 1 , dist = dist ))
                                            linear.set_damping( tmp_Rq, tmp_Rp )
                                        if "quadratic" in levels:
                                            quadratic = GaussianQuadrupoleList.from_string( Water.get_string_from_waters(  tmp_waters,
                                                max_l = max_l , pol = 22,  hyper = 1 , dist = dist ))
                                            quadratic.set_damping( tmp_Rq, tmp_Rp )
                                    try:
                                        zeroth.solve_scf()
                                    except UnboundLocalError:
                                        pass
                                    try:
                                        linear.solve_scf()
                                    except UnboundLocalError:
                                        pass
                                    try:
                                        quadratic.solve_scf()
                                    except UnboundLocalError:
                                        pass
                                    d_zeroth = []
                                    d_linear = []
                                    d_quadratic = []
                                    a_linear = []
                                    a_quadratic = []
                                    b_quadratic = []
                                    if "zeroth" in levels:
                                        d_zeroth = zeroth.total_dipole_moment( )
                                    if "linear" in levels:
                                        d_linear =  linear.total_dipole_moment( )
                                        a_linear = linear.alpha().diagonal()
                                    if "quadratic" in levels:
                                        d_quadratic = quadratic.total_dipole_moment( )
                                        a_quadratic = quadratic.alpha().diagonal()
                                        b_quadratic =  [ quadratic.beta()[ii, jj, kk] for ii, jj, kk in select]
                                    val = [ d_zeroth, d_linear, d_quadratic, a_linear, a_quadratic, b_quadratic ]
                                    self[ (r, tau, theta, rho1, rho2, rho3, max_l, dist, Rq_ind, Rp_ind )] = val
                                else:
                                    file_name = qm_method + "_" + i + ".out"


#By default read values obtained by QM, will be overridden by -noqm later if want rgs.model value
                                    if qm_method == "ccsd":
                                        qm_dipole, qm_alpha, qm_beta = self.get_props_ccsd( file_name )
                                    else:
                                        qm_dipole = self.get_qm_dipole( file_name )
                                        qm_alpha =  self.get_qm_alpha(  file_name )
                                        qm_beta =   self.get_qm_beta(   file_name )

                                    dipole = qm_dipole
                                    alpha  = np.einsum('ii->i', np.array(qm_alpha))
                                    beta = [ qm_beta[ii,jj,kk] for (ii, jj, kk) in select ]
                                    val = [ dipole, dipole, dipole, alpha, alpha, beta ]
#End of two blocks
                                    self[ (r, tau, theta, rho1, rho2, rho3, max_l, dist, Rq_ind, Rp_ind )] = val

    def get_matching_out_and_mol(self):
        tmp = []
        tmpOut = self.get_out_files()
        tmpFound = {}
        cnt = 0
        for j in self.get_mol_files():
            if os.path.isfile( os.path.join( os.getcwd(),  "hfqua_" + j.rstrip(".mol") + ".out" ) ):
                tmp.append( j.rstrip( ".mol" ))
                continue
            elif os.path.isfile( os.path.join( os.getcwd(),  "ccsd_" + j.rstrip(".mol") + ".out" ) ):
                tmp.append( j.rstrip( ".mol" ))
                continue
            elif os.path.isfile( os.path.join( os.getcwd(),  "dftqua_" + j.rstrip(".mol") + ".out" ) ):
                tmp.append( j.rstrip( ".mol" ))
                continue
            else:
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
        for i in open( fname ).readlines():
            if pat_xyz.match(i):
                f = pat_xyz.match(i).groups()
                matched = pat_xyz.match(i).groups()
                kwargs = { "element" :  matched[0], "x" : matched[1],
                        "y" : matched[2], "z" : matched[3] }
                tmpAtom = Atom( **kwargs )
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
        pat_alpha = re.compile(r'@.*QRLRVE:.*([XYZ])DIPLEN.*([XYZ])DIPLEN')
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

    def get_var(self):
        for i, key in enumerate(self.opts):
            if 'vary' in self.opts[key].keys():
                var = key
                break
        return var

#If data is in absolute value, format xvg string with qm method in labels
    def get_xvg_string( self,
            levels = ["zeroth"], 
            props = ["dipole"], 
            comps = ["Z"], 
            max_ls = [1], 
            rel = True,
            dists = [0],
            noqm = False,
            in_AA = False,
            out_AA = False,
            Rqs = [0.000001],
            Rps = [0.000001],
            qm_method = 'hfqua',
            out_degree = False,
            ):
        """ kwargs is a dictionary where user specifies which components, dipole models and properties
        will be returned and printed for a finished formatted .xvg xmgrace plot"""
        localCounter = 0
        string = ""

        var = self.get_var()

        for max_l in max_ls:
            for dist in dists:
                for Rp in Rps:
                    for Rq in Rqs:
                        x, y = self.get_x_and_y(
                                var = var,
                                in_AA = in_AA,
                                out_AA = out_AA,
                                max_l = max_l,
                                dist = dist,
                                Rq = "%.2f" %float(Rq),
                                Rp = "%.2f" %float(Rp),
                                out_degree = out_degree,
                                )
                        for level in levels:
                            for prop in props:
                                for comp in comps:
                                    try:
                                        ind1, ind2 =  index_dict[ (level, prop, comp) ]
                                        width = line_thick_dict[ level ][ prop ][ comp ]
                                        style = line_style_dict[ level ][ prop ][ comp ]
                                    except KeyError:
                                        print level, prop, comp
                                        print "Skipping (%s, %s, %s, )" %( level, prop, comp )
                                        continue

                                    if rel:
                                        string += '@TITLE "Relative errors as a function of %s"\n' \
                                            % Xms(var).make_greek() 
                                    else:
                                        if noqm:
                                            string += '@TITLE "Absolute values model"\n'
                                        else:
                                            string += '@TITLE "Absolute values for qm method %s"\n' \
                                                % qm_method

                                    subtitle = '@SUBTITLE "Using Rq, Rp: (%s, %s); max_l = %d;' \
                                        %( Rq, Rp, max_l )
                                    
                                    if dist:
                                        subtitle += 'Distributed;'
                                    else:
                                        subtitle += 'Oxygen-cent;'

                                    if var == 'r':
                                        subtitle += ' Constant %s: %s, %s: %s, %s: %s, %s: %s, %s: %s" \n'\
                                                %(
                                                  Xms("tau").make_greek(), self.opts["tau"]["constant"],
                                                  Xms("theta").make_greek(),   self.opts["theta"]["constant"],
                                                  Xms("rho1").make_greek(),  self.opts["rho1"]["constant"],
                                                  Xms("rho2").make_greek(),  self.opts["rho2"]["constant"],
                                                  Xms("rho3").make_greek(),  self.opts["rho3"]["constant"])
                                    if var == 'tau':
                                        subtitle += ' Constant %s: %s, %s: %s, %s: %s, %s: %s, %s: %s" \n'\
                                                %(
                                                  Xms("r").make_greek(), self.opts["r"]["constant"],
                                                  Xms("theta").make_greek(), self.opts["theta"]["constant"],
                                                  Xms("rho1").make_greek(), self.opts["rho1"]["constant"],
                                                  Xms("rho2").make_greek(), self.opts["rho2"]["constant"],
                                                  Xms("rho3").make_greek(), self.opts["rho3"]["constant"])
                                    if var == 'theta':
                                        subtitle += ' Constant %s: %s, %s: %s, %s: %s, %s: %s, %s: %s" \n'\
                                                %(
                                                  Xms("r").make_greek(), self.opts["r"]["constant"],
                                                  Xms("tau").make_greek(), self.opts["tau"]["constant"],
                                                  Xms("rho1").make_greek(), self.opts["rho1"]["constant"],
                                                  Xms("rho2").make_greek(), self.opts["rho2"]["constant"],
                                                  Xms("rho3").make_greek(), self.opts["rho3"]["constant"])
                                    if var == "rho1":
                                        subtitle += ' Constant %s: %s, %s: %s, %s: %s, %s: %s, %s: %s" \n'\
                                                %(
                                                  Xms("r").make_greek(), self.opts["r"]["constant"],
                                                  Xms("tau").make_greek(), self.opts["theta"]["constant"],
                                                  Xms("theta").make_greek(), self.opts["theta"]["constant"],
                                                  Xms("rho2").make_greek(), self.opts["rho2"]["constant"],
                                                  Xms("rho3").make_greek(), self.opts["rho3"]["constant"])
                                    if var == "rho2":
                                        subtitle += ' Constant %s: %s, %s: %s, %s: %s, %s: %s, %s: %s" \n'\
                                                %(
                                                  Xms("r").make_greek(), self.opts["r"]["constant"],
                                                  Xms("tau").make_greek(), self.opts["theta"]["constant"],
                                                  Xms("theta").make_greek(), self.opts["theta"]["constant"],
                                                  Xms("rho1").make_greek(), self.opts["rho1"]["constant"],
                                                  Xms("rho3").make_greek(), self.opts["rho3"]["constant"])
                                    if var == "rho3":
                                        subtitle += ' Constant %s: %s, %s: %s, %s: %s, %s: %s, %s: %s" \n'\
                                                %(Xms("r").make_greek(), self.opts["r"]["constant"],
                                                  Xms("tau").make_greek(), self.opts["tau"]["constant"],
                                                  Xms("theta").make_greek(), self.opts["theta"]["constant"],
                                                  Xms("rho1").make_greek(), self.opts["rho1"]["constant"],
                                                  Xms("rho2").make_greek(), self.opts["rho2"]["constant"])
                                                
                                    if args.no_subtitle:
                                        subtitle = '@SUBTITLE "\n' 
                                    else:
                                        string += subtitle

                                    string +=  '@VIEW 0.15, 0.10, 1.15, 0.85\n'
                                    string +=  '@LEGEND ON\n'
                                    string +=  '@LEGEND BOX ON\n'
                                    string +=  '@LEGEND BOX FILL OFF\n'
                                    string +=  '@LEGEND LOCTYPE VIEW\n'
                                    string +=  '@LEGEND 0.8, 0.5\n' 
                                    string +=  '@LEGEND CHAR SIZE 1.2\n' 

                                    xvg_label = self.input_style_to_xvg_output( level, prop, comp, dist, max_l,Rq,Rp)
                                    string +=  '@ s%d LEGEND "%s"\n' %( localCounter,xvg_label)
#                            if dist:
#                                string +=  '@ s%d LEGEND "%s, %s, %s, %s, %s"\n' %( localCounter, level, prop, comp, str(max_l), "LoProp") 
#                            else:
#                                string +=  '@ s%d LEGEND "%s"\n' %( localCounter, xvg_label )

                                    if out_AA:
                                        au = "Angstrom"
                                    else:
                                        au = "A.U."
                                    string +=  '@ XAXIS LABEL "%s \[%s\]"\n' % ( greek(var), au )

                                    if rel:
                                        string +=  '@ YAXIS LABEL "Relative error"\n' 
                                    else:
                                        string +=  '@ YAXIS LABEL "Absolute value \[A.U.\]"\n' 

                                    string +=  '@ XAXIS LABEL CHAR SIZE 1.5\n' 
                                    string +=  '@ YAXIS LABEL CHAR SIZE 1.5\n'
                                    string +=  '@ TITLE SIZE 2\n'
                                    string +=  '@ SUBTITLE SIZE 1.0\n'

#add the actual data x and y    

                                    for i in range(len( x )):
                                        try:
                                            string += "%s %.4f\n" %( x[i], y[i][ind1][ind2] )
                                        except IndexError:
                                            print x[i], ind1, ind2
                                            print y[i]
                                            raise SystemExit
                                        except TypeError:
                                            print ind1, ind2 
                                            print y[i][ ind1][ ind2 ]
                                            raise SystemExit
                                    string += '@ SORT s%d X ASCENDING\n' % localCounter
                                    string += '@ s%d LINEWIDTH %d\n' % (localCounter, width )
                                    string += '@ s%d LINESTYLE %d\n' % (localCounter, style )
                                    localCounter += 1
        return string

    def input_style_to_xvg_output( self, level, prop,
            component,
            dist,
            max_l,
            Rq = "0.00",
            Rp = "0.00",
            ):
        """level = zeroth, linear, quadratic
        prop = dipole, alpha, beta
        component = X, Y, Z
        max_l = 1, 2
        return nice label string for xmgrace"""
        
        Rq = "%.1f" %float(Rq)
        Rp = "%.1f" %float(Rp)

        d = " "
        if dist:
            d += "LoProp "

        t_level = "(0) "
        t_component = component
        if level == "linear":
            t_level = "(1) "
        elif level == "quadratic":
            t_level = "(2) "

        if prop == "dipole":
            t_prop = 'p'
        elif prop == "alpha":
            t_component *= 2
            t_prop = r'\xa\0'
        elif prop == "beta":
            t_component *= 2
            t_component += 'Z'
            t_prop = r'\xb\0'
        lab = t_level + t_prop + (r"\s%s\N " % t_component) + d
        lab += '(%s,%s)' %(Rq,Rp)
        return lab

    def get_zeroth_string( self, waters, to_au = True, dist = False ):
        """ Converts list of waters into Olav string for .pot"""

        if dist:
            pass
        else:
            for i in waters:
                if to_au:
                    i.to_au()
                i.center = [ i.o.x, i.o.y, i.o.z ]
            string_zeroth = "AU\n%d 1 0\n" %len(waters)
            for i in waters:
                string_zeroth += "%d %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n" %((i.number, i.center[0], 
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
        return string_zeroth
    def get_linear_string( self, waters ):
        """ Converts list of waters into Olav string for .pot"""
        string_linearizable = "AU\n%d 1 2 1\n" %len(waters)
        for i in waters:
            i.center = [ i.o.x, i.o.y, i.o.z ]
        for i in waters:
            string_linearizable += "%d %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n" %((i.number, i.center[0],
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

        return string_linearizable

    def get_string( self, waters, max_l = 1, pol = 22 , quadratic = 2, dist = False ):
        """ Converts list of waters into Olav string for quadraticlinearizable .pot"""
# If the properties are in distributed form, I. E. starts from Oxygen, then H in +x and H -x
        if dist:
            string = "AU\n%d %d %d %d\n" % ( len(waters)*3,
                    int(max_l), pol, quadratic )
            for i in waters:
                i.set_property_on_each_atom()
                string +=  i.o.potline( max_l = max_l, pol =pol, quadratic= quadratic, dist= dist )
                string +=  i.h1.potline(max_l = max_l, pol =pol, quadratic= quadratic, dist= dist )
                string +=  i.h2.potline(max_l = max_l, pol =pol, quadratic= quadratic, dist= dist )
            return string
        else:
            string = "AU\n%d %d %d %d\n" % ( len(waters),
                    max_l, pol, quadratic )
            for i in waters:
                string +=  "%d %.5f %.5f %.5f " %(
                        i.o.res_id, i.o.x, i.o.y, i.o.z) + i.Property.potline( max_l = max_l, pol =pol, quadratic= quadratic, dist= dist ) + '\n'
            return string

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
            tmp = Water()
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
            cnt += 1
            waters.append( tmp )
        return waters

    def get_props_ccsd( self, file_name ):

        mol_dip = np.zeros(3)
        alpha = np.zeros(  [3,3])
        beta = np.zeros(   [3,3,3])
        beta_dict = {}
        atoms = []
        lab = ["X", "Y", "Z"]

        pat_dipole = re.compile(r'Total Molecular Dipole Moment')
        pat_xyz = re.compile(r'^\s*(\w+)\s+(-*\d*\.+\d+)\s+(-*\d*\.+\d+)\s+(-*\d*\.+\d+) *$')
        pat_pol = re.compile(r'([XYZ])DIPLEN.*total.*:')

        pat_alpha= re.compile(r'([XYZ])DIPLEN.*([XYZ])DIPLEN.*')
        pat_beta=  re.compile(r'([XYZ])DIPLEN.*([XYZ])DIPLEN.*([XYZ])DIPLEN')

# Reading in dipole
        lines = open( file_name ).readlines()
        for i in range(len( lines )):
            if pat_dipole.search( lines[i] ):
                mol_dip[0] = lines[i+5].split()[1]
                mol_dip[1] = lines[i+6].split()[1]
                mol_dip[2] = lines[i+7].split()[1]

# Reading in Alfa 
        for i in open( file_name ).readlines():
            if pat_alpha.search( i ):
                if len(i.split()) < 8:
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
                    A = pat_alpha.search(i).groups(1)[0]
                    B = pat_alpha.search(i).groups(1)[1]
                    alpha[ lab.index( A ) , lab.index( B ) ]  = frac
                    if A == "X" and B == "Y":
                        alpha[ lab.index( B ) , lab.index( A ) ]  = frac
                    if A == "X" and B == "Z":
                        alpha[ lab.index( B ) , lab.index( A ) ]  = frac
                    if A == "Y" and B == "Z":
                        alpha[ lab.index( B ) , lab.index( A ) ]  = frac
#For Beta
        for i in open( file_name ).readlines():
            if pat_beta.search( i ):
                if len(i.split()) >8:
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

                    lab1 = pat_beta.search(i).groups(1)[0]
                    lab2 = pat_beta.search(i).groups(1)[1]
                    lab3 = pat_beta.search(i).groups(1)[2]

                    beta_dict[ lab1 + lab2 + lab3 ] = frac
        for i, l1 in enumerate(lab):
            for j, l2 in enumerate(lab):
                for k, l3 in enumerate(lab):
                    beta[i, j, k] = beta_dict[ l1 + l2 + l3 ]


        return mol_dip, alpha , beta

    def continue_if_not_relevant_var(self, r, tau, theta, rho1, rho2, rho3):
        if self.opts["r"].has_key( "constant" ):
            if r != self.opts["r"]["constant"]:
                return True
        if self.opts["tau"].has_key( "constant" ):
            if tau != self.opts["tau"]["constant"]:
                return True
        if self.opts["theta"].has_key( "constant" ):
            if theta != self.opts["theta"]["constant"]:
                return True
        if self.opts["rho1"].has_key( "constant" ):
            if rho1 != self.opts["rho1"]["constant"]:
                return True
        if self.opts["rho2"].has_key( "constant" ):
            if rho2 != self.opts["rho2"]["constant"]:
                return True
        if self.opts["rho3"].has_key( "constant" ):
            if rho3 != self.opts["rho3"]["constant"]:
                return True
        return False



class Xms( str ):

    def __init__(self, char):
        self.char = char
    def make_greek(self):
        if self.char not in ["r", "tau", "theta", "rho1", "rho2", "rho3"]:
            print "wrong char in XmgraceStyle class, exiting"; raise SystemExit
        if self.char == "r"    :return r"r"            
        if self.char == "tau":return r"\xt\0"  
        if self.char == "theta"  :return r"\xq\0" 
        if self.char == "rho1" :return r"\xr\0\s1\N" 
        if self.char == "rho2" :return r"\xr\0\s2\N" 
        if self.char == "rho3" :return r"\xr\0\s3\N" 

def greek( char ):
    if char == "r"    :return r"r"            
    if char == "tau":return r"\xt\0"  
    if char == "theta"  :return r"\xq\0" 
    if char == "rho1" :return r"\xr\0\s1\N" 
    if char == "rho2" :return r"\xr\0\s2\N" 
    if char == "rho3" :return r"\xr\0\s3\N" 
        
if __name__ == '__main__':

    A = argparse.ArgumentParser( add_help= True)

    A.add_argument( "-model"  , type = str, default = "gaussian", choices = [ "pointdipole", "quadrupole", "gaussian"] )

    A.add_argument( "-molecule"  , type = str, default = "water", choices = [ "water", "methanol"] )
    A.add_argument( "-mol_model"  , type = str, default = "tip3p", choices = [ "tip3p", "spc"] )

    A.add_argument( "-basis"  , type = str, default = "ANOPVDZ" )
    A.add_argument( "-Rq"  ,nargs = '*',  type = str, default = ["0.00001"] )
    A.add_argument( "-Rp"  ,nargs = '*',  type = str, default = ["0.00001"] )
    A.add_argument( "-R"  , nargs = '*',  type = str, default =  ["0.00001"] )

# RELATED TO WATER ERRORS
    A.add_argument( "-water_errors"  , default = False, action ='store_true' )

#Param analyzez 6 variables r, tau, theta, rho_{1,2,3}, so far only parameter
    A.add_argument( "-params"  , default = True )

    A.add_argument( "-dist"   , nargs = '*', type = int, default = [0] )
    A.add_argument( "-no_subtitle"  , default = False, action = 'store_true' )

    A.add_argument( "-qm"   , type = str,  default = "hfqua", dest = "qm_method" )
    A.add_argument( "-rel", action = 'store_true' , default = False)
    A.add_argument( "-in_AA", action = 'store_true' , default = False)
    A.add_argument( "-out_AA", action = 'store_true' , default = False)

#Only write the absolute value calculated with chosen model args.model
    A.add_argument( "-noqm", action = 'store_true' , default = False)

    A.add_argument( "-l", nargs = "*", default = ["zeroth"] )
    A.add_argument( "-p", nargs = "*", default = ["dipole"] )
    A.add_argument( "-c", nargs = "*", default = ["Z"] )

    A.add_argument( "-max_l", type = int, nargs='*', default = [1], choices = [1 , 2] )

    #A.add_argument( "-pol", type = int, default = 22, choices = [1 , 2, 22] )
    #A.add_argument( "-quadratic", type = int , default = 1, choices = [1] )
# TWO MOLS RELATED
    A.add_argument( "-two_mols", default = False , action = 'store_true' ) 
    A.add_argument( "-worst", default = False , action = 'store_true' ) 
    A.add_argument( "-h5", default = False , action = 'store_true' ) 
    A.add_argument( "-out_degree", default = False , action = 'store_true' ) 

    A.add_argument( "-o", dest = "output" , default = "tmp.xvg" , help = "Name of output xmgrace file" )

    A.add_argument( "-vary"     , default = "r" ,)
    A.add_argument( "-vary_r"   , default = False , action = 'store_true' ) 
    A.add_argument( "-vary_tau" , default = False , action = 'store_true' )
    A.add_argument( "-vary_theta"   , default = False , action = 'store_true' ) 
    A.add_argument( "-vary_rho1", default = False , action = 'store_true' ) 
    A.add_argument( "-vary_rho2", default = False , action = 'store_true' ) 
    A.add_argument( "-vary_rho3", default = False , action = 'store_true' ) 

    A.add_argument( "-r"   ,   type = str ) 
    A.add_argument( "-tau" , type = str ) 
    A.add_argument( "-theta"   , type = str ) 
    A.add_argument( "-rho1",   type = str ) 
    A.add_argument( "-rho2",   type = str ) 
    A.add_argument( "-rho3",   type = str ) 

    args = A.parse_args()
    args.c = map( lambda x: x.upper() , args.c )
    c = Calculator()


    if args.r:
        c.opts[ "r" ] = { "constant" : args.r }
    if args.tau:
        c.opts[ "tau" ] = { "constant" : args.tau }
    if args.theta:
        c.opts[ "theta" ] = { "constant" : args.theta }
    if args.rho1:
        c.opts[ "rho1" ] = { "constant" : args.rho1 }
    if args.rho2:
        c.opts[ "rho2" ] = { "constant" : args.rho2 }
    if args.rho3:
        c.opts[ "rho3" ] = { "constant" : args.rho3 }

    if args.vary_r:
        c.opts[ "r" ] = { "vary" : True }
        args.var = "r"
    if args.vary_tau:
        c.opts[ "tau" ] = { "vary" : True }
        args.var = "tau"
    if args.vary_theta:
        c.opts[ "theta" ] = { "vary" : True }
        args.var = "theta"
    if args.vary_rho1:
        c.opts[ "rho1" ] = { "vary" : True }
        args.var = "rho1"
    if args.vary_rho2:
        c.opts[ "rho2" ] = { "vary" : True }
        args.var = "rho2"
    if args.vary_rho3:
        c.opts[ "rho3" ] = { "vary" : True }
        args.var = "rho3"

    if args.water_errors:
        c.errors( molecule = args.molecule,
                model = args.mol_model )
    if args.two_mols:
        #c.two_mols_python(
        #        basis = args.basis,
        #        levels = args.l,
        #        props = args.p,
        #        comps = args.c,
        #        max_ls = args.max_l,
        #        dists = args.dist,
        #        rel = args.rel,
        #        h5 = args.h5,
        #        noqm = args.noqm,
        #        model = args.model,
        #        mol_model = args.mol_model,
        #        in_AA = args.in_AA,
        #        out_AA = args.out_AA,
        #        Rps = args.Rp,
        #        Rqs = args.Rq,
        #        worst = args.worst,
        #        )
        #raise SystemExit

        c.two_mols(basis = args.basis,
                levels = args.l,
                props = args.p,
                comps = args.c,
                max_ls = args.max_l,
                dists = args.dist,
                rel = args.rel,
                h5 = args.h5,
                noqm = args.noqm,
                model = args.model,
                mol_model = args.mol_model,
                in_AA = args.in_AA,
                out_AA = args.out_AA,
                Rqs = args.Rq,
                Rps = args.Rp,
                worst = args.worst,
                )

        string = c.get_xvg_string( 
                levels = args.l,
                props = args.p,
                comps = args.c,
                max_ls = args.max_l,
                in_AA = args.in_AA,
                out_AA = args.out_AA,
                rel = args.rel,
                dists = args.dist,
                noqm = args.noqm,
                Rqs = args.Rq,
                Rps = args.Rp,
                out_degree = args.out_degree,
                )
        open( args.output , 'w').write( string )
