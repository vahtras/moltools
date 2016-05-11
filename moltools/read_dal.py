#!/usr/bin/env python

__all__ = [ 'read_beta_hf_string' ] 

import os,sys, re, argparse, ctypes, multiprocessing, functools
import numpy as np
import math as m

#from particles import *
from matplotlib import pyplot as plt

from .molecules import Atom
from .template import Template

try:
    from applequist.gaussian import *
except ImportError:
    pass

a0 = 0.52917721092
lab = [ "X", "Y", "Z"]
charge_dic = {"H1": 1.0 ,"H2":1.0 , "C1":6.0, "C7":6.0, "H3":1.0,
            "H4":1.0, "H6": 1.0, "H8":1.0, 
            "H9":1.0, "H10": 1.0, "H12":1.0, 
            "O5":8.0, "O11": 8.0,
            "H": 1.0, "C": 6.0, "N": 7.0, "O": 8.0, "S": 16.0}
mass_dict = {"H": 1.008,  "C": 6.0, "N": 7.0, "O": 15.999, "S": 16.0}
freq_dict = {"0.0": "static","0.0238927": "1907_nm", "0.0428227" : "1064_nm",
        "0.0773571" : "589_nm" }
allowed_elements = ( 'H', 'O' )

def polar_to_cartesian( r, tau, theta):
    x, y, z = r* np.sin( theta )*np.cos( tau ) \
           , r* np.sin(  theta )*np.sin( tau )  \
           , r* np.cos(  theta ) 
    return x , y , z
def write_related( args ):

    if args.xyz.endswith(".pdb"):
        name = args.xyz.split(".")[0] + "_" + str(args.waters) + ".mol"
        waters = molecules.Water.read_waters( args.xyz ,
                in_AA = args.xAA, 
                out_AA = args.oAA,
                N_waters = args.waters )

    elif args.xyz.endswith( ".xyz" ):
        name = args.x.split(".")[0] + ".mol"

    f_ = open( name , "w" )

    if args.oAA:
        str_ = "Angstrom"
    else:
        str_ = ""
    f_.write( "ATOMBASIS\n\nComment\nmolecules.Atomtypes=2 Charge=0 Nosymm %s\n" %str_)

    if not args.wat:
        "Can't write to .mol file, didn't read water molecules"
        raise SystemExit

    hCnt = len(waters) * 2
    oCnt = len(waters)

    f_.write( "Charge=1.0 molecules.Atoms=%d Basis=cc-pVDZ\n" % hCnt)

    for i in waters:
        for j in i:
            if j.element == "H":
                f_.write( "%s   %.5f   %.5f   %.5f\n" %( j.element, j.x, j.y, j.z ))

    f_.write( "Charge=8.0 molecules.Atoms=%d Basis=cc-pVDZ\n" % oCnt)
    for i in waters:
        for j in i:
            if j.element == "O":
                f_.write( "%s   %.5f   %.5f   %.5f\n" %( j.element, j.x, j.y, j.z ))
    print "Finished writing mol files %s" %name
    raise SystemExit


def run_argparse( args ):
    A = argparse.ArgumentParser( )

# ----------------------------
# GENERIC VARIABLES
# ----------------------------
#
    A.add_argument("-dal", type= str, default = 'hflin' )
    A.add_argument("-mol", type= str, default = 'tip3p' )
    A.add_argument( "-dist", action = "store_true", default = False )

# ----------------------------
# READ ALPHA
# ----------------------------

    A.add_argument( "-alpha", type = str, )

# ----------------------------
# BETA ANALYSIS RELATED
# ----------------------------

    A.add_argument( "-beta_analysis_par", action = "store_true", default = False )
    A.add_argument( "-beta_analysis", action = "store_true", default = False )
    A.add_argument( "-freq", type = str, default = "0.0",
            choices = ["0.0", "0.0238927", "0.0428227", "0.0773571"] )
    A.add_argument( "-R", type = float, default = 0.000001)
    A.add_argument( "-beta",dest="beta", type = str,help="File that contains QUADRATIC response output with hyperpolarizabilities" ) 

    A.add_argument( "-in_AA", action = "store_true", default = False )
    A.add_argument( "-out_AA", action = "store_true", default = False )
    A.add_argument( "-basis", type= str, nargs = '*', default = "ANOPVDZ" )
    A.add_argument( "-beta_dal", type= str, default = "hfqua_" )
    A.add_argument( "-Ncpu", type= int, default = "4" )
    A.add_argument( "-N_waters", type= int, default = 15 )
    A.add_argument( "-model", default = "tip3p" )

# ----------------------------
# ALPHA ANALYSIS RELATED
# ----------------------------
#
    A.add_argument( "-alpha_analysis", action = "store_true", default = False )

    A.add_argument( "-nums", type = str, nargs = '*',
            default = map(str, range(1,10)) )

    A.add_argument( "-x", type = str, default = ["nums"],
            choices = ["snaps", "nums", "freqs"] )
    A.add_argument( "-y", type = str, default = ["yy"],
            choices = ["xx", "yy", "zz", "mean", "aniso"] )
    A.add_argument( "-freqs", type = str, nargs = '*',
            default = ['0.0', "0.0238927", "0.0428227", "0.0773571"]
            )
    A.add_argument( "-comps", type = str, nargs = '*', default = ["xx", "yy", "zz"],
            choices = ["xx", "yy", "zz", "mean", "aniso"])
    A.add_argument( "-snaps", type = str, nargs = '*',
            default = map(str, range(10)) )
    A.add_argument( "-eps_out", type = str )
    A.add_argument( "-template_freq", type = str,
            choices = ['0.0', "0.0238927", "0.0428227", "0.0773571"]
            )
    A.add_argument( "-hdf", action = "store_true", default = False )
# ----------------------------
# RELATED TO PLOT WINDOW APPEARANCE
# ----------------------------
#
    A.add_argument( "-ymin", type = float, default = -0.10 )
    A.add_argument( "-ymax", type = float, default = 0.10 )

# ----------------------------
# QM GENERATION RELATED
# ----------------------------

    A.add_argument( "-qm_generation", action = "store_true", default = False )

# ----------------------------
# QM ANALYSIS RELATED
# ----------------------------

    A.add_argument( "-qm_analysis", action = "store_true", default = False )

# ----------------------------
# QMMM GENERATION RELATED
# ----------------------------

    A.add_argument( "-qmmm_generation", action = "store_true", default = False )

    A.add_argument( "-potstyle", default = "QMMM",
            choices = ["QMMM", "PEQM"])

    A.add_argument( "-qm_waters", type = int, nargs = '*',
            default = [1] )
    A.add_argument( "-mm_waters", type = int, nargs = '*',
            default = [1] )
    A.add_argument( "-file_type", type = str, default = "pdb" )
    A.add_argument( "-tname", type = str, default = "TIP3P" )
    A.add_argument( "-tmethod", type = str, default = "HF" )
    A.add_argument( "-tbasis", type = str, default = "ANOPVDZ" )

#also share same arguments -snaps -freqs with -alpha_analysis

# ----------------------------
# QMMM ANALYSIS RELATED
# ----------------------------

    A.add_argument( "-qmmm_analysis", action = "store_true", default = False )
    A.add_argument( "-n_qm", type = str, nargs = '*',
            default = map(str, range(1,10)) )
    A.add_argument( "-n_mm", type = str, nargs = '*',
            default = map(str, range(1,101)) )
    A.add_argument( "-potfreqs", type = str, nargs = '*',
            default = ["0.0", "0.0238927", "0.0428227", "0.0773571"] )

# ----------------------------
# WRITE RELATED pdb to mol generation RELATED
# ----------------------------

    A.add_argument("-waters", type = int , default = 4, help = "how many waters to take closest to center atom, default: 4")

    A.add_argument("-v","--verbose", action='store_true' , default = False)

    A.add_argument("-write", nargs='*', default = [],  help = "Supply any which files to write from a selection: pot, xyz" )


    A.add_argument( "-xyz", dest="xyz", type = str, help = 'Coordinate file with water molecules for the output .pot file. [ xyz , pdb ]')

    A.add_argument( "-xAA", default = False ,action='store_true',
            help = 'Default coordinate type in AA or AU in -x input water coordinate file, default: False ')
    A.add_argument( "-oAA", default = False, action='store_true' , help='Default coordinate type AA or AU for -op output potential file, default: "AU"' )


    A.add_argument( "-tw", type = float, default = 0.0 )

    A.add_argument( "-wat", action = 'store_true' , default=  True )

    a = A.parse_args( args[1:] )
    return a

def is_ccsd( filename):
    """ Return true if the filename, which is DALTON .out file, is a quadratic ccsd calculation"""
    pat_ccsd = re.compile(r'FINAL CCSD RESULTS FOR THE FIRST HYPERPOLARIZABILITIES')
    for i in open(filename).readlines():
        if pat_ccsd.search( i ):
            return True
    return False

def read_alpha_hf( fstr, freq = '0.0', in_AA = False, freqs = 1 ):
# If freqs > 1, will return a tuple of all alphas for each frequency
#
# Reading in Alpha tensor
    fre = freq[0:7]
    pat_alpha = re.compile(r'([XYZ])DIPLEN.*([XYZ])DIPLEN.*= *(-?\d*\.{1}\d+D*-?\+?\d*)')
    pat_new_freq = re.compile(r'FREQUENCY.*SECOND ORDER')
    alpha = np.zeros( [3,3,] )
    lab = ['X', 'Y', 'Z', ]

# For every new frequency, will append this one and store alpha in the last
# element, otherwise, first column is first frequency by default
    freqlist = None

    lines = fstr.split('\n')

    for i in lines:
        if pat_new_freq.search( i ):
            if freqlist is None:
                freqlist = []
            freqlist.append( np.zeros( (3,3 )) )

        if pat_alpha.search( i ):
            matched = pat_alpha.search(i).groups()
            if "D" in matched[2]:
                frac = float( matched[2].replace("D","E") )
            else:
                frac = float( matched[2] )

            A = matched[0]
            B = matched[1]

            alpha[ lab.index( A ) , lab.index( B ) ]  = frac
            freqlist[-1][lab.index( A ), lab.index( B ) ] = frac

            if A == "X" and B == "Y":
                alpha[ lab.index( B ) , lab.index( A ) ]  = frac
                freqlist[-1][lab.index( B ), lab.index( A ) ] = frac

            if A == "X" and B == "Z":
                alpha[ lab.index( B ) , lab.index( A ) ]  = frac
                freqlist[-1][lab.index( B ), lab.index( A ) ] = frac

            if A == "Y" and B == "Z":
                alpha[ lab.index( B ) , lab.index( A ) ]  = frac
                freqlist[-1][lab.index( B ), lab.index( A ) ] = frac

    if freqs > 1:
        return freqlist

    return alpha 

def read_energy( fname, calctype = 'HF' ):
    """Return the energy from dalton .out file fname"""

    for line in open(fname).readlines():
        if re.compile(r'.*Final.*energy').match(line):
            return line.split()[-1]

def read_alpha_ccsd( fstr ):

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

    lines = fstr.split('\n')
# Reading in Alfa 
    for i in lines:
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
    return alpha



def read_beta_ccsd( fstr ):

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
    lines = fstr.split('\n')
    for i in range(len( lines )):
        if pat_dipole.search( lines[i] ):
            mol_dip[0] = lines[i+5].split()[1]
            mol_dip[1] = lines[i+6].split()[1]
            mol_dip[2] = lines[i+7].split()[1]

# Reading in Alfa 
    for i in lines:
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
    for i in lines:
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

                beta_dict[ "".join( [lab1 + lab2 + lab3]) ] = frac
    for i, l1 in enumerate(lab):
        for j, l2 in enumerate(lab):
            for k, l3 in enumerate(lab):
                beta[i, j, k] = beta_dict[ l1 + l2 + l3 ]

    return atoms, mol_dip, alpha , beta


def read_beta( fstr, freq = "0.0",  in_AA = False, out_AA = False ):
    nuc_dip = np.zeros(3)
    el_dip = np.zeros(3)
    alpha = np.zeros([3,3])
    beta = np.zeros([3,3,3])
    tmp = []
    atoms = []
    missing = {}
    exists = {}

# Reading in Beta tensor
    fre = str("%.5f" % float(freq))
    lab = ['X', 'Y', 'Z', ]
    pat_beta = re.compile(r'@ B-freq')
    lines = fstr.split('\n')
    for i in lines:
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


def read_beta_hf( fstr, freq = "0.0",  in_AA = False, out_AA = False ):
    nuc_dip = np.zeros(3)
    el_dip = np.zeros(3)
    alpha = np.zeros([3,3])
    beta = np.zeros([3,3,3])
    tmp = []
    atoms = []
    missing = {}
    exists = {}
    lab = ["X", "Y", "Z"]

    pat_Q = re.compile(r'Total charge of the molecule')
    pat_xyz = re.compile(r'^\s*(\w+)\s+(-*\d*\.+\d+)\s+(-*\d*\.+\d+)\s+(-*\d*\.+\d+) *$')
    pat_pol = re.compile(r'([XYZ])DIPLEN.*total.*:')
#Special xyz hack for camb3lyp output from akka dalton to find atoms
    pat_akka_xyz = re.compile(r'^\s*(\w+)\s+:\s+\d\s+x\s+(-*\d*\.+\d+)\s+\d\s+y\s+(-*\d*\.+\d+)\s+\d\s+z\s+(-*\d*\.+\d+) *$')

    pat_labels_xyz = re.compile(r'^\s*(\S+-+\S+)\s+(-*\d*\.+\d+)\s+(-*\d*\.+\d+)\s+(-*\d*\.+\d+) *$')
# Reading in dipole and charge
    lines = fstr.split( '\n' )
    for i in lines:
        if pat_Q.search( i ):
            Q = float(i.split()[-1])
        if pat_xyz.match(i):
            f = pat_xyz.match(i).groups()
            matched = pat_xyz.match(i).groups()
#Skip coordinates in out file that are for MM region from QMMM
            kwargs = { "AA": in_AA, "element" :  matched[0], "x" : matched[1],
                    "y" : matched[2], "z" : matched[3] }
            tmpAtom = molecules.Atom( **kwargs )
            atoms.append( tmpAtom )

        elif pat_akka_xyz.match(i):
            f = pat_akka_xyz.match(i).groups()
            matched = pat_akka_xyz.match(i).groups()
#Skip coordinates in out file that are for MM region from QMMM
            kwargs = { "AA": in_AA, "element" :  matched[0], "x" : matched[1],
                    "y" : matched[2], "z" : matched[3] }
            tmpAtom = molecules.Atom( **kwargs )
            atoms.append( tmpAtom )

        elif pat_labels_xyz.match(i):
            f = pat_labels_xyz.match(i).groups()
            matched = pat_labels_xyz.match(i).groups()
            lab = matched[0]
            if len(lab.split('-')) == 4:
                element = "H"
            else:
                element = lab.split('-')[2][0]
            kwargs = { "AA": in_AA, "element" :  element, "x" : matched[1],
                    "y" : matched[2], "z" : matched[3] }
            tmpAtom = molecules.Atom( **kwargs )
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


#Set center of nuceli charge to 0

    coc = sum([ x.r * charge_dic[x.element] for x in atoms ]) /\
            sum([ charge_dic[x.element] for x in atoms ])
    for i in atoms:
        nuc_dip += charge_dic[ i.element ] * (i.r - coc )

    if in_AA:
# Make sure center of charge is in Atomic units to give correct electronic dipole
        coc /= a0

# Reading in Alfa and Beta tensor

    fre = str("%.5f" % float(freq))
    pat_alpha = re.compile(r'@.*QRLRVE.*([XYZ])DIPLEN.*([XYZ])DIPLEN.*%s' %fre)
    alpha = np.zeros( [3,3,] )
    lab = ['X', 'Y', 'Z', ]

    for i in lines:
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
    for i in lines:
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
    N_el = sum([charge_dic[at.element] for at in atoms]) - Q
    tot_dip = el_dip - coc * N_el

    return atoms, tot_dip, alpha , beta

def read_props_qmmm( file_, freq = "0.0",  in_AA = False ):
    """ Same as read_beta_hf but skips coordinates not in allowd_elements
    """
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
    for i in open( file_ ).readlines():
        if pat_xyz.match(i):
            f = pat_xyz.match(i).groups()
            
            matched = pat_xyz.match(i).groups()
#Skip coordinates in out file that are for MM region from QMMM
            if matched[0] not in allowed_elements:
                continue

            kwargs = { "element" :  matched[0], "x" : matched[1],
                    "y" : matched[2], "z" : matched[3] }
            tmpAtom = molecules.Atom( **kwargs )
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

    fre = str("%.5f" % float(freq))
    pat_alpha = re.compile(r'@.*QRLRVE.*([XYZ])DIPLEN.*([XYZ])DIPLEN.*%s' %fre)
    alpha = np.zeros( [3,3,] )
    lab = ['X', 'Y', 'Z', ]

    for i in open( file_ ).readlines():
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
    for i in open( file_ ).readlines():
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
    if in_AA:
        nuc_dip /= a0
    tot_dip = nuc_dip - el_dip

    return atoms, nuc_dip - el_dip, alpha , beta

def main():
    """ 
    Program reads alpha and beta tensor and dipole moment from DALTON output

    """
    args = run_argparse( sys.argv )

    if args.alpha:
        a = read_alpha( args.alpha, )

    if args.beta_analysis:
        beta_analysis(args, basis = args.basis,
                dal = args.beta_dal, in_AA = args.in_AA,
                out_AA = args.out_AA,
                ncpu = args.Ncpu,
                N_waters = args.N_waters)

    if args.beta_analysis_par:
        run_beta_analysis_par( N_waters = args.N_waters,
                ncpu = args.Ncpu,
                model = args.model )

    if args.alpha_analysis:
        alpha_analysis(args)

    if args.qm_generation:
        qm_generation( 
                qm_waters = args.qm_waters,
                basis = args.basis
                )

    if args.qmmm_generation:
        qmmm_generation( 
                qm_waters = args.qm_waters,
                mm_waters = args.mm_waters,
                potfreqs = args.potfreqs,
                potstyle = args.potstyle,
                basis = args.basis)

    if args.qm_analysis:
        qm_analysis( in_AA = args.in_AA,
                out_AA = args.out_AA )

    if args.qmmm_analysis:
        qmmm_analysis( args )

    if args.write:
        write_related( args )
    
def read_dipole( file_, freq = "0.0",  in_AA = False, out_AA = False ):
    nuc_dip = np.zeros(3)
    el_dip = np.zeros(3)
    alpha = np.zeros([3,3])
    beta = np.zeros([3,3,3])
    tmp = []
    atoms = []
    missing = {}
    exists = {}
    lab = ["X", "Y", "Z"]

    pat_Q = re.compile(r'Total charge of the molecule')
    pat_xyz = re.compile(r'^\s*(\w+)\s+(-*\d*\.+\d+)\s+(-*\d*\.+\d+)\s+(-*\d*\.+\d+) *$')
    pat_pol = re.compile(r'([XYZ])DIPLEN.*total.*:')
#Special xyz hack for camb3lyp output from akka dalton to find atoms
    pat_akka_xyz = re.compile(r'^\s*(\w+)\s+:\s+\d\s+x\s+(-*\d*\.+\d+)\s+\d\s+y\s+(-*\d*\.+\d+)\s+\d\s+z\s+(-*\d*\.+\d+) *$')

    pat_labels_xyz = re.compile(r'^\s*((?!-)\S+)\s+(-*\d*\.+\d+)\s+(-*\d*\.+\d+)\s+(-*\d*\.+\d+) *$')
# Reading in dipole and charge
    for i in open( file_ ).readlines():
        if pat_Q.search( i ):
            Q = float(i.split()[-1])
        if pat_xyz.match(i):
            f = pat_xyz.match(i).groups()
            matched = pat_xyz.match(i).groups()
#Skip coordinates in out file that are for MM region from QMMM
            kwargs = { "AA": in_AA, "element" :  matched[0], "x" : matched[1],
                    "y" : matched[2], "z" : matched[3] }
            tmpAtom = Atom( **kwargs )
            atoms.append( tmpAtom )

        elif pat_akka_xyz.match(i):
            f = pat_akka_xyz.match(i).groups()
            matched = pat_akka_xyz.match(i).groups()
#Skip coordinates in out file that are for MM region from QMMM
            kwargs = { "AA": in_AA, "element" :  matched[0], "x" : matched[1],
                    "y" : matched[2], "z" : matched[3] }
            tmpAtom = Atom( **kwargs )
            atoms.append( tmpAtom )

        elif pat_labels_xyz.match(i):
            f = pat_labels_xyz.match(i).groups()
            matched = pat_labels_xyz.match(i).groups()
            lab = matched[0]
            if len(lab.split('-')) == 4:
                element = "H"
            else:
                try:
                    element = lab.split('-')[2][0]
                except IndexError as e:
                    warnings.warn( 'Occured when finding wrong pattern for .xyz in read_beta_hf_string ' )
                    continue
            kwargs = { "AA": in_AA, "element" :  element, "x" : matched[1],
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


#Set center of nuceli charge to 0
    coc = sum([ x.r * charge_dic[x.element] for x in atoms ]) /\
            sum([ charge_dic[x.element] for x in atoms ])
    for i in atoms:
        nuc_dip += charge_dic[ i.element ] * (i.r - coc )

    if in_AA:
# Make sure center of charge is in Atomic units to give correct electronic dipole
        coc /= a0
    N_el = sum([charge_dic[at.element] for at in atoms]) - Q
    tot_dip = el_dip - coc * N_el

    return tot_dip

def read_props_qmmm( file_, freq = "0.0",  in_AA = False ):
    """ Same as read_beta_hf but skips coordinates not in allowd_elements
    """
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
    for i in open( file_ ).readlines():
        if pat_xyz.match(i):
            f = pat_xyz.match(i).groups()
            
            matched = pat_xyz.match(i).groups()
#Skip coordinates in out file that are for MM region from QMMM
            if matched[0] not in allowed_elements:
                continue

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

# Reading in Alfa and Beta tensor

    fre = str("%.5f" % float(freq))
    pat_alpha = re.compile(r'@.*QRLRVE.*([XYZ])DIPLEN.*([XYZ])DIPLEN.*%s' %fre)
    alpha = np.zeros( [3,3,] )
    lab = ['X', 'Y', 'Z', ]

    for i in open( file_ ).readlines():
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
    for i in open( file_ ).readlines():
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
    if in_AA:
        nuc_dip /= a0
    tot_dip = nuc_dip - el_dip

    return atoms, nuc_dip - el_dip, alpha , beta





def read_beta_hf( file_, freq = "0.0",  in_AA = False, out_AA = False ):
    with open( file_ ) as f:
        return read_beta_hf_string( f.read(), freq = freq,
            in_AA = in_AA, out_AA = out_AA )

def read_beta_hf_string( string_, freq = "0.0",  in_AA = False, out_AA = False, akka = False):
    nuc_dip = np.zeros(3)
    el_dip = np.zeros(3)
    alpha = np.zeros([3,3])
    beta = np.zeros([3,3,3])
    tmp = []
    atoms = []
    missing = {}
    exists = {}
    lab = ["X", "Y", "Z"]

    pat_Q = re.compile(r'Total charge of the molecule')
    pat_xyz = re.compile(r'^\s*(\w+)\s+(-*\d*\.+\d+)\s+(-*\d*\.+\d+)\s+(-*\d*\.+\d+) *$')
    pat_pol = re.compile(r'([XYZ])DIPLEN.*total.*:')
#Special xyz hack for camb3lyp output from akka dalton to find atoms
    if akka:
        pat_akka_xyz = re.compile(r'^\s*(\w+)\s+:\s+\d\s+x\s+(-*\d*\.+\d+)\s+\d\s+y\s+(-*\d*\.+\d+)\s+\d\s+z\s+(-*\d*\.+\d+) *$')
    else:
        pat_akka_xyz = re.compile(r'^(?!a)a')

    pat_labels_xyz = re.compile(r'^\s*((?!-)\S+)\s+(-*\d*\.+\d+)\s+(-*\d*\.+\d+)\s+(-*\d*\.+\d+) *$')
# Reading in dipole and charge
    for i in string_.split('\n'):
        if pat_Q.search( i ):
            Q = float(i.split()[-1])
        if pat_xyz.match(i):
            f = pat_xyz.match(i).groups()
            matched = pat_xyz.match(i).groups()
#Skip coordinates in out file that are for MM region from QMMM
            kwargs = { "AA": in_AA, "element" :  matched[0], "x" : matched[1],
                    "y" : matched[2], "z" : matched[3] }
            tmpAtom = Atom( **kwargs )
            atoms.append( tmpAtom )

        elif pat_akka_xyz.match(i):
            print i
            print 'asdf'
            raise SystemExit
            f = pat_akka_xyz.match(i).groups()
            matched = pat_akka_xyz.match(i).groups()
#Skip coordinates in out file that are for MM region from QMMM
            kwargs = { "AA": in_AA, "element" :  matched[0], "x" : matched[1],
                    "y" : matched[2], "z" : matched[3] }
            tmpAtom = Atom( **kwargs )
            atoms.append( tmpAtom )

        elif pat_labels_xyz.match(i):
            f = pat_labels_xyz.match(i).groups()
            matched = pat_labels_xyz.match(i).groups()
            lab = matched[0]
            if len(lab.split('-')) == 4:
                element = "H"
            else:
                try:
                    element = lab.split('-')[2][0]
                except IndexError as e:
                    warnings.warn( 'Occured when finding wrong pattern for .xyz in read_beta_hf_string ' )
                    continue
            kwargs = { "AA": in_AA, "element" :  element, "x" : matched[1],
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
    remove = []
    for ind, at in enumerate(atoms[:-1]):
        for other in atoms[ind+1:]:
            if at.equal( other ):
                remove.append( other )
    for each in remove:
        atoms.remove( each )


#Set center of nuceli charge to 0
    coc = sum([ x.r * charge_dic[x.element] for x in atoms ]) /\
            sum([ charge_dic[x.element] for x in atoms ])

    for i in atoms:
        nuc_dip += charge_dic[ i.element ] * (i.r - coc )

    if in_AA and not out_AA:
# Make sure center of charge is in Atomic units to give correct electronic dipole
        coc /= a0

# Reading in Alfa and Beta tensor
    fre = str("%.5f" % float(freq))
    pat_alpha = re.compile(r'@.*QRLRVE.*([XYZ])DIPLEN.*([XYZ])DIPLEN.*%s' %fre)
    alpha = np.zeros( [3,3,] )
    lab = ['X', 'Y', 'Z', ]

    for i in string_.split('\n'):
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
    for i in string_.split('\n'):
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
    N_el = sum([charge_dic[at.element] for at in atoms]) - Q
    tot_dip = -el_dip + coc * N_el

    return atoms, tot_dip, alpha , beta




if __name__ == '__main__':
    main()
