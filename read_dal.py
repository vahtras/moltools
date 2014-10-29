#!/usr/bin/env python

import os,sys, re, argparse
import numpy as np
import fractions as fr
import math as m

import ut

from particles import *
from gaussian import *

from molecules import Atom, Water, Property, Cluster
from template import Template
from analysis import Analysis

from matplotlib import pyplot as plt

import h5py
#from calculator import *


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

def write_related( args ):

    if args.xyz.endswith(".pdb"):
        name = args.xyz.split(".")[0] + "_" + str(args.waters) + ".mol"
        waters = Water.read_waters( args.xyz ,
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
    f_.write( "ATOMBASIS\n\nComment\nAtomtypes=2 Charge=0 Nosymm %s\n" %str_)

    if not args.wat:
        "Can't write to .mol file, didn't read water molecules"
        raise SystemExit

    hCnt = len(waters) * 2
    oCnt = len(waters)

    f_.write( "Charge=1.0 Atoms=%d Basis=cc-pVDZ\n" % hCnt)

    for i in waters:
        for j in i:
            if j.element == "H":
                f_.write( "%s   %.5f   %.5f   %.5f\n" %( j.element, j.x, j.y, j.z ))

    f_.write( "Charge=8.0 Atoms=%d Basis=cc-pVDZ\n" % oCnt)
    for i in waters:
        for j in i:
            if j.element == "O":
                f_.write( "%s   %.5f   %.5f   %.5f\n" %( j.element, j.x, j.y, j.z ))
    print "Finished writing mol files %s" %name
    raise SystemExit

def qmmm_generation( args ):

    n_qm = args.qm_waters
    n_mm = args.mm_waters
    ending = args.file_type

    pdb_files = [ i for i in os.listdir(os.getcwd()) if i.endswith( "0.pdb" ) ]


    str_ = ""
    if args.dist:
        str_ += "_dist"

    for files in pdb_files:
        c = Cluster.get_water_cluster( files , in_AA = True, out_AA = False,
                N_qm =  n_qm , N_mm = n_mm )
        out_mol = "%s_%dqm_%dmm_%s%s.mol" % ( files.rstrip( '.' + ending ),
                n_qm, n_mm, args.potfreq, str_)
        out_pot = "%s_%dqm_%dmm_%s%s.pot" % ( files.rstrip( '.' + ending ),
                n_qm, n_mm, args.potfreq, str_)
        for wat in [mol for mol in c if mol.in_mm ]:
            if args.dist:
                kwargs_dict = Template().get( *("TIP3P", "HF", "PVDZ",
                    args.dist , args.potfreq ))
                for at in wat:
                    Property.add_prop_from_template( at, kwargs_dict )
            else:
                kwargs_dict = Template().get( *("TIP3P", "HF", "PVDZ",
                    args.dist, args.potfreq  ))
                for at in wat:
                    Property.add_prop_from_template( at, kwargs_dict )
            t1, t2, t3  = wat.get_euler()
            Property.transform_ut_properties( wat.h1.Property, t1, t2 ,t3)
            Property.transform_ut_properties( wat.h2.Property, t1, t2 ,t3)
            Property.transform_ut_properties( wat.o.Property,  t1, t2 ,t3)

        #print c.get_pe_pot_string()
        
        open( out_mol , 'w' ).write( c.get_qm_mol_string()     )
        open( out_pot , 'w' ).write( c.get_qmmm_pot_string()   )

        print "wrote: %s %s" %(out_mol, out_pot)
    raise SystemExit

def beta_analysis(args ):
#These are used to create string for olavs dipole list class using templates
    dipole = np.zeros( [3] )
    alpha = np.zeros( [3, 3] )
    beta = np.zeros( [3, 3, 3])

#To be read from -b hfqua_file.out

    dipole_qm = np.zeros( [3] )
    alpha_qm = np.zeros( [3, 3] )
    beta_qm = np.zeros( [3, 3, 3] )

    waters = np.zeros( [] )

#read if quadratic calculation is supplied

    if not args.beta:
        print "ERROR: Supply file containing QUADRATIC response -args.beta <FILE>"
        raise SystemExit

    if is_ccsd( args.beta ):
        atoms, dipole_qm , alpha_qm , beta_qm = read_beta_ccsd( args )
    else:
        atoms, dipole_qm , alpha_qm , beta_qm = read_beta_hf( args.beta, args.freq )
#Explicit printing to stdout for testing, only the model water from linear / quadratic calc is printed

    print "Dipole: "
    for i in range(3):
        print "%.5f" % (dipole_qm[ i ] )

    print "\nAlpha: "
    for i, j in enumerate ( ut.upper_triangular( 2 )) :
        print "%.5f" % (alpha_qm[ j ] )

    print "\nBeta: "
    for i, jk in enumerate ( ut.upper_triangular( 3 )) :
        print "%.5f" % (beta_qm[ jk ] )

#Read coordinates for water molecules where to put properties to
    if not args.xyz:
        print "Supply coordinate file -xyz FILE to read molecules from"
    waters = Water.read_waters( args.xyz , in_AA = args.xAA , out_AA = args.oAA, N_waters = args.waters)

    alpha_qm = Water.square_2_ut( alpha_qm )
    beta_qm = Water.square_3_ut( beta_qm )

# Read in rotation angles for each water molecule follow by transfer of dipole, alpha and beta to coordinates
    for wat in waters:
        if args.dist:
            kwargs_dict = Template().get( *("TIP3P", "HF", "PVDZ",
                args.dist , "0.0"))
            for at in wat:
                Property.add_prop_from_template( at, kwargs_dict )
        else:
            kwargs_dict = Template().get( *("TIP3P", "HF", "PVDZ",
                args.dist, "0.0") )
            for at in wat:
                Property.add_prop_from_template( at, kwargs_dict )
        t1, t2, t3  = wat.get_euler()
        Property.transform_ut_properties( wat.h1.Property, t1, t2 ,t3)
        Property.transform_ut_properties( wat.h2.Property, t1, t2 ,t3)
        Property.transform_ut_properties( wat.o.Property,  t1, t2 ,t3)
                
    static= GaussianQuadrupoleList.from_string( Water.get_string_from_waters( waters, pol = 0, hyper = 0, dist = args.dist , AA = args.oAA ))
    polar = GaussianQuadrupoleList.from_string( Water.get_string_from_waters( waters, pol = 2, hyper = 0, dist = args.dist , AA = args.oAA ))
    hyper = GaussianQuadrupoleList.from_string( Water.get_string_from_waters( waters, pol = 22,
        hyper = 1, dist = args.dist , AA = args.oAA ))

    hyper.set_damp( args.R , args.R  )

    static.solve_scf()
    polar.solve_scf()
    hyper.solve_scf()

    sd = static.total_dipole_moment( dist = args.dist )
    pd = polar.total_dipole_moment(dist = args.dist)
    hd = hyper.total_dipole_moment(dist = args.dist)

    pa =  Water.square_2_ut( polar.alpha() )
    ha =  Water.square_2_ut( hyper.alpha() )
    hb =  Water.square_3_ut( hyper.beta() )

    lab1 = ["X", "Y", "Z"]
    lab2 = ["XX", "XY", "XZ", "YY", "YZ", "ZZ"]
    lab3 = ["XXX", "XXY", "XXZ", "XYY", "XYZ", "XZZ", "YYY", "YYZ", "YZZ", "ZZZ"]

    print "Dip;\t QM,\t Zero,\t Linea,\t Quadratic"
    for i in range(3):
        print "%s:\t %.2f\t %.2f\t %.2f\t %.2f" % (lab1[i], dipole_qm[i], sd[i], pd[i], hd[i] )

    print "\nAlpha;\t QM,\t Linea,\t Quadratic"
    for i, j in enumerate ( ut.upper_triangular( 2 )) :
        print "%s:\t %.2f\t %.2f\t %.2f" % (lab2[i],\
                alpha_qm[i],\
                pa[i],\
                ha[i]  )

    print "\nBeta;\t QM,\t Quadratic"
    for i, jk in enumerate ( ut.upper_triangular( 3 )) :
        print "%s:\t %.2f\t %.2f" % (lab3[i], beta_qm[ i ], hb[i] )

    norm = np.sqrt( np.sum(  (beta_qm - hb)**2 ) )
    print "Square of the sum of component differences: %.5f"  % norm

    hb = Water.ut_3_square(hb)
    beta_qm = Water.ut_3_square( beta_qm )

    hb_z = np.einsum('ijj->i' ,hb )
    qm_z = np.einsum('ijj->i', beta_qm )

    c = Cluster()
    for i in waters:
        c.append(i)

    print "\n\nThe projected beta, a.k.a. beta parallel component"
    print "Quantum mechanical:"
    qm = np.dot( qm_z, dipole_qm ) / np.linalg.norm( dipole_qm )
    print qm
    print "Quadratic model: "
    mm =  np.dot( hb_z, hd ) / np.linalg.norm( hd )
    print mm
    print args.R
    dists = c.min_dist()
    print len(c), dists[0], dists[1], dists[2], dists[3] ,  qm/mm




def qmmm_analysis( args ):

    dipole = np.zeros( [3] )
    alpha = np.zeros( [3, 3] )

    dipole_qm = np.zeros( [3] )
    alpha_qm = np.zeros( [3, 3] )

    waters = Cluster()

# args.mol by default tip3p, change to spc or olav or whatever if other needed
    pat_snap = re.compile(r'%s(\d+)' %args.mol)

# Find all kinds of QM water / MM water combinations
# prefixed by _ and suffixed by qm / mm
# _12qm[text]_15mm[text]
    pat_qm_mm = re.compile(r'_(\d+)qm.*_(\d+)mm')

# Find frequencies of type -0.0444, -.552, 0.07772 
# Will match all floats in string

    pat_freq = re.compile(r'-?(\d*\.\d+)')

# Grab all files for alpha analysis, only plot the specified by args
# -nums [] -snaps [] -freqs []

    freqs = Water.unique([  pat_freq.search(i).group(1) for i in os.listdir(os.getcwd()) if i.endswith('.dal')])


    snaps = Water.unique([ pat_snap.search(i).group(1) for i in os.listdir(os.getcwd()) if i.endswith('.out')])

    N_qm = Water.unique([ pat_qm_mm.search(i).group(1) for i in os.listdir(os.getcwd()) if i.endswith('.out')])


    N_mm = Water.unique([ pat_qm_mm.search(i).group(2) for i in os.listdir(os.getcwd()) if i.endswith('.out')])
    
    if len(snaps) != len( args.snaps ):
        print "WARNING ; supplied args.snaps doesn't match calculated ones"
        print "will give wrong average over snapshots"

#Sort all calculated files by their type
    snaps.sort( key = lambda x: int(x))
    freqs.sort( key = lambda x: float(x))
    N_qm.sort( key = lambda x: int(x))
    N_mm.sort( key = lambda x: int(x))
    dists = ["", "_dist",]

# component dictionary, plots specific alpha components 
    f_to_ind = {"0.0":0, "0.0238927":1, "0.0428227":2, "0.0773571":3}
    comp_to_ind = {"xx":0, "yy":1, "zz":2, "mean":3, "aniso": 4}
    c_to_ind = {"xx":0, "yy":1, "zz":2, "mean":3, "aniso": 4}

# [snap, qm_water, mm_water, freqs, components, dist]
#
    x = np.zeros( (10, 9, 100, 4, 5, 2) )

    for n_qm in N_qm:
        if n_qm not in args.n_qm:
            continue
        for n_mm in N_mm:
            if n_mm not in args.n_mm:
                continue
            for snap in snaps:
                if snap not in args.snaps:
                    continue
                for freq in freqs:
                    if freq not in args.freqs:
                        continue
                    for dist in dists:

                        out = "_".join( [args.dal,"%s"%freq,"%s%s"%(args.mol,snap), "%sqm"%(n_qm), "%smm" %(n_mm), "%s%s.out" % (freq, dist)] )

                        mol = "_".join( ["%s%s"%(args.mol,snap), "%sqm"%(n_qm),
                            "%smm" %(n_mm), "%s%s.mol" %(freq, dist)] )

                        if not os.path.isfile( os.path.join( os.getcwd(), out)):
                            continue

                        if not os.path.isfile( os.path.join( os.getcwd(), mol)):
                            continue

                        snap_ind = int(snap)
                        qm_ind = int(n_qm) - 1
                        mm_ind = int(n_mm) - 1
                        freq_ind = f_to_ind[freq]
                        if dist == "_dist":
                            dist_ind = 1
                        else:
                            dist_ind = 0
#read from typical quadratic response, find QRLRVE with freq in parenthesis
                        alpha_qm = read_alpha( out )

# Relative error in xx, yy, zz components
                        a_qm = einsum('ii->i', alpha_qm ) 
                        x[ snap_ind, qm_ind, mm_ind,
                                f_to_ind[freq], :3 ,dist_ind ] = a_qm

# Relative error for mean alpha
                        a_qm = einsum('ii', alpha_qm ) / 3
                        x[ snap_ind, qm_ind, mm_ind, 
                                f_to_ind[freq], 3 ,dist_ind ] = a_qm

# Relative error for anisotropic alpha
                        a_qm = 0.5 * ( 3 * einsum('ij,ij', alpha_qm, alpha_qm ) - einsum('ii,jj', alpha_qm, alpha_qm ))
                        x[ snap_ind, qm_ind, mm_ind, 
                                f_to_ind[freq], 4 ,dist_ind ] = a_qm

#Average over all snapshots requested, WARNING may give wrong result if no match
    x1 = x.sum( axis = 0 ) / len ( args.snaps )

    if args.hdf:
        h5file = h5py.File("data.h5", 'w')
        h5file["qmmm"] = x
        h5file.close()
    raise SystemExit

# hardcoded for a 4 frequencies and 5 components 
    max_vec = np.zeros( ( 4, 5, ) )
    for i in args.freqs:
        for j in args.comps:
            max_vec[ f_to_ind[i], comp_to_ind[j] ] = max(map( abs, x1[ :,:, f_to_ind[i], comp_to_ind[j],: ] ))

    if args.dist:
        print "LoProp Maximum errors" 
    else:
        print "Non-LoProp Maximum errors" 

    print "xx\t\tyy\t\tzz\t\tmean\t\taniso"
    print max_vec

    lab = [r"$\alpha_{xx}$", r"$\alpha_{yy}$",
            r"$\alpha_{zz}$", r"$\bar{\alpha}$",
            r"$\left(\Delta\alpha\right)^{2}$"]

    title = r'Relative error $\frac{\alpha^{Model}-\alpha_{qm}}{\alpha_{qm}}$'
    sub = "Averaged over %d snapshots; " %len( args.snaps )

    if args.dist:
        sub += "LoProp "
        if args.template_freq:
            sub += "at $\omega$=%s" %args.template_freq
        else:
            sub += "from respective frequency" 

    fig = plt.figure()
    ax = plt.axes([.15,.1,.8,.7])
    plt.figtext(.5,.9,title, fontsize=24, ha='center')
    plt.figtext(.5,.82,sub ,fontsize=16,ha='center')
    ax.set_xlabel('Number of water molecules', size = 14)
    ax.set_ylabel('Relative Error', size = 14)

    ax.set_xlim( 1, 10 )
    ax.set_ylim( args.ymin, args.ymax )

    for i in args.freqs:
        for j in args.comps:
            ax.plot( range(1, len(nums)+1),
                    x1[:, f_to_ind[i] , comp_to_ind[j] ], label = r'%s $\omega = %s$' %( lab[ comp_to_ind[j] ], freqs[f_to_ind[i]] )
                    )

    leg = plt.legend()

    fig = plt.gcf()

    if args.eps_out:
        out = args.eps_out
    else:
        out = '%s_%dwat_%dsnaps' %( freq_dict[ args.freqs[0] ],len(nums), len(snaps) )

    if args.dist:
        out += "_dist.eps"
    else:
        out += ".eps"

    fig.savefig( out , format = 'eps')
    print "Wrote: %s" % out
    raise SystemExit




def alpha_analysis(args ):

#These are used to create string for olavs dipole list class using templates
    dipole = np.zeros( [3] )
    alpha = np.zeros( [3, 3] )

#To be read from -b hfqua_file.out

    dipole_qm = np.zeros( [3] )
    alpha_qm = np.zeros( [3, 3] )

    waters = Cluster()

    pat_snap_mol = re.compile(r'%s(\d+)_(\d+)' % args.mol)

# Grab all files for alpha analysis, only plot the specified by args
# -nums [] -snaps [] -freqs []

    freqs = [ i.split('_')[1].rstrip('.dal') for i in os.listdir(os.getcwd()) if i.endswith('.dal')]

    snaps = Water.unique([ pat_snap_mol.search(i).group(1) for i in os.listdir(os.getcwd()) if i.endswith('.out')])

    nums = Water.unique([ pat_snap_mol.search(i).group(2) for i in os.listdir(os.getcwd()) if i.endswith('.out')])

    if len(snaps) != len( args.snaps ):
        print "WARNING ; supplied args.snaps doesn't match calculated ones"

    snaps.sort()
    nums.sort()
    freqs.sort()

# component dictionary, plots specific alpha components 
    f_to_ind = {"0.0":0, "0.0238927":1, "0.0428227":2, "0.0773571":3}
    comp_to_ind = {"xx":0, "yy":1, "zz":2, "mean":3, "aniso": 4}

    err = Analysis()

#
# [snap, qm_water, mm_water, freqs, components, dist]
    x = np.zeros( (10, 9, 1, 4, 5, 2) )

    for snap in snaps:
        if snap not in args.snaps:
            continue
        for num in nums:
            if num not in args.nums:
                continue
            for freq in freqs:
                if freq not in args.freqs:
                    continue
                for dist in [0, 1]:
                    out = "_".join( [args.dal,"%s"%freq,"%s%s"%(args.mol,snap), "%s.out"%num] )
                    mol = "_".join( ["%s%s"%(args.mol,snap), "%s.mol"%num] )

                    if not os.path.isfile( os.path.join( os.getcwd(), out)):
                        continue
                    if not os.path.isfile( os.path.join( os.getcwd(), mol)):
                        continue
# Created index from snapshot and water number to put in data vector x

                    snap_ind = int(snap)
                    qm_ind = int(num) - 1
                    mm_ind = 0
                    freq_ind = f_to_ind[freq]
                    dist_ind = dist

#read from typical quadratic response, find QRLRVE with freq in parenthesis

                    atoms, dipole_qm , alpha_qm , beta_qm = read_beta_hf( out, freq )

#Read in Water molecules for where to put properties to

                    waters = Water.read_waters( mol , in_AA = args.xAA , out_AA = args.oAA, N_waters = num )

# Calculate rotation angles for each water molecule followed by a transfer of the dipole, alpha and beta 
                    for wat in waters:
                        if args.template_freq:
                            kwargs_dict = Template().get(  \
                                    *( args.tname , args.tmethod,
                                        args.tbasis, dist == 1, args.template_freq ))
                        else:
                            kwargs_dict = Template().get(  \
                                    *( args.tname , args.tmethod,
                                        args.tbasis, dist == 1, freq ))
                        for at in wat:
                            Property.add_prop_from_template( at, kwargs_dict )
                        t1, t2, t3  = wat.get_euler()
                        for at in wat:
                            Property.transform_ut_properties( at.Property, t1, t2, t3 )

                    polar = GaussianQuadrupoleList.from_string( Water.get_string_from_waters( waters, pol = 2, hyper = 0, AA = args.oAA ))

                    polar.solve_scf()

# Relative error in xx, yy, zz components
                    a_qm = einsum('ii->i', alpha_qm ) 
                    a_cl = einsum('ii->i', polar.alpha())
                    e_xx, e_yy, e_zz = (a_cl - a_qm)/ a_qm

                    x[ snap_ind, qm_ind, mm_ind,
                            f_to_ind[freq], :3 ,dist ] = e_xx, e_yy, e_zz

# Relative error for mean alpha
                    a_qm = einsum('ii', alpha_qm ) / 3
                    a_cl = einsum('ii', polar.alpha() ) / 3
                    e_mean = (a_cl - a_qm ) / a_qm
                    x[ snap_ind, qm_ind, mm_ind, 
                            f_to_ind[freq], 3 ,dist ] = e_mean

# Relative error for anisotropic alpha
                    a_qm = 0.5 * ( 3 * einsum('ij,ij', alpha_qm, alpha_qm ) - einsum('ii,jj', alpha_qm, alpha_qm ))
                    a_cl = 0.5 * ( 3 * einsum('ij,ij', polar.alpha(), polar.alpha() ) - einsum('ii,jj', polar.alpha(), polar.alpha() ))
                    e_aniso =  (a_cl - a_qm) / a_qm
                    x[ snap_ind, qm_ind , mm_ind,
                            f_to_ind[freq], 4 ,dist ] = e_aniso

# Put all errors in 
                    err[ ( snap, num, freq ) ] = \
                            [ e_xx, e_yy, e_zz, e_mean, e_aniso ]

    #raise SystemExit
    #x = np.zeros( [len(snaps), len( nums ), len(freqs ), 5 ] )

    #for i in range(len( snaps )):
    #    for j in range(len( nums )):
    #        for k in range(len( freqs )):
    #            for l in range( 5 ):
    #                try:
    #                    x[i, j, k, l ] = err[ (str(snaps[i]), str(nums[j]), str(freqs[k]) ) ][ l ]
    #                except KeyError:
    #                    x[i, j, k, l ] = 0.0

#Average of all snapshots

    print snaps
    x1 = x.sum( axis = 0 ) / len ( args.snaps )

# hardcoded for a 4 frequencies and 5 components 

    lab = [r"$\alpha_{xx}$", r"$\alpha_{yy}$",
            r"$\alpha_{zz}$", r"$\bar{\alpha}$",
            r"$\left(\Delta\alpha\right)^{2}$"]

    title = r'Relative error $\frac{\alpha^{Model}-\alpha_{qm}}{\alpha_{qm}}$'
    sub = "Averaged over %d snapshots; " %len( args.snaps )

    if args.dist:
        sub += "LoProp "
        if args.template_freq:
            sub += "at $\omega$=%s" %args.template_freq
        else:
            sub += "from respective frequency" 

    fig = plt.figure()
    ax = plt.axes([.15,.1,.8,.7])
    plt.figtext(.5,.9,title, fontsize=24, ha='center')
    plt.figtext(.5,.82,sub ,fontsize=16,ha='center')
    ax.set_xlabel('Number of water molecules', size = 14)
    ax.set_ylabel('Relative Error', size = 14)

    ax.set_xlim( 1, 10 )
    ax.set_ylim( args.ymin, args.ymax )

    lop_dict = { 0: "", 1: "LoProp"}
    for i in args.freqs:
        for j in args.comps:
            for k in range(2):
                ax.plot( range(1, 10),
                        x1[:,0, f_to_ind[i] , comp_to_ind[j], k ], label = r'%s $\omega = %s$, %s' %( lab[ comp_to_ind[j] ], freqs[f_to_ind[i]], lop_dict[k] )
                        )

    leg = plt.legend()

    fig = plt.gcf()

    if args.eps_out:
        out = args.eps_out
    else:
        out = '%s_%dwat_%dsnaps' %( freq_dict[ args.freqs[0] ],len(args.nums), len(args.snaps) )

    if args.dist:
        out += "_dist.eps"
    else:
        out += ".eps"

    if args.hdf:
        h5file = h5py.File("alpha.h5", 'w')
        h5file["alpha_dist"] = x
        h5file.close()

    fig.savefig( out , format = 'eps')
    print "Wrote: %s" % out
    raise SystemExit

def run_argparse( args ):
    A = argparse.ArgumentParser( )

# ----------------------------
# GENERIC VARIABLES
# ----------------------------
#
    A.add_argument("-dal", type= str, default = 'hfqua' )
    A.add_argument("-mol", type= str, default = 'tip3p' )
    A.add_argument( "-dist", action = "store_true", default = False )

# ----------------------------
# BETA ANALYSIS RELATED
# ----------------------------

    A.add_argument( "-beta_analysis", action = "store_true", default = False )
    A.add_argument( "-freq", type = str, default = "0.0",
            choices = ["0.0", "0.0238927", "0.0428227", "0.0773571"] )
    A.add_argument( "-R", type = float, default = 0.000001)
    A.add_argument( "-beta",dest="beta", type = str,help="File that contains QUADRATIC response output with hyperpolarizabilities" ) 

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
# QMMM GENERATION RELATED
# ----------------------------

    A.add_argument( "-qmmm_generation", action = "store_true", default = False )
    A.add_argument( "-qm_waters", type = int, default = 1 )
    A.add_argument( "-mm_waters", type = int, default = 1 )
    A.add_argument( "-file_type", type = str, default = "pdb" )
    A.add_argument( "-potfreq", type = str, default = "0.0",
            choices = ["0.0", "0.0238927", "0.0428227", "0.0773571"] )

#also share same arguments -snaps -freqs with -alpha_analysis

# ----------------------------
# QMMM ANALYSIS RELATED
# ----------------------------

    A.add_argument( "-qmmm_analysis", action = "store_true", default = False )
    A.add_argument( "-n_qm", type = str, nargs = '*',
            default = map(str, range(1,10)) )
    A.add_argument( "-n_mm", type = str, nargs = '*',
            default = map(str, range(1,101)) )

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


    A.add_argument( "-tname", type = str, default = "TIP3P" )
    A.add_argument( "-tmethod", type = str, default = "HF" )
    A.add_argument( "-tbasis", type = str, default = "PVDZ" )
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

def read_alpha( file_, freq = '0.0', in_AA = False ):
# Reading in Alpha tensor
    fre = freq[0:7]
    pat_alpha = re.compile(r'([XYZ])DIPLEN.*([XYZ])DIPLEN.*= *(-?\d*\.{1}\d+D*-?\+?\d*)')
    alpha = np.zeros( [3,3,] )
    lab = ['X', 'Y', 'Z', ]
    for i in open( file_ ).readlines():

        if pat_alpha.search( i ):

            matched = pat_alpha.search(i).groups()

            if "D" in matched[2]:
                frac = float( matched[2].replace("D","E") )
            else:
                frac = float( matched[2] )

            A = matched[0]
            B = matched[1]
            alpha[ lab.index( A ) , lab.index( B ) ]  = frac
            if A == "X" and B == "Y":
                alpha[ lab.index( B ) , lab.index( A ) ]  = frac
            if A == "X" and B == "Z":
                alpha[ lab.index( B ) , lab.index( A ) ]  = frac
            if A == "Y" and B == "Z":
                alpha[ lab.index( B ) , lab.index( A ) ]  = frac
    return alpha 

def read_beta_ccsd( args ):

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
    lines = open( args.beta ).readlines()
    for i in range(len( lines )):
        if pat_dipole.search( lines[i] ):
            mol_dip[0] = lines[i+5].split()[1]
            mol_dip[1] = lines[i+6].split()[1]
            mol_dip[2] = lines[i+7].split()[1]
            print mol_dip

# Reading in Alfa 
    for i in open( args.beta ).readlines():
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
    for i in open( args.beta ).readlines():
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


    return atoms, mol_dip, alpha , beta

def read_beta_hf( file_, freq = "0.0",  in_AA = False ):
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

def main():
    """ 
    Program reads alpha and beta tensor and dipole moment from DALTON output

    """
    args = run_argparse( sys.argv )

    if args.beta_analysis:
        beta_analysis(args)

    if args.alpha_analysis:
        alpha_analysis(args)

    if args.qmmm_generation:
        qmmm_generation( args )

    if args.qmmm_analysis:
        qmmm_analysis( args )

    if args.write:
        write_related( args )
    
if __name__ == '__main__':
    main()
