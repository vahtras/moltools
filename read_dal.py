#!/usr/bin/env python

import os,sys, re, argparse
import numpy as np
import fractions as fr
import math as m

import ut

from particles import *
from water import *
from template import Template


#from calculator import *


a0 = 0.52917721092
lab = [ "X", "Y", "Z"]
charge_dic = {"H1": 1.0 ,"H2":1.0 , "H": 1.0, "C": 6.0, "N": 7.0, "O": 8.0, "S": 16.0}
mass_dict = {"H": 1.008,  "C": 6.0, "N": 7.0, "O": 15.999, "S": 16.0}


def get_string( waters, max_l = 1, pol = 22 , hyper = 1, dist = False, AA = False ):
    """ Converts list of waters into Olav string for hyperpolarizable .pot"""
# If the properties are in distributed form, I. E. starts from Oxygen, then H in +x and H -x
    if AA:
        str_ = "AA"
    else:
        str_ = "AU"
    if dist:
        string = "%s\n%d %d %d %d\n" % ( str_, len(waters)*3,
                max_l, pol, hyper )
        for i in waters:
            i.set_property_on_each_atom()
            string +=  i.o.potline( max_l = max_l, pol =pol, hyper= hyper, dist= dist )
            string +=  i.h1.potline(max_l = max_l, pol =pol, hyper= hyper, dist= dist )
            string +=  i.h2.potline(max_l = max_l, pol =pol, hyper= hyper, dist= dist )
        return string
    else:
        string = "%s\n%d %d %d %d\n" % ( str_, len(waters),
                max_l, pol, hyper )
        for i in waters:
            string +=  "%d %.5f %.5f %.5f " %(
                    int(i.o.res_id), i.o.x, i.o.y, i.o.z) + i.Property.potline( max_l = max_l, pol =pol, hyper= hyper, dist= dist ) + '\n'
        return string

def hyperq(vec):
    tmp = []
    vec_new = []
    for i in vec:
        if i[1] not in tmp:
            tmp.append(i[1])
            vec_new.append( i ) 
    return vec_new

def run_argparse( args ):
    A = argparse.ArgumentParser( description = \
            "This program reads alpha and beta from dalton .out files, obtained for an ideal water molecule centered as oxygen at origo and hydrogens in the xz-plane.\n\n It also read coordinates of arbitrary water molecules and transforms the above read properties to their coordinate reference frame, and writes it to a .pot file." ,add_help= True)
    A.add_argument( "-a", dest = "alpha", type = str, help="File that contains LINEAR response output with polarizabilities" )
    A.add_argument( "-al", type = str, default = "22", help="Symmetry res_id alpha [1, 2, 3]" )
    A.add_argument( "-b",dest="beta", type = str,help="File that contains QUADRATIC response output with hyperpolarizabilities" ) 
    A.add_argument( "-bl", type = str, default = "1", help="Symmetry res_id beta [1, 2, 3]" )
    A.add_argument( "-x", type = str, help = 'Coordinate file with water molecules for the output .pot file. [ xyz , pdb ]')
    A.add_argument( "-xAA", default = False ,action='store_true',
            help = 'Default coordinate type in AA or AU in -x input water coordinate file, default: "AU"')
    A.add_argument( "-dl", type = str, default = "1", help="Symmetry res_id dipole [0, 1]" )
    A.add_argument( "-test", action = "store_true", default = False )
    A.add_argument( "-mon", action = "store_true", default = False )
    A.add_argument("-v","--verbose", action='store_true' , default = False)
    A.add_argument("-write", nargs='*', default = [],  help = "Supply any which files to write from a selection: pot, xyz" )
    A.add_argument("-waters", type = int , default = 4, help = "how many waters to take closest to center atom, default: 4")
    A.add_argument( "-op", type = str, default = "conf.pot", help='output name of the .pot file, default: "conf.pot"' )
    A.add_argument( "-oAA", default = False, action='store_true' , help='Default coordinate type AA or AU for -op output potential file, default: "AU"' )

    A.add_argument( "-ox", type = str, default = "conf.xyz", help='output name of the .xyz file, default: "conf.xyz"' )


    A.add_argument( "-com", help="Place point properties on center-of-mass instead of Oxygen",action = 'store_true', default = False)
    A.add_argument( "-template", action = 'store_true', help= "Activate Beta tensor reading from templates, options provided below for choices", default = False)

    A.add_argument( "-tname", type = str, default = "CENTERED", 
            help = "available templates: CENTERED (default), TIP3P, SPC, OLAV")

    A.add_argument( "-tmethod", type = str, default = "HF",
            help = "available methods: HF (default) , B3LYP")

    A.add_argument( "-tbasis", type = str, default = "PVDZ",
            help = "available choices: PVDZ (default, is actually cc-pVDZ), \
                    MIDDLE (ano type) [O: 5s3p2d, H: 3s1p]")

    a = A.parse_args( args[1:] )
    return a

def is_ccsd( filename):
    """ Return true if the filename, which is DALTON .out file, is a quadratic ccsd calculation"""
    pat_ccsd = re.compile(r'FINAL CCSD RESULTS FOR THE FIRST HYPERPOLARIZABILITIES')
    for i in open(filename).readlines():
        if pat_ccsd.search( i ):
            return True
    return False

def read_alpha( args ):
# Reading in Alpha tensor
    pat_alpha = re.compile(r'@ -<< ([XYZ])DIPLEN.*([XYZ])DIPLEN')
    alpha = np.zeros([3,3])
    lab = ["X", "Y", "Z"]
    for i in open( args.a ).readlines():
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

def read_beta_hf( args ):

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
    for i in open( args.beta ).readlines():
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

# Reading in Alfa and Beta tensor
    pat_alpha = re.compile(r'@.*QRLRVE:.*([XYZ])DIPLEN.*([XYZ])DIPLEN')
    for i in open( args.beta ).readlines():
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
    for i in open( args.beta ).readlines():
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
    if args.xAA:
        nuc_dip /= a0
    tot_dip = nuc_dip - el_dip

    return atoms, nuc_dip - el_dip, alpha , beta


def main():
    """ 
    Program reads alpha and beta tensor and dipole moment from DALTON output

    """
    args = run_argparse( sys.argv )

    f_waters = False

#These are used to create string for olavs dipole list class using templates
    dipole = np.zeros( [3] )
    alpha = np.zeros( [3, 3] )
    beta = np.zeros( [3, 3, 3])

#To be read from -b hfqua_file.out

    dipole_qm = np.zeros( [3] )
    alpha_qm = np.zeros( [3, 3] )
    beta_qm = np.zeros( [3, 3, 3] )

    waters = np.zeros( [] )

#read if linear calculation is supplied

    if args.alpha:
        dipole_qm = read_alpha( args )


#read if quadratic calculation is supplied
    if args.beta:
        if is_ccsd( args.beta ):
            atoms, dipole_qm , alpha_qm , beta_qm = read_beta_ccsd( args )
        else:
            atoms, dipole_qm , alpha_qm , beta_qm = read_beta_hf( args )
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


    if args.verbose:
        for i in atoms:
            print i

    if args.verbose:
        for i in range(len(dipole)):
            print "Dipole_%s: %f" %( lab[i], dipole_qm[i] )
        for i in range(len(alpha)):
            for j in range(len(alpha[i])):
                print "Alpha_%s%s: %f" %( lab[i], lab[j] , alpha_qm[i][j] )
        for i in range(3):
            for j in range(3):
                print "Alpha_%s%s: %f" %( lab[i], lab[j] , alpha_qm[i][j] )
        for i in range(len(beta)):
            for j in range(len(beta[i])):
                for k in range(len(beta[j])):
                    print "Beta_%s%s%s: %f" %( lab[i], lab[j] , lab[k] , beta_qm[i][j][k] )


    #if args.template:

#Read coordinates for water molecules where to put properties to
    f_waters = True
    if f_waters:
        waters = Water.read_waters( args.x , in_AA = args.xAA , out_AA = args.oAA, N_waters = args.waters)

    alpha_qm = Water.square_2_ut( alpha_qm )
    beta_qm = Water.square_3_ut( beta_qm )

# Read in rotation angles for each water molecule follow by transfer of dipole, alpha and beta to coordinates
    if f_waters:
        for i in waters:
            kwargs = Template().get_data( "OLAV", "HF", "PVDZ" )
            p = Property.from_template( **kwargs )

            t1, t2, t3  = i.get_euler()
            p.transform_ut_properties( t1, t2 ,t3 )
            i.Property = p

    static = PointDipoleList.from_string( get_string( waters, pol = 0, hyper = 0, dist = False , AA = args.oAA ))
    polar = PointDipoleList.from_string( get_string( waters, pol = 2, hyper = 0, dist = False , AA = args.oAA ))
    hyper = PointDipoleList.from_string( get_string( waters, dist = False , AA = args.oAA ))

    static.solve_scf()
    polar.solve_scf()
    hyper.solve_scf()

    sd = static.total_dipole_moment()
    pd = polar.total_dipole_moment()
    hd = hyper.total_dipole_moment()

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

    hb = Water.ut_3_to_square(hb)
    beta_qm = Water.ut_3_to_square( beta_qm )

    hb_z = np.einsum('ijj->i' ,hb )
    qm_z = np.einsum('ijj->i', beta_qm )

    print "\n\nThe projected beta, a.k.a. beta parallel component"
    print "Quantum mechanical:"
    print np.dot( qm_z, dipole_qm ) / np.linalg.norm( dipole_qm )
    print "Quadratic model: "
    print np.dot( hb_z, hd ) / np.linalg.norm( hd )
    
# Read in rotation angles for each water molecule follow by transfer of dipole, alpha and beta to coordinates
    if f_waters and args.mon:
        for i in waters:
            tmp1 = Templates().getBeta( "MON1" , args.tmethod, args.tbasis )
            dipole1 = np.array( tmp1[0] )
            alpha1 = np.array(  tmp1[1] )
            beta1 = np.array(   tmp1[2] )

            tmp2 = Templates().getBeta( "MON2" , args.tmethod, args.tbasis )
            dipole2 = np.array( tmp2[0] )
            alpha2 = np.array(  tmp2[1] )
            beta2 = np.array(   tmp2[2] )

            i.o.toAA()
            if i.o.x < 0.1:
                dipole = dipole2
                alpha = alpha2
                beta = beta2
            else:
                dipole = dipole1
                alpha = alpha1
                beta = beta1
            i.o.toAU()

            i.dipole = dipole
            i.alpha = alpha
            i.beta = beta
            i.get_euler()
            #i.transfer_dipole()
            #i.transfer_alpha()
            #i.transfer_beta()

#Write to file, potential

#Only write output file if input alpha, beta, and coordinates are given
    if "pot" in args.write:
        write_pot = True
    else:
        write_pot = False

    if "xyz" in args.write:
        write_xyz = True
    else:
        write_xyz = False

    if "mol" in args.write:
        write_mol = True
    else:
        write_mol = False
#Perform some tests
    if args.test:
        if f_waters:
            print "Performing tests of transfers"
            for i in waters:

                print "Water res_id %d \nSquare dipole:" % i.res_id
                print i.square_dipole()

                print "alpha trace: "
                print i.alpha_trace()

                print "Square beta:"
                print i.square_beta()

                print "alpha projected on dipole:"
                print i.alpha_par()


#Inline temporary code to generate static, polarizable and hyperpolarizable string
#Write the mol file for target cluster, if the corresponding qua_ file
#already exist( i.e. calculation has been run, perform analysis on those)

    if write_mol:

        if args.x.endswith(".pdb"):
            name = args.x.split(".")[0] + "_" + str(args.waters) + ".mol"

        elif args.x.endswith( ".xyz" ):
            name = args.x.split(".")[0] + ".mol"

        f_ = open( name , "w" )

        if args.oAA:
            str_ = "Angstrom"
        else:
            str_ = ""
        f_.write( "ATOMBASIS\n\nComment\nAtomtypes=2 Charge=0 Nosymm %s\n" %str_)

        if not f_waters:
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
        raise SystemExit

#Write the xyz file for target cluster
    if write_xyz:
        f_ = open( args.ox , "w" )
        f_.write( str( 3* len(waters)) + "\n\n" )
        if f_waters:
            for i in waters:
                for j in i.atomlist:
                    f_.write( "%s %.5f %.5f %.5f\n" %( j.element, j.x, j.y, j.z ))
# Read in QM dipole moment from dalton .out files if they exist for supplied .xyz/.pdb file
# lin_tip3p3_10.out is the corresponding out file for tip3p md configuration 3 using 10
# water molecules obtained as linear response with dipole moment
#
    #if args.x:
    #    if args.x.endswith( ".pdb" ):
    #        lin_outfile = "lin_" + args.x.split('.')[0]+"_" + str( args.waters ) + ".out"
    #        qua_outfile = "qua_" + args.x.split('.')[0]+"_" + str( args.waters ) + ".out"
    #    elif args.x.endswith( ".xyz" ):
    #        lin_outfile = "lin_" + args.x.split('.')[0] + ".out"
    #        qua_outfile = "qua_" + args.x.split('.')[0] + ".out"
    #    if os.path.isfile( lin_outfile ):
    #        atoms, qm_dipole = read_coords_and_dipole( args, custom_file = lin_outfile )
    #        print "Found QM dipole moment in %s" %lin_outfile
    #        print qm_dipole


    raise SystemExit

# Do olav calculations for the generated .pot file for the supplied .xyz/.pdb file:
    if not f_waters:
        raise SystemExit

    #if args.op:
    #    string  =  open( args.op ).read()

    if f_waters:
        static = PointDipoleList.from_string( string_static )
        polarizable = PointDipoleList.from_string( string_polarizable )
        hyperpolarizable = PointDipoleList.from_string( string_hyperpolarizable )

        static.solve_scf()
        polarizable.solve_scf()
        hyperpolarizable.solve_scf()

#Dipole section

    select = [ (0, 0, 2), (1, 1, 2), (2, 2, 2)]
    ref = [ dipole_qm, alpha_qm, beta_qm ]

    c =  Calculator()
    c.writeLog()


    reference = dipole_qm

    print '\n\nDipole: p_x, p_y, p_z\n'
    print "Quantum Mech: ", reference
    print "static dipole", static.total_dipole_moment()
    print "polarizable dipole", polarizable.total_dipole_moment()
    print "hyperpolarizable dipole", hyperpolarizable.total_dipole_moment()

    print "\n"
    print "Relative error Dipole static:" , \
            [(this-ref)/ref for this, ref in zip(
                    static.total_dipole_moment(),
                            reference )]
    print "Relative error Dipole polarizable:" , \
         [(this-ref)/ref for this, ref in zip(
                    polarizable.total_dipole_moment(),
                            reference )]
    print "Relative error Dipole hyperpolarizable:" , \
        [(this-ref)/ref for this, ref in zip(
                    hyperpolarizable.total_dipole_moment(),
                            reference )]
#Alpha section
    reference = alpha_qm.diagonal()
    print '\n\nAlfa: a_xx, a_yy, a_zz\n'
    print "Quantum Mech: ", reference
    print "static alpha", static.alpha().diagonal()
    print "polarizable alpha", polarizable.alpha().diagonal()
    print "hyperpolarizable alpha", hyperpolarizable.alpha().diagonal()

    print "\n"
    print "Relative error Alpha polarizable:" ,\
            [(this-ref)/ref for this, ref in zip(
                polarizable.alpha().diagonal(), 
                        reference )]
    print "Relative error Alpha hyperpolarizable:",\
            [(this-ref)/ref for this, ref  \
            in zip( hyperpolarizable.alpha().diagonal(), 
                        reference )]



#Beta section
#Relative error for xxz, yyz, zzz

    select = [ (0, 0, 2), (1, 1, 2), (2, 2, 2)]
    reference = [beta_qm[i, j, k] for i, j, k in select]
    print "Relative error for xxz, yyz, zzz"
    print '\n\nBeta: b_xxz, b_yyz, b_zzz\n'
    print "Quantum Mech: ", reference
    print "Static:", [static.beta()[i, j, k] for i, j, k in select]
    print "Polarizable:" , [polarizable.beta()[i, j, k] for i, j, k in select]
    print "Hyperpolarizable:" ,[hyperpolarizable.beta()[i, j, k] for i, j, k in select]

    print "\n"
    print "Relative error Beta:" ,[(this-ref)/ref for this, ref in zip([
                             hyperpolarizable.beta()[i, j, k] for i, j, k in select],
                        reference)]


if __name__ == '__main__':
    main()
