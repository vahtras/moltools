#!/usr/bin/env python
# -*- coding: utf-8 -*-

import itertools, re, sys, os, sys, argparse, warnings
from matplotlib import pyplot as plt
import numpy as np
from read_dal import o_filter, read_beta_hf
from molecules import *
import pandas as pd

from pd.gaussian import *

a0 = 0.52917721092

class MyBoolAction( argparse.Action):
    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        super(MyBoolAction, self).__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        print values
        mod = values
        #mod = map( lambda x: x != 0, values )
        setattr(namespace, self.dest, mod)

#Customize label name
def customize_label( yd, p, l, c, lop ):
    label = ""
    c_mod = c
    if p == 'a': c_mod = c*2
    elif p == 'b': c_mod = c*2 + 'z'

    if lop and yd in ['abs_cl', 'rel_cl_qm']:
        label += r"$\mathsf{%s}_{%s}^{(%s)-LoProp}$" %(to_latex(p), c_mod, l)
    elif yd in ['abs_cl', 'rel_cl_qm']:
        label += r"$\mathsf{%s}_{%s}^{(%s)}$" %(to_latex(p), c_mod, l)
    else:
        label += r"$\mathsf{%s}_{%s}$" %(to_latex(p), c_mod)
# For none and qm, set controlling vars that ensures no duplicates are plotted
    if yd == 'abs_none':
        label = '$(' + label.strip('$') + ')^{\mathsf{None}}$'
    if yd == 'abs_cl':
        label = '$(' + label.strip('$') + ')^{\mathsf{Apple}}$'
    if yd == 'abs_qm':
        label = '$(' + label.strip('$') + ')^{\mathsf{Q.M.}}$'
    if yd == 'rel_none_qm':
        label = '$\epsilon\left(' + label.strip('$') + r'\right)^{\mathsf{None}}$'
    if yd == 'rel_cl_qm':
        label = '$\epsilon\left(' + label.strip('$') + r'\right)^{\mathsf{Apple}}$'
        
    return label

def get_style( ydim, prop, level, loprop):
#Customize label style
    ls = '-'
#Color settings
    color = 'black'
    if level ==0:
        color = 'blue'
    if level ==1:
        color = 'red'
    if "none" in ydim:
        color = 'black'
    if "qm" in ydim:
        color = 'black'
#linestyle settings
    if ydim == "abs_none":
        ls = ':'
    elif ydim == "abs_cl":
        ls = '--'
    elif ydim == "abs_qm":
        ls = '-.'
    elif ydim == "rel_none_qm":
        ls = ':'
    elif ydim == "rel_cl_qm":
        ls = '-'
#Thick lines for point dipole, thin line for LoProp
    if loprop == 0:
        lw =2.0
    else:
        lw = 1.0

    return lw, ls, color


def to_latex(char):
    if char not in ["d","a","b","r", "theta", "tau", "rho1", "rho2", "rho3"]:
        print "wrong char in latex greek, exiting"; raise SystemExit
    if char == "d"    :return r"p"
    if char == "a"    :return r"\alpha"
    if char == "b"    :return r"\beta"
    if char == "r"    :return r"r"
    if char == "theta":return r"\theta"
    if char == "tau"  :return r"\tau"
    if char == "rho1" :return r"\rho_{1}"
    if char == "rho2" :return r"\rho_{2}"
    if char == "rho3" :return r"\rho_{3}"

def run_mpl_2(
        outs, 
        x_dim = "r",
        y_dims = ["abs_cl"],
        comps = ["x"],
        props = ["d"],
        levels = [0],
        mol_model = "TIP3P",
        models = ["gaussian"],
        max_ls = [1],
        method = "HF",
        loprop = [True, False],
        basis = "ANOPVDZ",
        freqs = ["0.0"],
        Rqs = [1e-9],
        Rps = [1e-9],
        in_AA = False,
        out_AA = False,
        title = "Fancy title",
        x_label = "X-axis",
        y_label = "Y-axis",
        verbose = False,
        save = False,
        out_angle = False,
        x_lim = None,
        y_lim = None,
        dal_suffix = "hfqua",
        ):


    comps = map( lambda x: x.lower(), comps )
    props = map( lambda x: x.lower(), props )

    ind_map =  dict( r = 0, tau = 1, theta = 2, rho1 = 3, rho2 = 4, rho3 = 5 )
    comp_map =  dict( x=0, y =1, z = 2)
    comp_map_beta =  dict( x=2, y =7, z = 9)
    prop_map =  dict( d =0, a = 1, b = 2)
    level_map =  {"0" : 0, "1" : 1, "2" : 2 }

    pat = re.compile(r'_(.*)-(.*)-(.*)-(.*)-(.*)-(.*).out')

# sets to True later when this property has been plotted to avoid duplicates
    abs_none_plotted = { key : value for key, value in zip(['x','y','z'],[False,False,False]) }
    abs_qm_plotted =  { key : value for key, value in zip(['x','y','z'],[False,False,False]) }
    rel_none_qm_plotted =  { key : value for key, value in zip(['x','y','z'],[False,False,False]) }

    fig, ax = plt.subplots()
    fig.subplots_adjust(right = 0.65 )
    x_lab = '$' + to_latex(x_label) + '$'
    for ind1, (yd, c, p, l, lop, f, max_l, rq, rp) in enumerate(itertools.product( y_dims, comps, props, levels, loprop, freqs,max_ls,  Rqs, Rps )):
        y_lab = y_label
        x = []
        y = []
        if yd == 'abs_none':
            if abs_none_plotted[ c ]:continue
        if yd == 'abs_qm':
            if abs_qm_plotted[ c ]:continue
        if yd == 'rel_none_qm':
            if rel_none_qm_plotted[ c ]:continue
        for fname in outs:
            x.append( map(float, pat.search(fname).groups())[ind_map[ x_dim ] ] )
            atoms, dipole_qm, alpha_qm, beta_qm = read_beta_hf(fname, 
                    freq = f, in_AA = in_AA, out_AA = out_AA)
            clus = Cluster.get_water_cluster( fname.split('_')[-1].replace('.out','.mol'),
                    in_AA = in_AA,
                    out_AA = out_AA,)

            clus.attach_properties( 
                    model = mol_model,
                    method = method,
                    basis = basis,
                    loprop = lop == 1,
                    freq = f )
            clus.set_qm_mm(5)
            #print fname
            #print clus.sum_property['beta'][ comp_map_beta[c] ]
            #print Rotator.square_3_ut(g.beta())[ comp_map_beta[c] ]
            #print Rotator.square_3_ut(beta_qm) [ comp_map_beta[c] ]
            #print '------------'
#Plot only dipole for level 0
            if l == 0:
                g = GaussianQuadrupoleList.from_string( clus.get_qmmm_pot_string( max_l = max_l, pol = 0, hyp = 0) )
                g.set_damping( rq, rp )
                g.solve_scf()
                if p == 'd':
                    cl = g.total_dipole_moment()
                    qm = dipole_qm
                    if yd == "abs_none":
                        y.append( clus.p[comp_map[c]] )
                    elif yd == "abs_qm":
                        y.append( qm[ comp_map[c]] )
                    elif yd == "abs_cl":
                        y.append( cl[ comp_map[c]] )
                    elif yd == "rel_cl_qm":
                        y.append( ((cl - qm )/qm) [ comp_map[c]] )
                    elif yd == "rel_none_qm":
                        y.append( ((clus.p - qm )/qm) [ comp_map[c]] )
#Plot only dipole and alpha for level 1
            if l == 1:
                g = GaussianQuadrupoleList.from_string( clus.get_qmmm_pot_string( max_l = max_l, pol = 2, hyp = 0) )
                g.set_damping( rq, rp )
                g.solve_scf()
                if p == 'd':
                    cl = g.total_dipole_moment()
                    qm = dipole_qm
                    if yd == "abs_none":
                        y.append( clus.p[comp_map[c]] )
                    elif yd == "abs_qm":
                        y.append( qm[ comp_map[c]] )
                    elif yd == "abs_cl":
                        y.append( cl[ comp_map[c]] )
                    elif yd == "rel_cl_qm":
                        y.append( ((cl - qm )/qm) [ comp_map[c]] )
                    elif yd == "rel_none_qm":
                        y.append( ((clus.p - qm )/qm) [ comp_map[c]] )
                if p == 'a':
                    cl = g.alpha().diagonal()
                    qm = np.einsum('ii->i', alpha_qm )
                    if yd == "abs_none":
                        y.append( Rotator.ut_2_square(clus.sum_property['alpha']).diagonal()[ comp_map[c]] )
                    elif yd == "abs_qm":
                        y.append( qm[ comp_map[c]] )
                    elif yd == "abs_cl":
                        y.append( cl[ comp_map[c]] )
                    elif yd == "rel_cl_qm":
                        y.append( ((cl - qm )/qm) [ comp_map[c]] )
                    elif yd == "rel_none_qm":
                        y.append( ((Rotator.ut_2_square(clus.sum_property['alpha']).diagonal() - qm )/qm) [ comp_map[c]] )
#plot all comps for level 2 
            if l == 2:
                g = GaussianQuadrupoleList.from_string( clus.get_qmmm_pot_string( max_l = max_l, pol = 22, hyp = 1) )
                g.set_damping( rq, rp )
                g.solve_scf()
                if p == 'd':
                    cl = g.total_dipole_moment()
                    qm = dipole_qm
                    if yd == "abs_none":
                        y.append( clus.p[comp_map[c]] )
                    elif yd == "abs_qm":
                        y.append( qm[ comp_map[c]] )
                    elif yd == "abs_cl":
                        y.append( cl[ comp_map[c]] )
                    elif yd == "rel_cl_qm":
                        y.append( ((cl - qm )/qm) [ comp_map[c]] )
                    elif yd == "rel_none_qm":
                        y.append( ((clus.p - qm )/qm) [ comp_map[c]] )
                if p == 'a':
                    cl = g.alpha().diagonal()
                    qm = np.einsum('ii->i', alpha_qm )
                    if yd == "abs_none":
                        y.append( Rotator.ut_2_square(clus.sum_property[ 'alpha' ]).diagonal()[comp_map[c]] )
                    elif yd == "abs_qm":
                        y.append( qm[ comp_map[c]] )
                    elif yd == "abs_cl":
                        y.append( cl[ comp_map[c]] )
                    elif yd == "rel_cl_qm":
                        y.append( ((cl - qm )/qm) [ comp_map[c]] )
                    elif yd == "rel_none_qm":
                        y.append( Rotator.ut_2_square(((clus.sum_property['alpha'] - qm )/qm)) [ comp_map[c]] )
                if p == 'b':
                    cl = Rotator.square_3_ut( g.beta() )
                    qm = Rotator.square_3_ut( beta_qm )
                    if yd == "abs_none":
                        y.append( clus.sum_property[ 'beta' ][comp_map_beta[c]] )
                    elif yd == "abs_qm":
                        y.append( qm[ comp_map_beta[c]] )
                    elif yd == "abs_cl":
                        y.append( cl[ comp_map_beta[c]] )
                    elif yd == "rel_cl_qm":
                        y.append( ((cl - qm )/qm) [ comp_map_beta[c]] )
                    elif yd == "rel_none_qm":
                        y.append( ((clus.sum_property['beta'] - qm )/qm) [ comp_map_beta[c]] )
        if y == []:
            continue

# To plot in radians /degrees, angstrom
#       
        if x_dim == "r" and not in_AA and out_AA:
            x = map( lambda y: y*a0, x )
            if not "Angstrom" in x_lab:
                x_lab += " in Angstrom"
        elif x_dim == "r" and not out_AA:
            if not "Atomic" in x_lab:
                x_lab += " in Atomic units"
        if x_dim != "r" and out_angle:
            x = map( lambda y: y/np.pi *180, x )
            if not "Degrees" in x_lab:
                x_lab += " in Degrees"
        elif x_dim != "r" and not out_angle:
            if not "Radians" in x_lab:
                x_lab += " in Radians"

# Do the plotting
        print "Plotting: " , yd, c, p, l, lop, f, rq, rp

        lw, ls, color = get_style( yd, p, l, lop )

        label =  customize_label( yd, p, l, c, lop )

# For none and qm, set controlling vars that ensures no duplicates are plotted
        if yd == 'abs_none':
            abs_none_plotted[c] = True
        if yd == 'abs_cl':
            pass
        if yd == 'abs_qm':
            abs_qm_plotted[c] = True
        if yd == 'rel_none_qm':
            rel_none_qm_plotted[c] = True
        if yd == 'rel_cl_qm':
            pass
        if "rel" in yd: y_lab = "Relative Error"

#Sort the data and plot
        srt = np.array(sorted(zip (x,y )))
        x, y = srt[:,0], srt[:,1]
        x, y = map( np.array, [x, y] )

        ax.plot( x, y, label = label ,color = color, lw=lw,linestyle = ls,)
        

#Default x and y limits, overwritten later by args
    if x_dim == 'r':
        if out_AA:
            ax.set_xlim(1.5,5.5)
        else:
            ax.set_xlim(3,10)
    else:
        if out_angle:
            ax.set_xlim(0,90)
        else:
            ax.set_xlim(0, np.pi)
    if x_lim is not None:
        ax.set_xlim( x_lim[0], x_lim[1] )
    if y_lim is not None:
        ax.set_ylim( y_lim[0], y_lim[1] )

    if "rel" in y_dims:
        ax.set_ylim(-1, 1 )

    warnings.filterwarnings('error')
    try:
        leg = ax.legend( bbox_to_anchor = (1.01,0.5), loc=6, fontsize = 17,
                bbox_transform= ax.transAxes )
    except Warning:
        print "No data to plot, exititing ... " ; raise SystemExit
#setting title, xlabel and ylabel
    title = ""
    title += " ".join( [ mol_model, method, basis] )
    ax.set_title( title )
    ax.set_xlabel( x_lab )
    ax.set_ylabel( y_lab )
    if save:
        plt.savefig( 'plot.pdf' , format = 'pdf' )
        raise SystemExit
    plt.show()

if __name__ == '__main__':
    A = argparse.ArgumentParser()
    A.add_argument("-x", type= str, default = 'r',
            choices = ['r', 'tau', 'theta', 'rho1', 'rho2', 'rho3'])
    A.add_argument("-y", nargs = '*', type= str, default = ['rel'],
            choices= ['abs_qm','abs_cl','abs_none',
                'rel_cl_qm','rel_none_qm'])
    A.add_argument("-mol_model", type=str , default = "tip3p" )
    A.add_argument("-method", type=str , default = "HF" )
    A.add_argument("-basis", type= str , default = "ANOPVDZ" )

    A.add_argument("-props", nargs = '*', type= str, default = ['d'] )
    A.add_argument("-levels", nargs = '*', type= int, default = [0] )
    A.add_argument("-comps", nargs = '*', type= str, default = ['z'] )
    A.add_argument("-loprop", nargs = '*', type= int , default = [0] )

    A.add_argument("-dal", type=str , default = "hfqua" )

#For filtering output files, only choose based on 6 constant valeus of parameters
    A.add_argument("-r_const", type=float , default = 0.0 )
    A.add_argument("-tau_const", type=float , default = 0.0 )
    A.add_argument("-theta_const", type=float , default = 0.0 )
    A.add_argument("-rho1_const", type=float , default = 0.0 )
    A.add_argument("-rho2_const", type=float , default = 0.0 )
    A.add_argument("-rho3_const", type=float , default = 0.0 )

    A.add_argument("-max_l", type=int , nargs = '*', default = [1] )
    A.add_argument("-Rq", type=float , nargs ='*', default = [1e-9] )
    A.add_argument("-Rp", type=float , nargs ='*', default = [1e-9] )
# Plotting title, etc
    A.add_argument("-title", type=str , default = r"Title" )
    A.add_argument("-x_label", type=str , default = r"X-axis" )
    A.add_argument("-y_label", type=str , default = r"Atomic units" )

    A.add_argument("-x_lim", type=float , nargs = 2 )
    A.add_argument("-y_lim", type=float , nargs = 2 )

    A.add_argument("-v","--verbose", action = 'store_true', default = False )
    A.add_argument("-save", action = 'store_true', default = False )
    A.add_argument("-out_AA", action = 'store_true', default = False )
    A.add_argument("-in_AA", action = 'store_true', default = False )
    A.add_argument("-out_angle", action = 'store_true', default = False )
    args = A.parse_args( sys.argv[1:] )

    o_files = o_filter( [f for f in os.listdir(os.getcwd()) if f.startswith( args.dal ) and f.endswith('.out')], 
            vary = args.x,
            r = args.r_const,
            tau = args.tau_const,
            theta = args.theta_const,
            rho1 = args.rho1_const,
            rho2 = args.rho2_const,
            rho3 = args.rho3_const,
            )
    run_mpl_2( o_files,
            x_dim = args.x,
            y_dims = args.y,
            props = args.props ,
            levels = args.levels,
            comps = args.comps,
            verbose = args.verbose,
            loprop = args.loprop,
            title = args.title,
            x_label = args.x,
            y_label = args.y_label,
            save = args.save,
            in_AA = args.in_AA,
            out_AA = args.out_AA,
            out_angle = args.out_angle,
            x_lim = args.x_lim,
            y_lim = args.y_lim,
            method = args.method.upper(),
            basis = args.basis.upper(),
            mol_model = args.mol_model.upper() ,
            dal_suffix = args.dal,
            Rqs = args.Rq,
            Rps = args.Rp,
            max_ls = args.max_l,
            )

