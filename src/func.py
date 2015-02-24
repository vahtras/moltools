#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re, sys, os
from matplotlib import pyplot as plt
import numpy as np
from read_dal import o_filter, read_beta_hf
from gaussian import *
from molecules import *
import pandas as pd

import itertools

o_files = [f for f in os.listdir( os.getcwd()) if f.endswith('.out') ]

def run_mpl_2(
        outs, 
        x_dim = "r",
        y_dims = ["none","rel","qm_abs"],
        comps = ["x"],
        props = ["d"],
        levels = [0],
        mol_model = "TIP3P",
        models = ["gaussian"],
        max_l = [1],
        method = "HF",
        loprop = [True, False],
        basis = "ANOPVDZ",
        freqs = ["0.0"],
        Rqs = [0.0001],
        Rps = [0.0001],
        in_AA = False,
        out_AA = False,
        ):
    comps = map( lambda x: x.lower(), comps )
    props = map( lambda x: x.lower(), props )

    ind_map =  dict( r = 0, tau = 1, theta = 2, rho1 = 3, rho2 = 4, rho3 = 5 )
    comp_map =  dict( x=0, y =1, z = 2)
    comp_map_beta =  dict( x=2, y =7, z = 9)
    prop_map =  dict( d =0, a = 1, b = 2)
    level_map =  {"0" : 0, "1" : 1, "2" : 2 }

    pat = re.compile(r'_(.*)-(.*)-(.*)-(.*)-(.*)-(.*).out')

    fig, ax = plt.subplots()
    for ind1, (yd, c, p, l, lop, f, rq, rp) in enumerate(itertools.product( y_dims, comps, props, levels, loprop, freqs, Rqs, Rps )):
        x = []
        y = []
        print yd, c, p, l, lop, f, rq, rp
        for fname in outs:
            x.append( map(float, pat.search(fname).groups())[ind_map[ x_dim ] ] )
            atoms, dipole_qm, alpha_qm, beta_qm = read_beta_hf(fname, 
                    freq = f, in_AA = in_AA, out_AA = out_AA)
            clus = Cluster.get_water_cluster( fname.split('_')[-1].replace('.out','.mol'),
                    in_AA = in_AA,
                    out_AA = out_AA,)
            clus.attach_properties( model = mol_model,
                    method = method,
                    basis = basis,
                    loprop = lop,
                    freq = f )
            clus.set_qm_mm(2)
            g = GaussianQuadrupoleList.from_string( clus.get_qmmm_pot_string() )
            g.set_damping( rq, rp )
            g.solve_scf()
            #print fname
            #print clus.sum_property['beta'][ comp_map_beta[c] ]
            #print Rotator.square_3_ut(g.beta())[ comp_map_beta[c] ]
            #print Rotator.square_3_ut(beta_qm) [ comp_map_beta[c] ]
            #print '------------'
            
            if l == 0:
                g = GaussianQuadrupoleList.from_string( clus.get_qmmm_pot_string( pol = 0) )
                g.solve_scf()
                if p == 'd':
                    y.append( g.total_dipole_moment() [comp_map[c]] )
        x , y = map( np.array, [x, y] )
        s = pd.Series( y, index = x )
        s.sort()
        s.plot( label = "-".join(map(str, [c,p,l,lop] )) )
        #ax.set_ylim(0,2)
        ax.legend(loc='best')
    plt.show()


if __name__ == '__main__':
    o_files = o_filter ( o_files , vary = "r" )
    run_mpl_2( o_files , y_dims = ["rel",], props = ["d",], comps = ['z'] )
