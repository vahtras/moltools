#!/usr/bin/env python
import os, re, math, sys, argparse, subprocess, numpy, tarfile
from daltools import one, mol, dens, prop, lr
from daltools.util import full, blocked, subblocked, timing
from own import *
from operator import attrgetter
from applequistbreader import *

def setup():
    """ Default arguments used for program """
    global S
    A = argparse.ArgumentParser()
    A.add_argument('-t','--tmpdir',
           default = os.getcwd() , help="Which directory to find < -p pdbfile >")
    A.add_argument('-p', '--pdbfile', type = str, default = None,
            help = "Which PDB file to read from")
    A.add_argument('-w','--write',
           nargs = '*',
           help = 'Toggle writing of files, options: [ mol, con, pro, pot ]')
    A.add_argument('--force', 
           default = False, action = 'store_true',
           help = "set reading of tar.gz files and writing a pot file")

    A.add_argument('-level', type = int, default = 2,
            choices = [1,2,3],
           help = "Level of concap, default = 2")
    A.add_argument('-one_mol', type = int)

    A.add_argument('-max_l', type = int, default = 0)

    A.add_argument('-cutoff', type = float, default = 10.0, help= \
            "Defines how many residues to include near QM region, supply -qm [ resId, [ resId , [..]] for qm region")
    A.add_argument('-qm', type = int, nargs='*')

    A.add_argument('-v' , '-verbose', dest = "verbose", action = 'store_true', default = False)
    A.add_argument('--pol', type = int, default = 1)
    A.add_argument('-write_xyz', nargs = '*', type = str)
    A.add_argument('-potskip', default = False, action = 'store_true' )

    args = A.parse_args()

    args.level = 2

    f = "1YFP_h.pdb"
    P = Pdbfile( f )
    S = P.read_protein()
    assert len( S ) == 1

    S.connect_residues_in_chains()
    S.find_sulfur_bridges()
    S.build_concaps_in_chains( args )
    S.build_bridges_in_chains( args )

def test_prolines_concaps_level2():

    """ At level 2 all concaps have 12 atoms"""
    for chain in S:
        for res in chain.ready_concaps:
            if res.res_name == "PRO":
                cnt = 0
                for atom in res:
                    print res, atom
                    cnt += 1
                assert cnt == 12


def test_cro_concap_level2():
    for chain in S:
        for res in chain.ready_concaps:
            if res.res_name == "CRO":
                cnt = 0
                for atom in res:
                    print res, atom
                    cnt += 1
                assert cnt == 12

def test_concaps_level2():
    for chain in S:
        for res in chain.ready_concaps:
            if res.res_name != "PRO":
                cnt = 0
                for atom in res:
                    print res, atom
                    cnt += 1
                assert cnt == 12

def test_bridge_concaps_level2():
    for chain in S:
        for res in chain.ready_concaps:
            if res.Bridge:
                cnt = 0
                for atom in res:
                    print res, atom
                    cnt += 1
                assert cnt == 12

def test_bridges_level2():
    """ At level 1 all bridges have 10 atoms"""

    for chain in S:
        for res in chain.ready_bridges:
            cnt = 0
            for atom in res:
                print res, atom
                cnt += 1
            assert cnt == 10

if __name__ == '__main__':
    setup()
