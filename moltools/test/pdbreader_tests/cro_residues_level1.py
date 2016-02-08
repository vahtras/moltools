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

    args.level = 1

    f = "1YFP_h.pdb"
    P = Pdbfile( f )
    S = P.read_protein()
    assert len( S ) == 1

    S.connect_residues_in_chains()
    S.find_sulfur_bridges()
    S.build_residues_in_chains( args )
    S.build_bridges_in_chains( args )

def test_proline_level1():

    """ At level 1 all prolines have 20 atoms"""
    for chain in S:
        for res in chain.ready_residues:
            if res.n_term:
                continue
            if res.res_name == "PRO":
                cnt = 0
                for atom in res:
                    print res, atom
                    cnt += 1
                assert cnt == 20

def test_glycine_level1():
    """ At level 1 all residues have 13 atoms"""
    for chain in S:
        for res in chain.ready_residues:
            if res.n_term:
                continue
            if res.res_name == "GLY":
                cnt = 0
                for atom in res:
                    print res, atom
                    cnt += 1
                assert cnt == 13
def test_alanine_level1():
    """ At level 1 all alanines have 16 atoms"""
    for chain in S:
        for res in chain.ready_residues:
            if res.n_term:
                continue
            if res.res_name == "ALA":
                cnt = 0
                for atom in res:
                    print res, atom
                    cnt += 1
                assert cnt == 16
def test_serine_level1():
    """ At level 1 all serines have 17 atoms"""
    for chain in S:
        for res in chain.ready_residues:
            if res.n_term:
                continue
            if res.res_name == "SER":
                cnt = 0
                for atom in res:
                    print res, atom
                    cnt += 1
                assert cnt == 17
def test_threonine_level1():
    """ At level 1 all threonines have 20 atoms"""
    for chain in S:
        for res in chain.ready_residues:
            if res.n_term:
                continue
            if res.res_name == "THR":
                cnt = 0
                for atom in res:
                    print res, atom
                    cnt += 1
                assert cnt == 20

def test_leucine_level1():
    """ At level 1 all leucine have 25 atoms"""
    for chain in S:
        for res in chain.ready_residues:
            if res.n_term:
                continue
            if res.res_name == "LEU":
                cnt = 0
                for atom in res:
                    print res, atom
                    cnt += 1
                assert cnt == 25

def test_isoleucine_level1():
    """ At level 1 all isoleucine have 25 atoms"""
    for chain in S:
        for res in chain.ready_residues:
            if res.res_name == "ILE":
                cnt = 0
                for atom in res:
                    print res, atom
                    cnt += 1
                if res.c_term:
                    assert cnt == 23
                elif res.n_term:
                    assert cnt == 24
                else:
                    assert cnt == 25

def test_valine_level1():
    """ At level 1 all valiness have 22 atoms"""
    for chain in S:
        for res in chain.ready_residues:
            if res.n_term:
                continue
            if res.res_name == "VAL":
                cnt = 0
                for atom in res:
                    print res, atom
                    cnt += 1
                assert cnt == 22

def test_nonbridge_cys_level1():
    """ At level 1 all non-bridge cys have 17 atoms"""
    for chain in S:
        for res in chain.ready_residues:
            if res.res_name == "CYS":
                if not res.Bridge:
                    cnt = 0
                    for atom in res:
                        print res, atom
                        cnt += 1
                    assert cnt == 17

def test_bridge_cys_level1():
    """ At level 1 all bridge cys have 18 atoms"""
    for chain in S:
        for res in chain.ready_residues:
            if res.res_name == "CYS":
                if res.Bridge:
                    cnt = 0
                    for atom in res:
                        print res, atom
                        cnt += 1
                    assert cnt == 18

def test_methionine_level1():
    """ At level 1 all methionine have 23 atoms"""
    for chain in S:
        for res in chain.ready_residues:
            if res.n_term:
                continue
            if res.res_name == "MET":
                cnt = 0
                for atom in res:
                    print res, atom
                    cnt += 1
                assert cnt == 23

def test_phenyl_level1():
    """ At level 1 all phenyl have 26 atoms"""
    for chain in S:
        for res in chain.ready_residues:
            if res.n_term:
                continue
            if res.res_name == "PHE":
                cnt = 0
                for atom in res:
                    print res, atom
                    cnt += 1
                assert cnt == 26

def test_tyrosine_level1():
    """ At level 1 all tyrosine have 27 atoms"""
    for chain in S:
        for res in chain.ready_residues:
            if res.n_term:
                continue
            if res.res_name == "TYR":
                cnt = 0
                for atom in res:
                    print res, atom
                    cnt += 1
                assert cnt == 27

def test_trypto_level1():
    """ At level 1 all trypto have 30 atoms"""
    for chain in S:
        for res in chain.ready_residues:
            if res.n_term:
                continue
            if res.res_name == "TRP":
                cnt = 0
                for atom in res:
                    print res, atom
                    cnt += 1
                assert cnt == 30

def test_aspartic_level1():
    """ At level 1 all aspartic have 18 atoms as deprotonated"""
    for chain in S:
        for res in chain.ready_residues:
            if res.n_term:
                continue
            if res.res_name == "ASP":
                cnt = 0
                for atom in res:
                    print res, atom
                    cnt += 1
                assert cnt == 18

def test_glutamic_level1():
    """ At level 1 all glutamic have 21 atoms as deprotonated"""
    for chain in S:
        for res in chain.ready_residues:
            if res.n_term:
                continue
            if res.res_name == "GLU":
                cnt = 0
                for atom in res:
                    print res, atom
                    cnt += 1
                assert cnt == 21

def test_asparagine_level1():
    """ At level 1 all asparagine have 20 atoms """
    for chain in S:
        for res in chain.ready_residues:
            if res.n_term:
                continue
            if res.res_name == "ASN":
                cnt = 0
                for atom in res:
                    print res, atom
                    cnt += 1
                assert cnt == 20

def test_glutamine_level1():
    """ At level 1 all glutamine have 23 atoms """
    for chain in S:
        for res in chain.ready_residues:

            if res.res_name == "GLN":
                cnt = 0
                for atom in res:
                    print res, atom
                    cnt += 1
                if res.c_term:
                    assert cnt == 21
                elif res.n_term:
                    assert cnt == 22
                else:
                    assert cnt == 23

def test_histidine_level1():
    """ At level 1 all histidine have 23 atoms """
    for chain in S:
        for res in chain.ready_residues:
            if res.res_name == "HIS":
                cnt = 0
                for atom in res:
                    print res, atom
                    cnt += 1
                if res.c_term:
                    assert cnt == 21
                elif res.n_term:
                    assert cnt == 22
                else:
                    assert cnt == 23

def test_lysine_level1():
    """ At level 1 all lysine have 28 atoms """
    for chain in S:
        for res in chain.ready_residues:
            if res.res_name == "LYS":
                cnt = 0
                for atom in res:
                    print res, atom
                    cnt += 1

                if res.c_term:
                    assert cnt == 26

                elif res.n_term:
                    assert cnt == 27

                else:
                    assert cnt == 28

def test_arginine_level1():
    """ At level 1 all arginine have 30 atoms """
    for chain in S:
        for res in chain.ready_residues:
            if res.res_name == "ARG":
                cnt = 0
                for atom in res:
                    print res, atom
                    cnt += 1
                if res.c_term:
                    assert cnt == 28
                elif res.n_term:
                    assert cnt == 29
                else:
                    assert cnt == 30

def test_cro_level1():
    """ At level 1 cro will have 36 atoms """
    for chain in S:
        for res in chain.ready_residues:
            if res.res_name == "CRO":
                cnt = 0
                for atom in res:
                    print res, atom
                    cnt += 1
                assert cnt == 36

if __name__ == '__main__':
    setup()
