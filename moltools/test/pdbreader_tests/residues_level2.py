#!/usr/bin/env python
import os, re, math, sys, argparse, subprocess, numpy, tarfile
from daltools import one, mol, dens, prop, lr
from daltools.util import full, blocked, subblocked, timing
from operator import attrgetter
from applequistbreader import *
import unittest

FILE = os.path.join(os.path.dirname(__file__), 'snapshot1.pdb')

class TestResiduesLevel1( unittest.TestCase ):
    def setUp(self):
        """ Default arguments used for program 
        equivalent of argparser in pdbreader """

        self.S = System.read_protein_from_file( FILE )
        assert len( self.S ) == 1

        [res.connect_residues( ) for res in self.S ]
        [res.build_residues( level = 2 ) for res in self.S ]
#        [res.build_concaps(  ) for res in self.S ]
#        [res.build_bridges(  ) for res in self.S ]

    def test_proline_level1( self, ):
        """ At level 2 all prolines have 26 atoms"""
        for chain in self.S:
            for res in chain.ready_residues:
                if res.n_term:
                    continue
                if res.res_name == "PRO":
                    cnt = 0
                    for atom in res:
                        print res, atom
                        cnt += 1
                    assert cnt == 26

    def test_glycine_level1( self, ):
        """ At level 2 all residues have 19 atoms"""
        for chain in self.S:
            for res in chain.ready_residues:
                if res.n_term:
                    continue
                if res.res_name == "GLY":
                    cnt = 0
                    for atom in res:
                        print res, atom
                        cnt += 1
                    assert cnt == 19
    def test_alanine_level1( self, ):
        """ At level 2 all alanines have 22 atoms"""
        for chain in self.S:
            for res in chain.ready_residues:
                if res.n_term:
                    continue
                if res.res_name == "ALA":
                    cnt = 0
                    for atom in res:
                        print res, atom
                        cnt += 1
                    assert cnt == 22
    def test_serine_level1( self, ):
        """ At level 2 all serines have 23 atoms"""
        for chain in self.S:
            for res in chain.ready_residues:
                if res.n_term:
                    continue
                if res.res_name == "SER":
                    cnt = 0
                    for atom in res:
                        print res, atom
                        cnt += 1
                    assert cnt == 23
    def test_threonine_level1( self, ):
        """ At level 2 all threonines have 26 atoms"""
        for chain in self.S:
            for res in chain.ready_residues:
                if res.n_term:
                    continue
                if res.res_name == "THR":
                    cnt = 0
                    for atom in res:
                        print res, atom
                        cnt += 1
                    assert cnt == 26

    def test_leucine_level1( self, ):
        """ At level 2 all leucine have 31 atoms"""
        for chain in self.S:
            for res in chain.ready_residues:
                if res.n_term:
                    continue
                if res.res_name == "LEU":
                    cnt = 0
                    for atom in res:
                        print res, atom
                        cnt += 1
                    assert cnt == 31

    def test_isoleucine_level1( self, ):
        """ At level 2 all isoleucine have 31 atoms"""
        for chain in self.S:
            for res in chain.ready_residues:
                if res.n_term:
                    continue
                if res.res_name == "ILE":
                    cnt = 0
                    for atom in res:
                        print res, atom
                        cnt += 1
                    assert cnt == 31

    def test_valine_level1( self, ):
        """ At level 2 all valiness have 28 atoms"""
        for chain in self.S:
            for res in chain.ready_residues:
                if res.n_term:
                    continue
                if res.res_name == "VAL":
                    cnt = 0
                    for atom in res:
                        print res, atom
                        cnt += 1
                    assert cnt == 28

    def test_nonbridge_cys_level1( self, ):
        """ At level 2 all non-bridge cys have 23 atoms"""
        for chain in self.S:
            for res in chain.ready_bridges:
                if res.res_name == "CYS":
                    if not res.Bridge:
                        cnt = 0
                        for atom in res:
                            print res, atom
                            cnt += 1
                        assert cnt == 23

    def test_bridge_cys_level1( self, ):
        """ At level 2 all bridge cys have 24 atoms"""
        for chain in self.S:
            for res in chain.ready_bridges:
                if res.res_name == "CYS":
                    if res.Bridge:
                        cnt = 0
                        for atom in res:
                            print res, atom
                            cnt += 1
                        assert cnt == 24

    def test_methionine_level1( self, ):
        """ At level 2 all methionine have 29 atoms"""
        for chain in self.S:
            for res in chain.ready_residues:
                if res.n_term:
                    continue
                if res.res_name == "MET":
                    cnt = 0
                    for atom in res:
                        print res, atom
                        cnt += 1
                    assert cnt == 29

    def test_phenyl_level1( self, ):
        """ At level 2 all phenyl have 32 atoms"""
        for chain in self.S:
            for res in chain.ready_residues:
                if res.n_term:
                    continue
                if res.res_name == "PHE":
                    cnt = 0
                    for atom in res:
                        print res, atom
                        cnt += 1
                    assert cnt == 32

    def test_tyrosine_level1( self, ):
        """ At level 2 all tyrosine have 33 atoms"""
        for chain in self.S:
            for res in chain.ready_residues:
                if res.n_term:
                    continue
                if res.res_name == "TYR":
                    cnt = 0
                    for atom in res:
                        print res, atom
                        cnt += 1
                    assert cnt == 33

    def test_trypto_level1( self, ):
        """ At level 2 all trypto have 36 atoms"""
        for chain in self.S:
            for res in chain.ready_residues:
                if res.n_term:
                    continue
                if res.res_name == "TRP":
                    cnt = 0
                    for atom in res:
                        print res, atom
                        cnt += 1
                    assert cnt == 36

    def test_aspartic_level1( self, ):
        """ At level 2 all aspartic have 24 atoms as deprotonated"""
        for chain in self.S:
            for res in chain.ready_residues:
                if res.n_term:
                    continue
                if res.res_name == "ASP":
                    cnt = 0
                    for atom in res:
                        print res, atom
                        cnt += 1
                    assert cnt == 24

    def test_glutamic_level1( self, ):
        """ At level 2 all glutamic have 27 atoms as deprotonated"""
        for chain in self.S:
            for res in chain.ready_residues:
                if res.n_term:
                    continue
                if res.res_name == "GLU":
                    cnt = 0
                    for atom in res:
                        print res, atom
                        cnt += 1
                    assert cnt == 27

    def test_asparagine_level1( self, ):
        """ At level 2 all asparagine have 26 atoms """
        for chain in self.S:
            for res in chain.ready_residues:
                if res.n_term:
                    continue
                if res.res_name == "ASN":
                    cnt = 0
                    for atom in res:
                        print res, atom
                        cnt += 1
                    assert cnt == 26

    def test_glutamine_level1( self, ):
        """ At level 2 all glutamine have 29 atoms """
        for chain in self.S:
            for res in chain.ready_residues:

                if res.res_name == "GLN":
                    cnt = 0
                    for atom in res:
                        print res, atom
                        cnt += 1
                    if res.c_term:
                        assert cnt == 24
                    elif res.n_term:
                        assert cnt == 28
                    else:
                        assert cnt == 29

    def test_histidine_level1( self, ):
        """ At level 2 all histidine have 23 atoms """
        for chain in self.S:
            for res in chain.ready_residues:
                if res.res_name == "HIS":
                    cnt = 0
                    for atom in res:
                        print res, atom
                        cnt += 1
                    if res.c_term:
                        assert cnt == 27
                    elif res.n_term:
                        assert cnt == 28
                    else:
                        assert cnt == 29

    def test_lysine_level1( self, ):
        """ At level 2 all lysine have 28 atoms """
        for chain in self.S:
            for res in chain.ready_residues:
                if res.res_name == "LYS":
                    cnt = 0
                    for atom in res:
                        print res, atom
                        cnt += 1
                    if res.c_term:
                        assert cnt == 32
                    elif res.n_term:
                        assert cnt == 33
                    else:
                        assert cnt == 34

    def test_arginine_level1( self, ):
        """ At level 2 all arginine have 30 atoms """
        for chain in self.S:
            for res in chain.ready_residues:
                if res.res_name == "ARG":
                    cnt = 0
                    for atom in res:
                        print res, atom
                        cnt += 1
                    if res.c_term:
                        assert cnt == 34
                    elif res.n_term:
                        assert cnt == 35
                    else:
                        assert cnt == 36

if __name__ == '__main__':
    unittest.main()
