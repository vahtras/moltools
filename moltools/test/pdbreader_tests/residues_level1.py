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

        level = 1

        self.ch = System.read_protein_from_file( FILE )
        for chains in self.ch:
            chains.connect_residues()
            for res in chains:
                res.gather_ready( r = True, level = 1 )

    def test_proline_level1( self, ):
        """ At level 1 all prolines have 20 atoms"""
        for chain in self.ch:
            for res in chain:
                if res.res_name == "PRO":
                    if res.n_term:
                        assert len( res.ready ) == 19
                    elif res.c_term:
                        assert len( res.ready ) == 18
                    else:
                        assert len( res.ready ) == 20

    def test_glycine_level1( self, ):
        """ At level 1 all residues have 13 atoms"""
        for chain in self.ch:
            for res in chain:
                if res.res_name == "GLY":
                    if res.n_term:
                        assert len(res.ready) == 12
                    elif res.c_term:
                        assert len(res.ready) == 11
                    else:
                        assert len(res.ready) == 13
#    def test_alanine_level1( self, ):
#        """ At level 1 all alanines have 16 atoms"""
#        for chain in self.S:
#            for res in chain.ready_residues:
#                if res.n_term:
#                    continue
#                if res.res_name == "ALA":
#                    cnt = 0
#                    for atom in res:
#                        print res, atom
#                        cnt += 1
#                    assert cnt == 16
#    def test_serine_level1( self, ):
#        """ At level 1 all serines have 17 atoms"""
#        for chain in self.S:
#            for res in chain.ready_residues:
#                if res.n_term:
#                    continue
#                if res.res_name == "SER":
#                    cnt = 0
#                    for atom in res:
#                        print res, atom
#                        cnt += 1
#                    assert cnt == 17
#    def test_threonine_level1( self, ):
#        """ At level 1 all threonines have 20 atoms"""
#        for chain in self.S:
#            for res in chain.ready_residues:
#                if res.n_term:
#                    continue
#                if res.res_name == "THR":
#                    cnt = 0
#                    for atom in res:
#                        print res, atom
#                        cnt += 1
#                    assert cnt == 20
#
#    def test_leucine_level1( self, ):
#        """ At level 1 all leucine have 25 atoms"""
#        for chain in self.S:
#            for res in chain.ready_residues:
#                if res.n_term:
#                    continue
#                if res.res_name == "LEU":
#                    cnt = 0
#                    for atom in res:
#                        print res, atom
#                        cnt += 1
#                    assert cnt == 25
#
#    def test_isoleucine_level1( self, ):
#        """ At level 1 all isoleucine have 25 atoms"""
#        for chain in self.S:
#            for res in chain.ready_residues:
#                if res.n_term:
#                    continue
#                if res.res_name == "ILE":
#                    cnt = 0
#                    for atom in res:
#                        print res, atom
#                        cnt += 1
#                    assert cnt == 25
#
#    def test_valine_level1( self, ):
#        """ At level 1 all valiness have 22 atoms"""
#        for chain in self.S:
#            for res in chain.ready_residues:
#                if res.n_term:
#                    continue
#                if res.res_name == "VAL":
#                    cnt = 0
#                    for atom in res:
#                        print res, atom
#                        cnt += 1
#                    assert cnt == 22
#
#    def test_nonbridge_cys_level1( self, ):
#        """ At level 1 all non-bridge cys have 17 atoms"""
#        for chain in self.S:
#            for res in chain.ready_bridges:
#                if res.res_name == "CYS":
#                    if not res.Bridge:
#                        cnt = 0
#                        for atom in res:
#                            print res, atom
#                            cnt += 1
#                        assert cnt == 17
#
#    def test_bridge_cys_level1( self, ):
#        """ At level 1 all bridge cys have 18 atoms"""
#        for chain in self.S:
#            for res in chain.ready_bridges:
#                if res.res_name == "CYS":
#                    if res.Bridge:
#                        cnt = 0
#                        for atom in res:
#                            print res, atom
#                            cnt += 1
#                        assert cnt == 18
#
#    def test_methionine_level1( self, ):
#        """ At level 1 all methionine have 23 atoms"""
#        for chain in self.S:
#            for res in chain.ready_residues:
#                if res.n_term:
#                    continue
#                if res.res_name == "MET":
#                    cnt = 0
#                    for atom in res:
#                        print res, atom
#                        cnt += 1
#                    assert cnt == 23
#
#    def test_phenyl_level1( self, ):
#        """ At level 1 all phenyl have 26 atoms"""
#        for chain in self.S:
#            for res in chain.ready_residues:
#                if res.n_term:
#                    continue
#                if res.res_name == "PHE":
#                    cnt = 0
#                    for atom in res:
#                        print res, atom
#                        cnt += 1
#                    assert cnt == 26
#
#    def test_tyrosine_level1( self, ):
#        """ At level 1 all tyrosine have 27 atoms"""
#        for chain in self.S:
#            for res in chain.ready_residues:
#                if res.n_term:
#                    continue
#                if res.res_name == "TYR":
#                    cnt = 0
#                    for atom in res:
#                        print res, atom
#                        cnt += 1
#                    assert cnt == 27
#
#    def test_trypto_level1( self, ):
#        """ At level 1 all trypto have 30 atoms"""
#        for chain in self.S:
#            for res in chain.ready_residues:
#                if res.n_term:
#                    continue
#                if res.res_name == "TRP":
#                    cnt = 0
#                    for atom in res:
#                        print res, atom
#                        cnt += 1
#                    assert cnt == 30
#
#    def test_aspartic_level1( self, ):
#        """ At level 1 all aspartic have 18 atoms as deprotonated"""
#        for chain in self.S:
#            for res in chain.ready_residues:
#                if res.n_term:
#                    continue
#                if res.res_name == "ASP":
#                    cnt = 0
#                    for atom in res:
#                        print res, atom
#                        cnt += 1
#                    assert cnt == 18
#
#    def test_glutamic_level1( self, ):
#        """ At level 1 all glutamic have 21 atoms as deprotonated"""
#        for chain in self.S:
#            for res in chain.ready_residues:
#                if res.n_term:
#                    continue
#                if res.res_name == "GLU":
#                    cnt = 0
#                    for atom in res:
#                        print res, atom
#                        cnt += 1
#                    assert cnt == 21
#
#    def test_asparagine_level1( self, ):
#        """ At level 1 all asparagine have 20 atoms """
#        for chain in self.S:
#            for res in chain.ready_residues:
#                if res.n_term:
#                    continue
#                if res.res_name == "ASN":
#                    cnt = 0
#                    for atom in res:
#                        print res, atom
#                        cnt += 1
#                    assert cnt == 20
#
#    def test_glutamine_level1( self, ):
#        """ At level 1 all glutamine have 23 atoms """
#        for chain in self.S:
#            for res in chain.ready_residues:
#
#                if res.res_name == "GLN":
#                    cnt = 0
#                    for atom in res:
#                        print res, atom
#                        cnt += 1
#                    if res.c_term:
#                        assert cnt == 21
#                    elif res.n_term:
#                        assert cnt == 22
#                    else:
#                        assert cnt == 23
#
#    def test_histidine_level1( self, ):
#        """ At level 1 all histidine have 23 atoms """
#        for chain in self.S:
#            for res in chain.ready_residues:
#                if res.res_name == "HIS":
#                    cnt = 0
#                    for atom in res:
#                        print res, atom
#                        cnt += 1
#                    if res.c_term:
#                        assert cnt == 21
#                    elif res.n_term:
#                        assert cnt == 22
#                    else:
#                        assert cnt == 23
#
#    def test_lysine_level1( self, ):
#        """ At level 1 all lysine have 28 atoms """
#        for chain in self.S:
#            for res in chain.ready_residues:
#                if res.res_name == "LYS":
#                    cnt = 0
#                    for atom in res:
#                        print res, atom
#                        cnt += 1
#                    if res.c_term:
#                        assert cnt == 26
#                    elif res.n_term:
#                        assert cnt == 27
#                    else:
#                        assert cnt == 28
#
#    def test_arginine_level1( self, ):
#        """ At level 1 all arginine have 30 atoms """
#        for chain in self.S:
#            for res in chain.ready_residues:
#                if res.res_name == "ARG":
#                    cnt = 0
#                    for atom in res:
#                        print res, atom
#                        cnt += 1
#                    if res.c_term:
#                        assert cnt == 28
#                    elif res.n_term:
#                        assert cnt == 29
#                    else:
#                        assert cnt == 30
#
if __name__ == '__main__':
    unittest.main()
