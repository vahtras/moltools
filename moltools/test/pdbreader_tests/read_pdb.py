#!/usr/bin/env python
import os, re, math, sys, argparse, subprocess, numpy, tarfile
from daltools import one, mol, dens, prop, lr
from daltools.util import full, blocked, subblocked, timing
from operator import attrgetter
import unittest

from applequistbreader import *

class ReadPdb( unittest.TestCase ):

    def setUp(self):
        self.f = "snapshot1.pdb"

    def test_pdbfile(self):
        """ Test loading a pdb file """
        Pdbfile( self.f )

    def test_system_init(self):
        """ Test creating system from the pdb file """
        P = Pdbfile( self.f )
        S = P.read_protein( )

    def test_find_all_sulfur_bridges(self):
        """ Test to find all sulfur bridges in 4COX protein

        should be 5 in total, i.e. 10 residues should be bridges
        
        """
        P = Pdbfile( self.f )
        S = P.read_protein()
        S.find_sulfur_bridges()
        
        for chain in S:
            bridges = 0
            for acid in chain:
                if acid.Bridge:
                    bridges += 1
            assert bridges == 10

if __name__ == '__main__':
    unittest.main()
