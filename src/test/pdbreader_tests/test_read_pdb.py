#!/usr/bin/env python
import os, re, math, sys, argparse, subprocess, numpy, tarfile
from daltools import one, mol, dens, prop, lr
from daltools.util import full, blocked, subblocked, timing
from operator import attrgetter
from pdbreader import *
import unittest

FILE = os.path.join(os.path.dirname(__file__), 'snapshot1.pdb')

@unittest.skip('Skip due to being too slow')
class ReadPdb( unittest.TestCase ):

    def setUp(self):

        self.S = System.read_protein_from_file ( FILE )

        assert len( self.S ) == 1

    def test_find_all_sulfur_bridges(self):
        """ Test to find all sulfur bridges in 4COX protein

        should be 5 in total, i.e. 10 residues should be bridges
        
        """
        self.S.find_sulfur_bridges()
        
        for chain in self.S:
            bridges = 0
            for acid in chain:
                if acid.Bridge:
                    bridges += 1
            assert bridges == 10

if __name__ == '__main__':
    unittest.main()
