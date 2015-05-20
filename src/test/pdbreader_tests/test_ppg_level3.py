#!/usr/bin/env python
from pdbreader import System
import unittest, os

FILE = os.path.join(os.path.dirname(__file__), 'ppg.pdb')

@unittest.skip('Skip due to being too slow')
class TestConcapsLevel1( unittest.TestCase ):
    def setUp(self):
        """ Default arguments used for program 
        equivalent of argparser in pdbreader """

        self.ch = S = System.read_protein_from_file ( FILE )

        for chain in self.ch:
            chain.connect_residues()
            for res in chain:
                res.gather_ready( c = True, level = 3 )
                res.gather_ready( r = True, level = 3 )

    def test_prolines_concaps_level3( self, ):
        """ 

        At level 3 the N-terminal proline should have:

        -) 12 from its own.
        -) 10 from next residue 
        -) 1 from next next ( capped XH )
        12 + 10 + 1 = 23
                Because of the sequence PPG, Proline will always be 33 as concap

        """

        for chain in self.ch:
            for res in chain:
                if res.n_term and res.res_name == 'PRO':
                    assert len( res.con ) == 23

#    def test_concaps_level3( self, ):
#        """ At level 3 all concaps have 12 atoms"""
#        for chain in self.ch:
#            for res in chain:
#                if res.res_name != "PRO":
#                    if res.c_term:
#                        continue
#                    assert len( res.con ) == 12
#
    def test_readys_level_3( self, ):
        """At level 3 the regular proline should have:

        -) 14 from its own.
            --) 8 from Next + NextNext if it is a glycine.
            --) 11 from Next + NextNext if it is a glycine.

            -) 8 from Prev + PrevPrev if it is glycine
            -) 11 from Prev + PrevPrev if it is glycine

        14 + 8 + 11 = 33"""
        for chain in self.ch:
            for res in chain:
                if res.res_name == 'PRO':
                    if res.Next and res.Prev:
                        if res.Next.Next and res.Prev.Prev:
                            assert len( res.ready ) == 33
                if res.res_name == 'GLY':
                    if res.Next and res.Prev:
                        if res.Next.Next and res.Prev.Prev:
                            assert len( res.ready ) == 29

#        For ready made residues:
#        For N-terminal, total = N-4
#        For C-terminal, total = N-5
#
#        Res | N
#        --------
#        PRO | 26
#        HYP | 27
#        GLY | 19
#        ALA | 22
#
#        """
#        for chain in self.ch:
#            for res in chain:
#                if res.res_name == "PRO":
#                    if res.n_term:
#                        assert len( res.ready ) == 22
#                    elif res.c_term:
#                        assert len( res.ready ) == 21
#                    else:
#                        assert len( res.ready ) == 26
#                elif res.res_name == "HYP":
#                    assert len( res.ready ) == 27
#                elif res.res_name == "GLY":
#                    if res.c_term:
#                        assert len( res.ready ) == 15
#                    else:
#                        assert len( res.ready ) == 19
#                elif res.res_name == "ALA":
#                    assert len( res.ready ) == 22

if __name__ == '__main__':
    unittest.main(  )
