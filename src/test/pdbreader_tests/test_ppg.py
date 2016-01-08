import os, unittest

FILE = os.path.join(os.path.dirname(__file__), 'ppg.pdb')

from moltools import *
from nose.plugins.attrib import attr

@attr(speed = 'slow' )
class TestConcapsLevel1( unittest.TestCase ):
    def setUp(self):
        """ Default arguments used for program 
        equivalent of argparser in pdbreader """

        self.ch = S = System.from_pdb_string ( open(FILE).read() )

        for chain in self.ch:
            chain.connect_residues()
            for res in chain:
                res.gather_ready( c = True, level = 1 )
                res.gather_ready( r = True, level = 1 )

    def test_prolines_concaps_level1( self, ):
        """ At level 1 all concaps have 6 atoms"""
        for chain in self.ch:
            for res in chain:
                if res.res_name == "PRO":
                    if res.c_term:
                        continue
                    assert len( res.concap ) == 6

    def test_concaps_level1( self, ):
        """ At level 1 all concaps have 6 atoms"""
        for chain in self.ch:
            for res in chain:
                if res.res_name != "PRO":
                    if res.c_term:
                        continue
                    assert len( res.concap ) == 6

    def test_collagen_level_1( self, ):
        """ At level 1, assert how many atoms each concap and ready residue has.
        
        For ready made residues:
        For N-terminal, total = N-2
        For C-terminal, total = N-1

        Res | N
        --------
        PRO | 20
        HYP | 21
        GLY | 13
        ALA | 16

        """
        for chain in self.ch:
            for res in chain:
                if res.res_name == "PRO":
                    if res.n_term:
                        assert len( res.ready ) == 19
                    elif res.c_term:
                        assert len( res.ready ) == 18
                    else:
                        assert len( res.ready ) == 20
                elif res.res_name == "HYP":
                    assert len( res.ready ) == 21
                elif res.res_name == "GLY":
                    if res.c_term:
                        assert len( res.ready ) == 12
                    else:
                        assert len( res.ready ) == 13
                elif res.res_name == "ALA":
                    assert len( res.ready ) == 16

#    def test_bridge_concaps_level1( self, ):
#        """ At level 1 all concaps have 6 atoms"""
#        for chain in self.S:
#            for res in chain.ready_concaps:
#                if res.Bridge:
#                    cnt = 0
#                    for atom in res:
#                        print res, atom
#                        cnt += 1
#                    assert cnt == 6
#
#    def test_bridges_level1( self, ):
#        """ At level 1 all bridges have 4 atoms"""
#
#        for chain in self.S:
#            for res in chain.ready_bridges:
#                cnt = 0
#                for atom in res:
#                    print res, atom
#                    cnt += 1
#                assert cnt == 4
#
if __name__ == '__main__':
    unittest.main(  )
