import unittest, os

from nose.plugins.attrib import attr

from nose.tools import nottest
from moltools import Property
from moltools.pdbreader import *

FILE = os.path.join(os.path.dirname(__file__), 'ppg.pdb')
first_res_file = os.path.join(os.path.dirname(__file__), 'A-PRO1.p' )
second_res_file = os.path.join(os.path.dirname(__file__), 'A-PRO2.p' )
first_con_file = os.path.join(os.path.dirname(__file__), 'A-PRO1-con.p')

#@unittest.skip('Skip due to being too slow')
@attr(speed = 'slow' )
class TestConcapsLevel1( unittest.TestCase ):

    def setUp(self):
        """ Default arguments used for program 
        equivalent of argparser in pdbreader """

        self.ch = System.from_pdb_string( open(FILE).read() )[0]
        self.ch.connect_residues()
        for res in self.ch:
            res.gather_ready( r = True, level = 1 )
            res.gather_ready( c = True, level = 1 )

        f_res = self.ch[0].ready
        f_con = self.ch[0].concap
        s_res = self.ch[1].ready
        s_con = self.ch[1].concap

        self.f_res = f_res
        self.f_con = f_con
        self.s_res = s_res
        self.s_con = s_con

        self.first = self.ch[0]


    def test_relevants(self):
        """First and last residue only have 2 relevant residues"""
        tot_res = reduce( lambda a, x: a + len(x.get_relevant_residues()) , self.ch, 0 )
        tot_con = reduce( lambda a, x: a + len(x.get_relevant_concaps()) , self.ch, 0 )
        print tot_res, tot_con
        print len(self.ch) * 2
        assert tot_res == len(self.ch) * 3 - 2
        assert tot_con == len(self.ch) * 2 - 2

#    def test_all_in_first_and_second(self):
#        """All atoms not in the concaps should remain their property in first res
#        
#        The concapped ones should add/subtract from other residues ones encountered
#        in"""
#        first = self.ch[0]
#        first_res = Residue.load( first_res_file )
#        second_res = Residue.load( second_res_file )
#        first_con = Residue.load( first_con_file )
#        first.ready = first_res
#        first.con = first_con
#        first.Next.ready = second_res
#        second = first.Next
#        first._res_id = first_con._res_id = first_res._res_id = 1
#        second._res_id = second_res._res_id = 2
#        first.mfcc_props()
#        second.mfcc_props()
#
##For the first residue
#        ats1 = ['N', 'H1', 'H2', "CA", 'HA',
#                "CB", "HB1", "HB2",
#                "CG", "HG1", "HG2", 
#                "CD", "HD1", "HD2", ]
#        
#        ats_ref = map(lambda x: first.ready.get_atom_by_label( '1-PRO-%s' %x), ats1 )
#        prop_ref = map(lambda x: x.Property.copy(), ats_ref )
#
#        ats_new = map(lambda x: first.get_atom_by_label( '1-PRO-%s' %x), ats1 )
#        prop_new = map(lambda x: x.Property.copy(), ats_new )
#        for new, ref in zip( prop_new, prop_ref ):
#            np.testing.assert_almost_equal( new['charge'], ref['charge'] )
#            np.testing.assert_almost_equal( new['dipole'], ref['dipole'] )
#            np.testing.assert_almost_equal( new['alpha'], ref['alpha'] )
#            np.testing.assert_almost_equal( new['quadrupole'], ref['quadrupole'] )
#            np.testing.assert_almost_equal( new['beta'], ref['beta'] )
#
##For the second residue
#        ats2 = ["CA", 'HA',
#                "CB", "HB1", "HB2",
#                "CG", "HG1", "HG2", 
#                "CD", "HD1", "HD2",]
#        
#        ats_ref = map(lambda x: second.ready.get_atom_by_label( '2-PRO-%s' %x), ats2 )
#        prop_ref = map(lambda x: x.Property.copy(), ats_ref )
#
#        ats_new = map(lambda x: second.get_atom_by_label( '2-PRO-%s' %x), ats2 )
#        prop_new = map(lambda x: x.Property.copy(), ats_new )
#        for new, ref in zip( prop_new, prop_ref ):
#            np.testing.assert_almost_equal( new['charge'], ref['charge'] )
#            np.testing.assert_almost_equal( new['dipole'], ref['dipole'] )
#            np.testing.assert_almost_equal( new['alpha'], ref['alpha'] )
#            np.testing.assert_almost_equal( new['quadrupole'], ref['quadrupole'] )
#            np.testing.assert_almost_equal( new['beta'], ref['beta'] )



#    def test_con_in_first(self):
#        first = self.ch[0]
#        first_res = Residue.load( first_res_file )
#        second_res = Residue.load( second_res_file )
#        first_con = Residue.load( first_con_file )
#        first.ready = first_res
#        first.con = first_con
#        first.Next.ready = second_res
#        second = first.Next
#        first._res_id = first_con._res_id = first_res._res_id = 1
#        second._res_id = second_res._res_id = 2
#        first.mfcc_props()
#        second.mfcc_props()
#
## For Carbon C
#        c_res_1 = first_res.get_atom_by_label( '1-PRO-C' )
#        c_res_2 = second_res.get_atom_by_label( '1-PRO-C' )
#        c_con_1 = first_con.get_atom_by_label( '1-PRO-C' )
#
#        hx_res_2 = second_res.get_atom_by_label( '1-PRO-CA-XH' )
#        hx_con_1 = first_con.get_atom_by_label( '1-PRO-CA-XH' )
#
#        ats = [ c_res_1, c_res_2, hx_res_2, hx_con_1, c_con_1 ]
#        p1 = reduce( lambda a,x:a+x.Property, ats[:3], Property())
#        p2 = reduce( lambda a,x:a+x.Property, ats[3:], Property())
#        prop_ref = p1 - p2
#
#        prop_new = first.get_atom_by_label( '1-PRO-C' ).Property
#        np.testing.assert_almost_equal( prop_ref['charge'], prop_new['charge'] )
#        np.testing.assert_almost_equal( prop_ref['alpha'], prop_new['alpha'] )
#        np.testing.assert_almost_equal( prop_ref['beta'], prop_new['beta'] )
#        np.testing.assert_almost_equal( prop_ref['dipole'], prop_new['dipole'] )
#        np.testing.assert_almost_equal( prop_ref['quadrupole'], prop_new['quadrupole'] )
#
#
## For nitrogen N
#        n_res_1 = first_res.get_atom_by_label( '2-PRO-N' )
#        n_res_2 = second_res.get_atom_by_label( '2-PRO-N' )
#        n_con_1 = first_con.get_atom_by_label( '2-PRO-N' )
#
#        hx1_res_1 = first_res.get_atom_by_label( '2-PRO-CA-XH' )
#        hx2_res_1 = first_res.get_atom_by_label( '2-PRO-CD-XH' )
#
#        hx1_con_1 = first_con.get_atom_by_label( '2-PRO-CA-XH' )
#        hx2_con_1 = first_con.get_atom_by_label( '2-PRO-CD-XH' )
#
#        ats = [ n_res_1, n_res_2, hx1_res_1, hx2_res_1,
#                n_con_1, hx1_con_1, hx2_con_1 ]
#        p1 = reduce( lambda a,x: a + x.Property, ats[:4], Property())
#        p2 = reduce( lambda a,x: a + x.Property, ats[4:], Property())
#        prop_ref = p1 - p2
#
#        prop_new = second.get_atom_by_label( '2-PRO-N' ).Property
#        np.testing.assert_almost_equal( prop_ref['charge'], prop_new['charge'] )
#        np.testing.assert_almost_equal( prop_ref['alpha'], prop_new['alpha'] )
#        np.testing.assert_almost_equal( prop_ref['beta'], prop_new['beta'] )
#        np.testing.assert_almost_equal( prop_ref['dipole'], prop_new['dipole'] )
#        np.testing.assert_almost_equal( prop_ref['quadrupole'], prop_new['quadrupole'] )
#



    def test_first(self):
        assert self.f_res.Next == self.ch[1]
        assert self.f_res.Next.ready == self.s_res
        assert self.f_res.Next.concap == self.s_con
        assert self.f_res == self.ch[0].ready
        assert self.f_con == self.ch[0].concap


    def test_relevant_first(self):
        resi = self.ch[0].get_relevant_residues(  )
        ats = reduce( lambda a,x : a+x, resi )
        assert len(ats) == 39

    @nottest
    def test_mfcc_first(self):
        logging.warning("Test not implemented")
        first = self.ch[0]
        second = self.ch[1]

        first_res = self.ch[0].ready
        first_con = self.ch[0].concap
        second_res = self.ch[1].ready
# Set up properties which should lead to C = 0.8, N = 0.3 in the end

        first_res.get_atom_by_label( "2-PRO-N" ).Property['charge'] = np.array( -0.5)
        first_res.get_atom_by_label( "1-PRO-C" ).Property['charge'] = np.array( -0.2)
        for at in first_res.get_dummy_h():
            at.Property['charge'] = np.array( 0.2 )

        first_con.get_atom_by_label( "2-PRO-N" ).Property['charge'] = np.array( -0.4)
        first_con.get_atom_by_label( "1-PRO-C" ).Property['charge'] = np.array( -0.6)
        for at in first_con.get_dummy_h():
            at.Property['charge'] = np.array( 0.1 )

        second_res.get_atom_by_label( "2-PRO-N" ).Property['charge'] = np.array(0.2)
        second_res.get_atom_by_label( "1-PRO-C" ).Property['charge'] = np.array(0.3)
        for at in second_res.get_dummy_h():
            at.Property['charge'] = np.array( 0.2 )

        first.mfcc_props()
        second.mfcc_props()

        np.testing.assert_almost_equal( first.get_atom_by_label( "1-PRO-C" ).Property['charge'] , 0.8, decimal = 7 )
        np.testing.assert_almost_equal( second.get_atom_by_label( "2-PRO-N" ).Property['charge'] , 0.3, decimal = 7 )

if __name__ == '__main__':
    unittest.main(  )
