#!/usr/bin/env python
import os, re, sys, argparse, tarfile, ctypes, multiprocessing, pickle, logging

#from mayavi import mlab

from mpl_toolkits.mplot3d import Axes3D
from copy import deepcopy

#from daltools import one, mol, dens, prop, lr
#from daltools.util import full, blocked, subblocked, timing
from loprop.loprop import MolFrag, penalty_function, shift_function

import molecules
import utilz
import matplotlib.pyplot as plt

import numpy as np

res_dict = {'ALA':'A', 'VAL':'V', 'ILE':'I','LEU':'L','MET':'M',
        'PHE':'F','TYR':'Y','TRP':'W','SER':'S','THR':'T','ASN':'N', 'CRO':'X1',
        'CRO1':'X1', 'CRO2':"X2",'CRO3':'X3','CRO4':'X4',
        'GLN':'Q','CYS':'C','CH6': 'X1' ,'GLY':'G','PRO':'P','ARG':'R','HIS':'H',
        'LYS':'K','ASP':'D','GLU':'E','SEC':'U','PYL':'U', 'HIP':'Z',
        'HIE':'H','CYX':'C','HSE':'H','HID':'H','HSD':'H','HSP':'H2',"TIP3": 'T3',
        'HIP':'H2','HYP':'PX', 'MOL' : 'X', 'WAT' : 'W1', 'SOL' : 'W1' }
chargeDict = {'ARG':1, "LYS" : 1, "ASP":-1, "GLU":-1,
            'R':1 , 'K':1 ,'H':0, 'H2':1 , 'E':-1 , 'D':-1, 'X2': 1}
proline_dict = { "PRO" : "P", "HYP" : "PX" }
custom_dict = { "CRO2" : "X2", "MOL" : "X3" }

color_dict = { "H" : (1, 1, 1),
        "C" : (0.5,0,0),
        "N" : (0, 0, 0.5),
        "O" : (1, 0, 0),
        "S" : (1,1,0) }

scale_factor_dict = { "H" : 0.5 , "C" : 1.0 , "N" : 1.0 , "O" : 1.0 , "S" : 1.2 }

mass_dict = { "H" : 1.0 , "C" : 12.0 , "N" : 14.0 , "O" : 16.0 , "S" : 32.0 }

pat_xyz = re.compile(r'^\s*(\w|-)+\s+(-*\d*.+\d+)\s+(-*\d*.+\d+)\s+(-*\d*.+\d+) *$')



def all_residues_from_pdb_file( _file,
        in_AA = True,
        out_AA = True):
    with open( _file ) as f:
        return all_residues_from_pdb_string( f.read(),
                in_AA = in_AA,
                out_AA = out_AA)

def all_residues_from_pdb_string( _string,
        in_AA = True,
        out_AA = True):
    """Will return all residues in pdbfile, residues in same chain will belong
    to the same chain type"""
    pat = re.compile(r'^ATOM|^HETATM')
    pat_meta = re.compile(r'^TITLE|^REMARK|^CRYST1|^MODEL')
    text = [f for f in _string.split('\n') if pat.match(f)]
    meta_text = [f for f in _string.split('\n') if pat_meta.match(f)]

    atoms = []
    res_ids = []
    chain_ids = []
    chain_dict = {}
    for i, line in enumerate( text ):
        res_name = text[i][17:21].strip()
        res_id = int( text[i][22:26].strip() )
        pdb_name = text[i][11:16].strip()
        element = pdb_name[0]
        chain_id = text[i][21:22].strip()
        if chain_id == "":
            chain_id = "X"
        x = text[i][30:38].strip()
        y = text[i][38:46].strip()
        z = text[i][46:54].strip()
        x, y, z = map( float, [x, y, z] )

        atoms.append( NewAtom( x = x, y = y, z = z,
            element = element,
            pdb_name = pdb_name,
            res_id = res_id,
            chain_id = chain_id,
            AA = in_AA,
            res_name = res_name,
            ))

        res_ids.append( res_id )
        chain_ids.append( chain_id )
        if chain_id in chain_dict:
            chain_dict[chain_id].append( res_id )
        else:
            chain_dict[chain_id] = []

    res_ids = utilz.unique( res_ids )
    chain_ids = utilz.unique( chain_ids )

    res = ([NewResidue([a for a in atoms if (a.res_id == r and a.chain_id == c)], AA = in_AA) for c in chain_ids for r in res_ids if r in chain_dict[c] ])

    for each in res:
        each.res_id = each[0].res_id
        each.res_name = each[0].res_name
    

    if in_AA and not out_AA:
        for each in res:
            each.to_AU()

    return res, meta_text
#
def uniq( inp ):
    output = []
    rest = []
    for x in inp:
        if x not in output:
            output.append(x)
        else:
            rest.append(x)
    return output

def almost_eq( a, b, thr = 1.0e-4 ):
    if abs( a - b ) < thr:
        return True
    return False





def get_rep_2( con_list, rep_list, pat_rep, pat_con ):
    tmp_atomlist = []
    for at in rep_list:
        if pat_rep.match( at.pdb_name ):
            tmp_atom = at.copy()
            tmp_atom.label += "-XH"
            tmp_atom.element = "H"
            H = tmp_atom.r
            for at2 in con_list:
                if at2.pdb_name == pat_con[at.pdb_name][0]:
                    real_at = at2.copy()
            C = real_at
            dist = Pattern().pdb_dist_dict[ C.pdb_name ]
            C = C.r
            tmp_atom.x, tmp_atom.y, tmp_atom.z = C + (H - C) * dist / np.linalg.norm(H - C)
            tmp_atomlist.append( tmp_atom )
    return tmp_atomlist

def get_replacement( atom_list, pat_rep, pat_con, same = True ):
    tmp_atomlist = []
    for at in atom_list:
        if pat_rep.match( at.pdb_name ):
            tmp_atom = at.copy()
            tmp_atom.label += "-XH"
            tmp_atom.element = "H"
            H = tmp_atom.r
            for at2 in at.Molecule:
                if at2.pdb_name == pat_con[at.pdb_name][0]:
                    real_at = at2.copy()
            C = real_at
            dist = Pattern().pdb_dist_dict[ C.pdb_name ]
            C = C.r
            tmp_atom.x, tmp_atom.y, tmp_atom.z = C + (H - C) * dist / np.linalg.norm(H - C)
            tmp_atomlist.append( tmp_atom )
    return tmp_atomlist

def get_matching( atom_list, pat_add ):
    tmp = []
    for at in atom_list:
        if pat_add.match( at.pdb_name ):
            tmp.append( at.copy() )
    return tmp
def run_argparse( args ):
    A = argparse.ArgumentParser()

########################################
#
# MOL FILES BASIS SET RELATED
# 
########################################

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

#Defaults to 1. 2 and 3 not done yet
    A.add_argument('-level', type = int, default = 1,
            choices = [1,2,3],
           help = "Level of concap, default = 2")
    A.add_argument('-one_mol', type = int)

    A.add_argument('-qmmm_xyz', default = False, action = 'store_true' )

    A.add_argument('-max_l', type = int, default = 0)
    A.add_argument('-co', '--cutoff', dest = 'cutoff',
            type = float, default = 5.0, help= \
            "Defines cutoff radius for residues to include near QM region, supply -qm [ res_id, [ res_id , [..]] for qm region")
    A.add_argument('-qm', default = [107], type = int, nargs='*')
    A.add_argument('-v' , '-verbose', dest = "verbose", action = 'store_true', default = False)
    A.add_argument('--pol', type = int, default = 1)
    A.add_argument('-write_xyz', nargs = '*', type = str)
    A.add_argument('-potskip', default = False, action = 'store_true' )
    args = A.parse_args( args[1:] )

    if not args.pdbfile:
        print 'you need to add "-p <pdbfile>" after invoking pdbreader'
        raise SystemExit
    return args

class Pattern( dict ):
    """
        Patterns determining which atoms to include from p = previous, n = next, b = bridge
        t = this,

        init argument is an integer determining level of capping, (1, 2, 3 ) implemented (see below)

        XH is the label of a hydrogen which has replaced previously heavy atom, and scaled
        with distance according to the heavy atom types distance in self.pdb_dist_dict

        0: Cap with single hydrogens (not impl.)

        1: Include COH - res - NH2 groups from side residues

        2: Include CH3-COH - res - NH-CH3  groups from side residues

        3: Include H2N-CA-COH - res - NH-CA-COH  groups from side residues
                      /  \         /  \
                    HA   XH       HA1 HA2

        4: Include H2N-CA-COH - res - NH-CA-COH  groups from side residues + side residue chains
                      /  \              / \
                    HA   CB-HB1       HA1 HA2
                        / \
                       HB2 XH 

        Abbreviations :

        cus ; Customly defined residues such as gfp CRO
        pro ; Proline
        reg ; All other normal residues
    """

    #def null_pattern():
    #    p = Pattern()
    #    for a, b, c, d, e in itertools.product( ['reg', 'cus', 'pro' ],
    #            ['res', 'con'],
    #            ['p', 't', 'n', 'pp', 'tt', 'nn', 'pp_p', 'nn_n'],
    #            ['add', 'rep', 'con'],
    #            [ 1, 2, 3]):
# Ne#gative lookahead that will never match http://stackoverflow.com/questions/1723182/a-regex-that-will-never-be-matched-by-anything
    #        p[ (a, b, c, d, e) ] = re.compile(r'(?!x)x')
    #    return p

    def __init__(self):
        self.pat_xyz = re.compile(r'^\s*(\w|-)+\s+(-*\d*.+\d+)\s+(-*\d*.+\d+)\s+(-*\d*.+\d+) *$')
        self.pdb_dist_dict = { "CA" : 1.0, "CD" : 1.0, 
                "SG" : 1.0, "N" : 1.0, "C": 1.0,
                "N1" : 1.0, "C3" : 1.0, "CA1" : 1.0, "CA3" : 1.0,
                "CB" : 1.0, "C1": 1.0 , "N3" : 1.0,
                "CB1" :1.0, "CB2" :1.0,
                }

#Residues level 1 to add
        self[ ( 'reg', 'res', 'p', 'add', 1 ) ] = re.compile(r'C$|O$')
        self[ ( 'cus', 'res', 'p', 'add', 1 ) ] = re.compile(r'C3$|O3$')
        self[ ( 'pro', 'res', 'p', 'add', 1 ) ] = re.compile(r'C$|O$')

        self[ ( 'reg', 'res', 't', 'add', 1 ) ] = self.reg_t = re.compile(r'.*')
        self[ ( 'cus', 'res', 't', 'add', 1 ) ] = self.cus_t = re.compile(r'.*')
        self[ ( 'pro', 'res', 't', 'add', 1 ) ] = self.pro_t = re.compile(r'.*')

        self[ ( 'reg', 'res', 'n', 'add', 1 ) ] = self.reg_n = re.compile(r'N$|HN$|H$')
        self[ ( 'cus', 'res', 'n', 'add', 1 ) ] = self.cus_n = re.compile(r'N1$|HN$')
        self[ ( 'pro', 'res', 'n', 'add', 1 ) ] = self.pro_n = re.compile(r'N$')

        self[ ( 'reg', 'res', 'b', 'add', 1 ) ]= re.compile(r'SG$')

#Residues level 2 to add
        self[ ( 'reg', 'res', 'p', 'add', 2 ) ] = re.compile(r'HA1$|HA2$|HA$|CA$|C$|O$')
        self[ ( 'cus', 'res', 'p', 'add', 2 ) ] = re.compile(r'C3$|O3$|CA3$|HA31$|HA32$')
        self[ ( 'pro', 'res', 'p', 'add', 2 ) ] = re.compile(r'HA1$|HA2$|HA$|CA$|C$|O$')
        self[ ( 'reg', 'res', 't', 'add', 2 ) ] = re.compile(r'.*')
        self[ ( 'cus', 'res', 't', 'add', 2 ) ] = re.compile(r'.*')
        self[ ( 'pro', 'res', 't', 'add', 2 ) ] = re.compile(r'.*')
        self[ ( 'reg', 'res', 'n', 'add', 2 ) ] = re.compile(r'N$|H$|HN$|CA$|HA$|HA1$|HA2$')
        self[ ( 'cus', 'res', 'n', 'add', 2 ) ] = re.compile(r'N1$|HN$|CA1$|HA1$|HA2$')
        self[ ( 'pro', 'res', 'n', 'add', 2 ) ] = re.compile(r'N$|H$|HN$|CA$|HA$|HA1$|HA2$')
        self[ ( 'reg', 'res', 'b', 'add', 2 ) ]= re.compile(r'SG$|CB$|HB1$|HB2$')

#Residues level 3 to add
        self[ ( 'reg', 'res', 'p', 'add', 3 ) ] = re.compile(r'N$|HN$|H$|CB$|HB1$|HB2$|HB$|HA1$|HA2$|HA$|CA$|C$|O$')
        #self[ ( 'cus', 'res', 'p', 'add', 3 ) ] = re.compile(r'C3$|O3$|CA3$|HA31$|HA32$')
        self[ ( 'pro', 'res', 'p', 'add', 3 ) ] = re.compile(r'N$|H1$|H2$|HN$|H$|CB$|HB1$|HB2$|HB$|HA1$|HA2$|HA$|CA$|C$|O$')
        self[ ( 'reg', 'res', 't', 'add', 3 ) ] = re.compile(r'.*')
        self[ ( 'cus', 'res', 't', 'add', 3 ) ] = re.compile(r'.*')
        self[ ( 'pro', 'res', 't', 'add', 3 ) ] = re.compile(r'.*')
        self[ ( 'reg', 'res', 'n', 'add', 3 ) ] = re.compile(r'OC1$|OC2$|CB$|HB1$|HB2$|C$|O$|N$|H$|HN$|CA$|HA$|HA1$|HA2$')
        #self[ ( 'cus', 'res', 'n', 'add', 3 ) ] = re.compile(r'N1$|HN$|CA1$|HA1$|HA2$')
        self[ ( 'pro', 'res', 'n', 'add', 3 ) ] = re.compile(r'OC1$|OC2$|C$|O$|CB$|HB1$|HB2$|N$|H$|HN$|CA$|HA$|HA1$|HA2$')
        self[ ( 'reg', 'res', 'b', 'add', 3 ) ]= re.compile(r'SG$|CB$|HB1$|HB2$')




#Concaps level 1 to add
        self[ ( 'reg', 'con', 'p', 'add', 1 ) ] = re.compile(r'DUMMY')
        self[ ( 'cus', 'con', 'p', 'add', 1 ) ] = re.compile(r'DUMMY')
        self[ ( 'pro', 'con', 'p', 'add', 1 ) ] = re.compile(r'DUMMY')
        self[ ( 'reg', 'con', 't', 'add', 1 ) ] = re.compile(r'C$|O$')
        self[ ( 'cus', 'con', 't', 'add', 1 ) ] = re.compile(r'C3$|O3$')
        self[ ( 'pro', 'con', 't', 'add', 1 ) ] = re.compile(r'C$|O$')
        self[ ( 'reg', 'con', 'n', 'add', 1 ) ] = re.compile(r'N$|HN$|H$')
        self[ ( 'cus', 'con', 'n', 'add', 1 ) ] = re.compile(r'N1$|HN$')
        self[ ( 'pro', 'con', 'n', 'add', 1 ) ] = re.compile(r'N$')
        self[ ( 'reg', 'con', 'b', 'add', 1 ) ]= re.compile(r'SG$')

#Concaps level 2 to add
        self[ ( 'reg', 'con', 'p', 'add', 2 ) ] = re.compile(r'DUMMY')
        self[ ( 'cus', 'con', 'p', 'add', 2 ) ] = re.compile(r'DUMMY')
        self[ ( 'pro', 'con', 'p', 'add', 2 ) ] = re.compile(r'DUMMY')
        self[ ( 'reg', 'con', 't', 'add', 2 ) ] = re.compile(r'HA1$|HA2$|HA$|CA$|C$|O$')
        self[ ( 'cus', 'con', 't', 'add', 2 ) ] = re.compile(r'C3$|O3$|CA3$|HA31$|HA32$')
        self[ ( 'pro', 'con', 't', 'add', 2 ) ] = re.compile(r'HA1$|HA2$|HA$|CA$|C$|O$')
        self[ ( 'reg', 'con', 'n', 'add', 2 ) ] =  re.compile(r'N$|H$|HN$|CA$|HA$|HA1$|HA2$')
        self[ ( 'cus', 'con', 'n', 'add', 2 ) ] =  re.compile(r'N1$|HN$|CA1$|HA1$|HA2$')
        self[ ( 'pro', 'con', 'n', 'add', 2 ) ] =  re.compile(r'N$|H$|CA$|HA$|HA1$|HA2$')
        self[ ( 'reg', 'con', 'b', 'add', 2 ) ] =  re.compile(r'SG$|CB$|HB1$|HB2$')

#Concaps level 3 to add
        self[ ( 'reg', 'con', 'p', 'add', 3 ) ] = re.compile(r'DUMMY')
        self[ ( 'cus', 'con', 'p', 'add', 3 ) ] = re.compile(r'DUMMY')
        self[ ( 'pro', 'con', 'p', 'add', 3 ) ] = re.compile(r'DUMMY')
        self[ ( 'reg', 'con', 't', 'add', 3 ) ] = re.compile(r'CA$|HA1$|HA2$|HA$|N$|HN$|H$|CB$|HB1$|HB2$|C$|O$')
        #self[ ( 'cus', 'con', 't', 'add',31 ) ] = re.compile(r'C3$|O3$')
        self[ ( 'pro', 'con', 't', 'add', 3 ) ] = re.compile(r'H1$|H2$|CA$|HA1$|HA2$|HA$|N$|HN$|H$|CB$|HB1$|HB2$|C$|O$')
        self[ ( 'reg', 'con', 'n', 'add', 3 ) ] = re.compile(r'C$|O$|CA$|HA$|HA1$|HA2$|N$|HN$|H$')
        self[ ( 'cus', 'con', 'n', 'add', 3 ) ] = re.compile(r'OC1$|OC2$|N1$|HN$')
        self[ ( 'pro', 'con', 'n', 'add', 3 ) ] = re.compile(r'OC1$|OC2$|C$|O$|CB$|HB1$|HB2$|N$|CA$|HA$|HA1$|HA2$')
        self[ ( 'reg', 'con', 'b', 'add', 3 ) ]= re.compile(r'SG$')



#Which heavy atoms to replace with hydrogens, the key is to be replaced,
#the list corresponds to various sites that should get a hydrogen connected to the removed
#p_p denotes that both are situated in previous residue

#Residues level 1
        self[ ( 'reg', 'res', 'pp', 'con', 1 ) ] = { "CA" : ["C"] }
        self[ ( 'cus', 'res', 'pp', 'con', 1 ) ] = { "CA3" : ["C3"] }
        self[ ( 'pro', 'res', 'pp', 'con', 1 ) ] = { "CA" : ["C"] }
        self[ ( 'reg', 'res', 'tt', 'con', 1 ) ] = { }
        self[ ( 'cus', 'res', 'tt', 'con', 1 ) ] = { }
        self[ ( 'pro', 'res', 'tt', 'con', 1 ) ] = { }
        self[ ( 'reg', 'res', 'nn', 'con', 1 ) ] = { "CA" : ["N"] }
        self[ ( 'cus', 'res', 'nn', 'con', 1 ) ] = { "CA1" : ["N1"] }
        self[ ( 'pro', 'res', 'nn', 'con', 1 ) ] = { "CA" : ["N"] ,  "CD" : ["N"], "CD2" : ["N"] }
        self[ ( 'reg', 'res', 'bb', 'con', 1 ) ] = { "CB" : ["SG"] }
        self[ ( 'reg', 'res', 'nn_n', 'con', 1 ) ] = {}

        #Dummy patterns
        self[ ( 'pro', 'res', 'nn_n', 'con', 1 ) ] = { }
        self[ ( 'pro', 'res', 'p_t', 'con',  1 ) ] = { }
        self[ ( 'pro', 'res', 'pp_p', 'con', 1 ) ] = { }
        self[ ( 'reg', 'res', 'nn_n', 'con', 1 ) ] = { }
        self[ ( 'reg', 'res', 'p_t', 'con',  1 ) ] = { }
        self[ ( 'reg', 'res', 'pp_p', 'con', 1 ) ] = { }

#Residues level 2
        self[ ( 'reg', 'res', 'pp', 'con', 2 ) ] = { "CB": ["CA"], "N" : ["CA"] }
        self[ ( 'cus', 'res', 'pp', 'con', 2 ) ] = { "N3" : ["CA3"] }
        self[ ( 'pro', 'res', 'pp', 'con', 2 ) ] = { "CB" : ["CA"], "N" : ["CA"] }
        self[ ( 'reg', 'res', 'tt', 'con', 2 ) ] = { }
        self[ ( 'cus', 'res', 'tt', 'con', 2 ) ] = { }
        self[ ( 'pro', 'res', 'tt', 'con', 2 ) ] = { }
        self[ ( 'reg', 'res', 'nn', 'con', 2 ) ] = { "CB" : ["CA"], "C" : ["CA"]  }
        self[ ( 'cus', 'res', 'nn', 'con', 2 ) ] = { "CB1" : ["CA1"], "C1" : ["CA1"] }
        self[ ( 'pro', 'res', 'nn', 'con', 2 ) ] = { "C" : ["CA"] ,  "CB" : ["CA"], "CD" : ["N"], "CD2" : ["N"] }
        self[ ( 'reg', 'res', 'bb', 'con', 2 ) ] = { "CA" : ["CB"] }

        #Dummy patterns
        self[ ( 'pro', 'res', 'nn_n', 'con', 2 ) ] = { }
        self[ ( 'pro', 'res', 'p_t', 'con', 2 ) ] = { }
        self[ ( 'pro', 'res', 'pp_p', 'con', 2 ) ] = { }
        self[ ( 'reg', 'res', 'nn_n', 'con', 2 ) ] = { }
        self[ ( 'reg', 'res', 'p_t', 'con', 2 ) ] = { }
        self[ ( 'reg', 'res', 'pp_p', 'con', 2 ) ] = { }

#Residues level 3
        self[ ( 'reg', 'res', 'pp', 'con', 3 ) ] = { "CG": ["CB"], "SG" : ["CB"],
                "OG" : ["CB"] }
        self[ ( 'cus', 'res', 'pp', 'con', 3 ) ] = { "N3" : ["CA3"] }
        self[ ( 'pro', 'res', 'pp', 'con', 3 ) ] = { "CG" : ["CB"], "CD2" : ["N"],
                "CD" : ["N"] }

        self[ ( 'reg', 'res', 'p_t', 'con', 3 ) ] = { }
        self[ ( 'pro', 'res', 'p_t', 'con', 3 ) ] = { }

        self[ ( 'reg', 'res', 'pp_p', 'con', 3 ) ] = { "C": ["N"] }
        self[ ( 'pro', 'res', 'pp_p', 'con', 3 ) ] = { "C": ["N"] }

        self[ ( 'reg', 'res', 'tt', 'con', 3 ) ] = { }
        self[ ( 'cus', 'res', 'tt', 'con', 3 ) ] = { }
        self[ ( 'pro', 'res', 'tt', 'con', 3 ) ] = { }
        self[ ( 'reg', 'res', 'nn', 'con', 3 ) ] = { "CG" : ["CB"], "SG" : ["CB"],
                "OG" : ["CB"] }
        self[ ( 'cus', 'res', 'nn', 'con', 3 ) ] = { "CB1" : ["CA1"], "C1" : ["CA1"] }
        self[ ( 'pro', 'res', 'nn', 'con', 3 ) ] = { "CG" : ["CB"],
                "CD" : ["N"], "CD2" : ["N"] }

        self[ ( 'reg', 'res', 'nn_n', 'con', 3 ) ] = { "N" : ["C"] }
        self[ ( 'pro', 'res', 'nn_n', 'con', 3 ) ] = { "N" : ["C"] }

        self[ ( 'reg', 'res', 'bb', 'con', 3 ) ] = { "CA" : ["CB"] }



#Concap level 1
        self[ ( 'reg', 'con', 'pp', 'con', 1 ) ] = { "CA" : ["C"] }
        self[ ( 'cus', 'con', 'pp', 'con', 1 ) ] = { "CA3" : ["C3"] }
        self[ ( 'pro', 'con', 'pp', 'con', 1 ) ] = { "CA" : ["C"] }
        self[ ( 'reg', 'con', 'tt', 'con', 1 ) ] = { "CA" : ["C"] }
        self[ ( 'cus', 'con', 'tt', 'con', 1 ) ] = { "CA3" : ["C3"] }
        self[ ( 'pro', 'con', 'tt', 'con', 1 ) ] = { "CA" : ["C"] }
        self[ ( 'reg', 'con', 'nn', 'con', 1 ) ] = { "CA" : ["N"] }
        self[ ( 'cus', 'con', 'nn', 'con', 1 ) ] = { "CA1" : ["N1"] }
        self[ ( 'pro', 'con', 'nn', 'con', 1 ) ] = { "CA" : ["N"] ,  "CD" : ["N"], "CD2" : ["N"] }
        self[ ( 'reg', 'con', 'bb', 'con', 1 ) ] = { "CB" : ["SG"] }
        self[ ( 'reg', 'con', 'nn_n', 'con', 1 ) ] = { }

        #Extra dummy patterns
        self[ ( 'pro', 'con', 'nn_n', 'con', 1 ) ] = { }
        self[ ( 'pro', 'con', 'p_t', 'con', 1 ) ] = { }
        self[ ( 'reg', 'con', 'nn_n', 'con', 1 ) ] = { }
        self[ ( 'reg', 'con', 'p_t', 'con', 1 ) ] = { }

#Concap level 2
        self[ ( 'reg', 'con', 'pp', 'con', 2 ) ] = { "CB": ["CA"], "N" : ["CA"] }
        self[ ( 'cus', 'con', 'pp', 'con', 2 ) ] = { "N3" : ["CA3"] }
        self[ ( 'pro', 'con', 'pp', 'con', 2 ) ] = { "CB" : ["CA"], "N" : ["CA"] }
        self[ ( 'reg', 'con', 'tt', 'con', 2 ) ] = { "CB": ["CA"], "N" : ["CA"] }
        self[ ( 'cus', 'con', 'tt', 'con', 2 ) ] = { "N3" : ["CA3"] }
        self[ ( 'pro', 'con', 'tt', 'con', 2 ) ] = { "CB" : ["CA"], "N" : ["CA"] }
        self[ ( 'reg', 'con', 'nn', 'con', 2 ) ] = { "CB" : ["CA"], "C" : ["CA"]  }
        self[ ( 'cus', 'con', 'nn', 'con', 2 ) ] = { "CB1" : ["CA1"], "C1" : ["CA1"] }
        self[ ( 'pro', 'con', 'nn', 'con', 2 ) ] = { "C" : ["CA"] ,  "CB" : ["CA"], "CD" : ["N"], "CD2" : ["N"] }
        self[ ( 'reg', 'con', 'bb', 'con', 2 ) ] = { "CA" : ["CB"] }

        #Extra dummy patterns
        self[ ( 'pro', 'con', 'nn_n', 'con', 2 ) ] = { }
        self[ ( 'pro', 'con', 'p_t', 'con', 2 ) ] = { }
        self[ ( 'reg', 'con', 'nn_n', 'con', 2 ) ] = { }
        self[ ( 'reg', 'con', 'p_t', 'con', 2 ) ] = { }

#Concap level 3
        self[ ( 'reg', 'con', 'pp', 'con', 3 ) ] = { }
        self[ ( 'cus', 'con', 'pp', 'con', 3 ) ] = { }
        self[ ( 'pro', 'con', 'pp', 'con', 3 ) ] = { }

        self[ ( 'reg', 'con', 'p_t', 'con', 3 ) ] = { "C": ["N"] }
        self[ ( 'pro', 'con', 'p_t', 'con', 3 ) ] = { "C": ["N"] }

        self[ ( 'reg', 'con', 'tt', 'con', 3 ) ] = { "CG" : ["CB"] }
        self[ ( 'cus', 'con', 'tt', 'con', 3 ) ] = { "N3" : ["CA3"] }
        self[ ( 'pro', 'con', 'tt', 'con', 3 ) ] = { "CG" : ["CB"], "CD" : ["N"],
                "CD2" : ["N"] }
        self[ ( 'reg', 'con', 'nn', 'con', 3 ) ] = { "CG" : ["CB"]  }
        self[ ( 'cus', 'con', 'nn', 'con', 3 ) ] = { "CB1" : ["CA1"], "C1" : ["CA1"] }
        self[ ( 'pro', 'con', 'nn', 'con', 3 ) ] = { "CG" : ["CB"], "CD" : ["N"], "CD2" : ["N"] }
        self[ ( 'reg', 'con', 'nn_n', 'con', 3 ) ] = { "N" : ["C"] }
        self[ ( 'pro', 'con', 'nn_n', 'con', 3 ) ] = { "N" : ["C"] }

        self[ ( 'reg', 'con', 'bb', 'con', 3 ) ] = { "CA" : ["CB"] }



#Custom heavy to transform to hydrogen and scale them in distance
#Residues level 1
        self[ ( 'reg', 'res', 'pp', 'rep', 1 ) ]  = re.compile(r'CA$')
        self[ ( 'cus', 'res', 'pp', 'rep', 1 ) ]  = re.compile(r'CA3$')
        self[ ( 'pro', 'res', 'pp', 'rep', 1 ) ]  = re.compile(r'CA$')
        self[ ( 'reg', 'res', 'tt', 'rep', 1 ) ] = re.compile(r'BALLONY')
        self[ ( 'cus', 'res', 'tt', 'rep', 1 ) ] = re.compile(r'BALLONY')
        self[ ( 'pro', 'res', 'tt', 'rep', 1 ) ] = re.compile(r'BALLONY')
        self[ ( 'reg', 'res', 'nn', 'rep', 1 ) ]  = re.compile(r'CA$')
        self[ ( 'cus', 'res', 'nn', 'rep', 1 ) ]  = re.compile(r'CA1$')
        self[ ( 'pro', 'res', 'nn', 'rep', 1 ) ]  = re.compile(r'CA$|CD$|CD2$')
        self[ ( 'reg', 'res', 'bb', 'rep', 1 ) ]  = re.compile(r'CB$')
        self[ ( 'reg', 'res', 'nn_n', 'rep', 1 ) ] = re.compile(r'DUMMY')

        #Dummy patterns
        self[ ( 'pro', 'res', 'nn_n', 'rep', 1 ) ] = re.compile(r'DUMMY')
        self[ ( 'pro', 'res', 'p_t', 'rep', 1 ) ] = re.compile(r'DUMMY')
        self[ ( 'pro', 'res', 'pp_p', 'rep', 1 ) ] = re.compile(r'DUMMY')
        self[ ( 'reg', 'res', 'nn_n', 'rep', 1 ) ] = re.compile(r'DUMMY')
        self[ ( 'reg', 'res', 'p_t', 'rep', 1 ) ] = re.compile(r'DUMMY')
        self[ ( 'reg', 'res', 'pp_p', 'rep', 1 ) ] = re.compile(r'DUMMY')

#Residues level 2
        self[ ( 'reg', 'res', 'pp', 'rep', 2 ) ] = re.compile(r'N$')
        self[ ( 'cus', 'res', 'pp', 'rep', 2 ) ] = re.compile(r'N3$')
        self[ ( 'pro', 'res', 'pp', 'rep', 2 ) ] = re.compile(r'N$|CB$')
        self[ ( 'reg', 'res', 'tt', 'rep', 2 ) ] = re.compile(r'BALLONY')
        self[ ( 'cus', 'res', 'tt', 'rep', 2 ) ] = re.compile(r'BALLONY')
        self[ ( 'pro', 'res', 'tt', 'rep', 2 ) ] = re.compile(r'BALLONY')
        self[ ( 'reg', 'res', 'nn', 'rep', 2 ) ] = re.compile(r'C$')
        self[ ( 'cus', 'res', 'nn', 'rep', 2 ) ] = re.compile(r'C2$|CB1')
        self[ ( 'pro', 'res', 'nn', 'rep', 2 ) ] = re.compile(r'CB$|CD$|CD2$|C$')
        self[ ( 'reg', 'res', 'bb', 'rep', 2 ) ] = re.compile(r'CA$')

        #Dummy patterns
        self[ ( 'pro', 'res', 'nn_n', 'rep', 2 ) ] = re.compile(r'DUMMY')
        self[ ( 'pro', 'res', 'p_t', 'rep', 2 ) ] = re.compile(r'DUMMY')
        self[ ( 'pro', 'res', 'pp_p', 'rep', 2 ) ] = re.compile(r'DUMMY')
        self[ ( 'reg', 'res', 'nn_n', 'rep', 2 ) ] = re.compile(r'DUMMY')
        self[ ( 'reg', 'res', 'p_t', 'rep', 2 ) ] = re.compile(r'DUMMY')
        self[ ( 'reg', 'res', 'pp_p', 'rep', 2 ) ] = re.compile(r'DUMMY')

#Residues level 3
        self[ ( 'reg', 'res', 'pp', 'rep', 3 ) ] = re.compile(r'CG$')
        self[ ( 'cus', 'res', 'pp', 'rep', 3 ) ] = re.compile(r'N3$')
        self[ ( 'pro', 'res', 'pp', 'rep', 3 ) ] = re.compile(r'CG$|CD$')
        self[ ( 'reg', 'res', 'pp_p', 'rep', 3 ) ] = re.compile(r'C$')
        self[ ( 'pro', 'res', 'pp_p', 'rep', 3 ) ] = re.compile(r'C$')


        self[ ( 'reg', 'res', 'tt', 'rep', 3 ) ] = re.compile(r'BALLONY')
        self[ ( 'cus', 'res', 'tt', 'rep', 3 ) ] = re.compile(r'BALLONY')
        self[ ( 'pro', 'res', 'tt', 'rep', 3 ) ] = re.compile(r'BALLONY')
        self[ ( 'reg', 'res', 'nn', 'rep', 3 ) ] = re.compile(r'CG$')
        self[ ( 'cus', 'res', 'nn', 'rep', 3 ) ] = re.compile(r'C2$|CB1')
        self[ ( 'pro', 'res', 'nn', 'rep', 3 ) ] = re.compile(r'CG$|CD$|CD2$')

        self[ ( 'reg', 'res', 'p_t', 'rep', 3 ) ] = re.compile(r'BALLONY')
        self[ ( 'pro', 'res', 'p_t', 'rep', 3 ) ] = re.compile(r'BALLONY')
        self[ ( 'reg', 'res', 'nn_n', 'rep', 3 ) ] = re.compile(r'N$')
        self[ ( 'pro', 'res', 'nn_n', 'rep', 3 ) ] = re.compile(r'N$')

        self[ ( 'reg', 'res', 'bb', 'rep', 3 ) ] = re.compile(r'CA$')

#Concaps level 1
        self[ ( 'reg', 'con', 'pp', 'rep', 1 ) ] = re.compile(r'CA$')
        self[ ( 'cus', 'con', 'pp', 'rep', 1 ) ] = re.compile(r'CA3$')
        self[ ( 'pro', 'con', 'pp', 'rep', 1 ) ] = re.compile(r'CA$')
        self[ ( 'reg', 'con', 'tt', 'rep', 1 ) ] = re.compile(r'CA$')
        self[ ( 'cus', 'con', 'tt', 'rep', 1 ) ] = re.compile(r'CA3$')
        self[ ( 'pro', 'con', 'tt', 'rep', 1 ) ] = re.compile(r'CA$')
        self[ ( 'reg', 'con', 'nn', 'rep', 1 ) ] = re.compile(r'CA$')
        self[ ( 'cus', 'con', 'nn', 'rep', 1 ) ] = re.compile(r'CA1$')
        self[ ( 'pro', 'con', 'nn', 'rep', 1 ) ] = re.compile(r'CA$|CD$|CD2$')
        self[ ( 'reg', 'con', 'bb', 'rep', 1 ) ] = re.compile(r'CB$')
        self[ ( 'reg', 'con', 'nn_n', 'rep', 1 ) ] = re.compile(r'DUMMY')

        #Extra dummy patterns
        self[ ( 'pro', 'con', 'nn_n', 'rep', 1 ) ] = re.compile(r'DUMMY')
        self[ ( 'pro', 'con', 'p_t', 'rep', 1 ) ] = re.compile(r'DUMMY')
        self[ ( 'reg', 'con', 'nn_n', 'rep', 1 ) ] = re.compile(r'DUMMY')
        self[ ( 'reg', 'con', 'p_t', 'rep', 1 ) ] = re.compile(r'DUMMY')

#Concaps level 2
        self[ ( 'reg', 'con', 'pp', 'rep', 2 ) ] = re.compile(r'N$')
        self[ ( 'cus', 'con', 'pp', 'rep', 2 ) ] = re.compile(r'N3$')
        self[ ( 'pro', 'con', 'pp', 'rep', 2 ) ] = re.compile(r'N$|CB$')
        self[ ( 'reg', 'con', 'tt', 'rep', 2 ) ] = re.compile(r'N$')
        self[ ( 'cus', 'con', 'tt', 'rep', 2 ) ] = re.compile(r'N3$')
        self[ ( 'pro', 'con', 'tt', 'rep', 2 ) ] = re.compile(r'N$|CB$')
        self[ ( 'reg', 'con', 'nn', 'rep', 2 ) ] = re.compile(r'C$')
        self[ ( 'cus', 'con', 'nn', 'rep', 2 ) ] = re.compile(r'C2$|CB1')
        self[ ( 'pro', 'con', 'nn', 'rep', 2 ) ] = re.compile(r'CB$|CD$|CD2$|C$')
        self[ ( 'reg', 'con', 'bb', 'rep', 2 ) ] = re.compile(r'CA$')

        #Extra dummy patterns
        self[ ( 'pro', 'con', 'nn_n', 'rep', 2 ) ] = re.compile(r'DUMMY')
        self[ ( 'pro', 'con', 'p_t', 'rep', 2 ) ] = re.compile(r'DUMMY')
        self[ ( 'reg', 'con', 'nn_n', 'rep', 2 ) ] = re.compile(r'DUMMY')
        self[ ( 'reg', 'con', 'p_t', 'rep', 2 ) ] = re.compile(r'DUMMY')

#Concaps level 3
        self[ ( 'reg', 'con', 'pp', 'rep', 3 ) ] = re.compile(r'N$')
        self[ ( 'cus', 'con', 'pp', 'rep', 3 ) ] = re.compile(r'N3$')
        self[ ( 'pro', 'con', 'pp', 'rep', 3 ) ] = re.compile(r'N$|CB$')

        self[ ( 'reg', 'con', 'p_t', 'rep', 3 ) ] = re.compile(r'C$')
        self[ ( 'pro', 'con', 'p_t', 'rep', 3 ) ] = re.compile(r'C$')

        self[ ( 'reg', 'con', 'tt', 'rep', 3 ) ] = re.compile(r'CG$')
        self[ ( 'cus', 'con', 'tt', 'rep', 3 ) ] = re.compile(r'N3$')
        self[ ( 'pro', 'con', 'tt', 'rep', 3 ) ] = re.compile(r'CD$|CD2$|CG$')
        self[ ( 'reg', 'con', 'nn', 'rep', 3 ) ] = re.compile(r'CG$')
        self[ ( 'cus', 'con', 'nn', 'rep', 3 ) ] = re.compile(r'C2$|CB1')
        self[ ( 'pro', 'con', 'nn', 'rep', 3 ) ] = re.compile(r'CG$|CD$|CD2$')

        self[ ( 'reg', 'con', 'nn_n', 'rep', 3 ) ] = re.compile(r'N$')
        self[ ( 'pro', 'con', 'nn_n', 'rep', 3 ) ] = re.compile(r'N$')

        self[ ( 'reg', 'con', 'bb', 'rep', 3 ) ] = re.compile(r'CA$')



    def get( self, res_name ="reg", res_type = "res",  what_pos ="t", what_todo = 'add', level = 1 ):
        return self[ (res_name, res_type, what_pos, what_todo,level) ]

class NewAtom( molecules.Atom ):
    def __init__(self, *args, **kwargs):
        self._chain_id = None
        self._name = None
        self._res_name = None
        super( NewAtom, self ).__init__( *args, **kwargs )

        if kwargs != {}:
            setattr( self, "_chain_id",  kwargs.get( "chain_id", None ) )
            setattr( self, "_res_name",  kwargs.get( "res_name", None ) )

    @property
    def res_name(self):
        if self._res_name:
            return self._res_name
        return "XXX"
    @res_name.setter
    def res_name(self,val):
        self._res_name = val

    @property
    def chain_id(self):
        if self._chain_id:
            return self._chain_id
        if self.Molecule:
            if self.Molecule.Chain:
                return self.Molecule.Chain.chain_id
        return None

class Atom( molecules.Atom ):

    def __init__(self, *args, **kwargs):
        super( Atom, self ).__init__( *args, **kwargs )
        self.pdb_name = None

        self._label = None
        self.residue = None

        self._res_name = None

        self.in_qm_region = False
        self.in_qmmm_border = False

    @property
    def order_nr(self):
        return self.residue.index( self )

    def dist_to_atom(self, other):
        r = np.sqrt( (float(self.x) - float(other.x))**2 +(float(self.y) - float(other.y))**2 +(float(self.z) - float(other.z))**2)
        return  r

    def copy( self):
        """Copy Atom method"""
        new = Atom()
        new.x, new.y, new.z = self.r.copy()
        new.pdb_name = self.pdb_name
        new.element = self.element
        new._label = self.label
#To keep original res_id info not overiding when adding to new residue
        new._res_id = self._res_id
        new.Property = self.Property.copy_property()
        return new

    def __str__(self):
        return self.pdb_name

    def transfer_props(self, other):
        other.Props += self.Props
        self.Props.set_to_zero()

    def get_closest( self, cutoff = 1.0, residues = 0 ):
        """Given cutoff and residues, returns a list of closest atoms
        to this one.
        
        By default only check this residue, include more if
        atom happends to be on a bordering region between two residues"""


        alist = [a for a in self.residue if a is not self and 
                self.dist_to_atom(a) < cutoff ]
        alist.sort( key = lambda x: x.dist_to_atom(self ) )

        return alist

    def setXyz(self, array):
        self.x = array[0]
        self.y = array[1]
        self.z = array[2]

    def transfer_own_props_to_list(self, atm_list):
        """ given a list atm_list, will evenly distribute properties
        in self.Props to each atom in atm_list"""

        self.Props /= len(atm_list)

        for atm in atm_list:
            atm.Props += self.Props

    def xyz(self):
        return " ".join( [self.element] + map( str, self.r() )  )

class NewResidue( molecules.Molecule ):
    def __init__(self, *args, **kwargs):
        self._chain_id = None
        self._snapshot = None
        self._res_id = None
        self._res_name = None
        self._Chain = None
        self._Cluster = None
        self._Bridge = None
        self._Next = None
        self._Prev = None
        self.is_concap = False
        self.is_ready = False
        self.is_bridge = False

        self.label_dict = {}

        super( NewResidue, self ).__init__( *args, **kwargs )

        self.concap = None
        self.ready = None

        self.c_term = False
        self.n_term = False

    @property
    def Prev(self):
        return self._Prev
    @property
    def Next(self):
        return self._Next

    @property
    def Bridge(self):
        if self._Bridge:
            return self._Bridge
    @Bridge.setter
    def Bridge(self, val):
        self._Bridge = val

    @property
    def snapshot(self):
        if self._snapshot:
            return self._snapshot
        if self.Chain:
            return self.Chain.snapshot

#NewResidue
    def __str__(self):
        base = "-".join( [self.Chain.chain_id, self.res_name + str(self.res_id)] )
        if self.is_concap:
            base += '-con'
        if self.is_ready:
            base += '-ready'
        return base

#NewResidue
    def add_atom( self, atom):
        if isinstance( atom, molecules.Atom ):
            self.append( atom )
            self.label_dict[ atom.label ] = atom
            atom.Molecule = self

#NewResidue
    def get_mol(self, basis = ['ano-1 2','ano-1 4 3 1', 'ano-2 5 4 1' ] ):
        charge = 0

        elem = [x.element for x in self]
        atomtypes = len( uniq( elem ) )

        if "H" in elem:
            h_atoms = [ x for x in self if x.element == "H" ]
        if "C" in elem:
            c_atoms = [ x for x in self if x.element == "C" ]
        if "N" in elem:
            n_atoms = [ x for x in self if x.element == "N" ]
        if "O" in elem:
            o_atoms = [ x for x in self if x.element == "O" ]
        if "S" in elem:
            s_atoms = [ x for x in self if x.element == "S" ]

        if not self.is_concap:
            if self.n_term:
                charge += 1
            elif self.c_term:
                charge -= 1
            if res_dict[ self.res_name ] in chargeDict:
                charge += chargeDict[ res_dict[ self.res_name] ]
#Temporary fix to give level 3 concaps negative and positive even if they are not
#n or c terminal
        if self._level == 3:
            charge = 0
            for at in self:
                if at.pdb_name == "H1":
                    charge += 1
                    break
                if at.pdb_name == "OC1":
                    charge -= 1
                    break
                

        string = ""
        string += 'ATOMBASIS\n'
        string += 'The studied system is %s\n' % ( self ) 
        string += 'File generated by pdbreader 4.0 -- //By Ignat Harczuk --\n'
        string += 'Atomtypes=%s Charge=%s Angstrom Nosymm\n'%( str(atomtypes ) , str(charge) )
        
        if "H" in elem:
            string += 'Charge=1.0 Atoms=%d Basis=%s\n'%( len(h_atoms), basis[0])
            for k in h_atoms:
                string += '{0:15s}{1:10.5f}{2:10.5f}{3:10.5f}\n'.format( *tuple([ k.label, float(k.x), float(k.y), float(k.z) ]) )
        if "C" in elem:
            string += 'Charge=6.0 Atoms=%d Basis=%s\n'%( len(c_atoms), basis[1])
            for k in c_atoms:
                string += '{0:15s}{1:10.5f}{2:10.5f}{3:10.5f}\n'.format( *tuple([ k.label, float(k.x), float(k.y), float(k.z) ]) )
        if "N" in elem:
            string += 'Charge=7.0 Atoms=%d Basis=%s\n'%( len(n_atoms), basis[1])
            for k in n_atoms:
                string += '{0:15s}{1:10.5f}{2:10.5f}{3:10.5f}\n'.format( *tuple([ k.label, float(k.x), float(k.y), float(k.z) ]) )
        if "O" in elem:
            string += 'Charge=8.0 Atoms=%d Basis=%s\n'%( len(o_atoms), basis[1])
            for k in o_atoms:
                string += '{0:15s}{1:10.5f}{2:10.5f}{3:10.5f}\n'.format( *tuple([ k.label, float(k.x), float(k.y), float(k.z) ]) )
        if "S" in elem:
            string += 'Charge=16.0 Atoms=%d Basis=%s\n'%( len(s_atoms), basis[2])
            for k in s_atoms:
                string += '{0:15s}{1:10.5f}{2:10.5f}{3:10.5f}\n'.format( *tuple([ k.label, float(k.x), float(k.y), float(k.z) ]) )
        return string.rstrip( '\n' )



#Property of NewResidue
    @property
    def Chain(self):
        if self._Chain:
            return self._Chain
        if self.Cluster:
            return self.Cluster
        return None
    @Chain.setter
    def Chain(self, val):
        self._Chain = val

    @property
    def res_id(self):
        if self._res_id:
            return self._res_id
        tmp_id = self[0].res_id
        for at in self:
            try:
                assert tmp_id == at.res_id
            except AssertionError:
                logging.error( "No _res_id in NewResidue and not all atoms in same residue")
        return tmp_id

    def get_atom(self, pdb_name):
        #try:
        for i in self:
            if i.pdb_name == pdb_name:
                return i.copy()

    @res_id.setter
    def res_id(self, val):
        self._res_id = val

    @property
    def res_name(self):
        if self._res_name:
            return self._res_name
        return "YYY"
    @res_name.setter
    def res_name(self, val):
        self._res_name = val





    @property
    def chain_id(self):
        if self.Chain:
            return self.Chain.chain_id
        tmp_ch = self[0].chain_id
        for at in self:
            try:
                assert tmp_ch == at.chain_id
            except AssertionError:
                logging.error( "No Chain object or _chain_id in NewResidue and not all atoms have same chain_id")
        return tmp_ch

#NewResidue Method
    def copy_info(self):
        new = NewResidue()
        new.c_term = self.c_term
        new.n_term = self.n_term
        new.AA = self.AA

        new.in_qm_region = self.in_qm
        new.in_mm_region = self.in_mm

        new._res_id = self.res_id
        new._res_name = self.res_name
        new._Next =   self.Next
        new._Prev =   self.Prev
        new._Bridge = self.Bridge
        new._Chain =  self.Chain
        new._Cluster =  self.Cluster
        return  new



#NewResidue Method
    def gather_ready( self, 
            residue = False, r = False,
            concap = False, c = False,
            bridge = False, b = False,
            level = 1 ):
        p = Pattern()
        if residue or r:
            tmp_residue = self.copy_info()
            tmp_residue.is_ready = True
            res_type = "res"
        elif concap or c:
            tmp_residue = self.copy_info()
            tmp_residue = self.copy_info() 
            tmp_residue.is_concap = True
            res_type = "con"
        elif bridge or b:
            tmp_residue = self.copy_info() 
            tmp_residue = Residue()
            tmp_residue.AA = self.AA
            tmp_residue.is_bridge = True
            res_type = "bri"
        else:
            return
        tmp_residue._level = level

        p_rep_pp_p = p.get()
        
#To set defaults for pattern for this residue
        if self.res_name in proline_dict:
            res_name = "pro"
        elif self.res_name in custom_dict:
            res_name = "cus"
        else:
            res_name = "reg"

        p_con_pp_p = {}
        p_rep_pp_p = re.compile( r'(?!x)x' )
# Patterns to just add atoms, previos depends on type of previous
        if self.Prev:
            if self.Prev.res_name in proline_dict:
                p_add_p = p.get( res_type = res_type, 
                        res_name = "pro",
                        level = level,
                        what_pos ='p', 
                        what_todo = 'add')
                p_rep_pp = p.get( res_type = res_type, 
                        res_name = "pro",
                        level = level,
                        what_pos ='pp', 
                        what_todo = 'rep')
                p_con_pp = p.get( res_type = res_type, 
                        res_name = "pro",
                        level = level,
                        what_pos ='pp', 
                        what_todo = 'con')
                p_con_p_t = p.get( res_type = res_type, 
                        res_name = "pro",
                        level = level,
                        what_pos ='p_t', 
                        what_todo = 'con')
                p_rep_p_t = p.get( res_type = res_type, 
                        res_name = "pro",
                        level = level,
                        what_pos ='p_t', 
                        what_todo = 'rep')

                if self.Prev.Prev and not c:
                    p_rep_pp_p = p.get( res_type = res_type, 
                        res_name = "pro",
                        level = level,
                        what_pos ='pp_p', 
                        what_todo = 'rep')
                    p_con_pp_p = p.get( res_type = res_type, 
                        res_name = "pro",
                        level = level,
                        what_pos ='pp_p', 
                        what_todo = 'con')

            elif self.Prev.res_name in custom_dict:
                p_add_p = p.get( res_type = res_type, 
                        res_name = "cus",
                        level = level,
                        what_pos ='p', 
                        what_todo = 'add')
                p_rep_pp = p.get( res_type = res_type, 
                        res_name = "cus",
                        level = level,
                        what_pos ='pp', 
                        what_todo = 'rep')
                p_con_pp = p.get( res_type = res_type, 
                        res_name = "cus",
                        level = level,
                        what_pos ='pp', 
                        what_todo = 'con')
            else:
                p_add_p = p.get( res_type = res_type, 
                        res_name = "reg",
                        level = level,
                        what_pos ='p', 
                        what_todo = 'add')
                p_rep_pp = p.get( res_type = res_type, 
                        res_name = "reg",
                        level = level,
                        what_pos ='pp', 
                        what_todo = 'rep')
                p_con_pp = p.get( res_type = res_type, 
                        res_name = "reg",
                        level = level,
                        what_pos ='pp', 
                        what_todo = 'con')
                p_con_p_t = p.get( res_type = res_type, 
                        res_name = "reg",
                        level = level,
                        what_pos ='p_t', 
                        what_todo = 'con')

                p_rep_p_t = p.get( res_type = res_type, 
                        res_name = "reg",
                        level = level,
                        what_pos ='p_t', 
                        what_todo = 'rep')

                if self.Prev.Prev and not c:
                    p_rep_pp_p = p.get( res_type = res_type, 
                        res_name = "reg",
                        level = level,
                        what_pos ='pp_p', 
                        what_todo = 'rep')
                    p_con_pp_p = p.get( res_type = res_type, 
                        res_name = "reg",
                        level = level,
                        what_pos ='pp_p', 
                        what_todo = 'con')



# Set adding for this residue, should depend on self
        p_add_t = p.get( res_type = res_type, 
                res_name = res_name,
                level = level,
                what_pos ='t', 
                what_todo = 'add')
        if self.Next:
            if self.Next.res_name in proline_dict:
                p_add_n = p.get( res_type = res_type, 
                        res_name = "pro",
                        level = level,
                        what_pos ='n', 
                        what_todo = 'add')
                p_rep_nn = p.get( res_type = res_type, 
                        res_name = "pro",
                        level = level,
                        what_pos ='nn', 
                        what_todo = 'rep')
                p_con_nn = p.get( res_type = res_type, 
                        res_name = "pro",
                        level = level,
                        what_pos ='nn', 
                        what_todo = 'con')

                if self.Next.Next:
                    p_con_nn_n = p.get( res_type = res_type, 
                        res_name = "pro",
                        level = level,
                        what_pos ='nn_n', 
                        what_todo = 'con')
                    p_rep_nn_n = p.get( res_type = res_type, 
                        res_name = "pro",
                        level = level,
                        what_pos ='nn_n', 
                        what_todo = 'rep')

            elif self.Next.res_name in custom_dict:
                p_add_n = p.get( res_type = res_type, 
                        res_name = "cus",
                        level = level,
                        what_pos ='n', 
                        what_todo = 'add')
                p_rep_nn = p.get( res_type = res_type, 
                        res_name = "cus",
                        level = level,
                        what_pos ='nn', 
                        what_todo = 'rep')
                p_con_nn = p.get( res_type = res_type, 
                        res_name = "cus",
                        level = level,
                        what_pos ='nn', 
                        what_todo = 'con')
            else:
                p_add_n = p.get( res_type = res_type, 
                        res_name = "reg",
                        level = level,
                        what_pos ='n', 
                        what_todo = 'add')
                p_rep_nn = p.get( res_type = res_type, 
                        res_name = "reg",
                        level = level,
                        what_pos ='nn', 
                        what_todo = 'rep')
                p_con_nn = p.get( res_type = res_type, 
                        res_name = "reg",
                        level = level,
                        what_pos ='nn', 
                        what_todo = 'con')
                if self.Next.Next:
                    p_con_nn_n = p.get( res_type = res_type, 
                        res_name = "reg",
                        level = level,
                        what_pos ='nn_n', 
                        what_todo = 'con')
# Patterns to replace these atoms, e.g. the CA from previous residue at level = 1
                    p_rep_nn_n = p.get( res_type = res_type, 
                        res_name = "reg",
                        level = level,
                        what_pos ='nn_n', 
                        what_todo = 'rep')
# Patterns to replace these atoms, e.g. the CA from previous residue at level = 1
#
        p_rep_tt = p.get( res_type = res_type, 
                res_name = res_name,
                level = level,
                what_pos ='tt', 
                what_todo = 'rep')
# Patterns to connect, governs which atoms to connect the H that replaced
# The atom given in pattern replace above
#
# e.g. the H that replaced the CA at level 1 should connect to N, and scaled in 
# distance
        p_con_tt = p.get( res_type = res_type, 
                res_name = res_name,
                level = level,
                what_pos ='tt', 
                what_todo = 'con')

        if self.n_term:
            at_r3 = []
            at_t = get_matching( self, p_add_t )
            at_n = get_matching( self.Next, p_add_n )
            at_r1 = get_replacement( self, p_rep_tt, p_con_tt, )

            at_r2 = get_replacement( self.Next, p_rep_nn, p_con_nn, )
            if self.Next.Next:
                at_r3 = get_rep_2( self.Next, self.Next.Next, p_rep_nn_n, p_con_nn_n )
            for at in at_t + at_n + at_r1 + at_r2 + at_r3:
                tmp_residue.add_atom( at )

        elif self.c_term:
            if concap or c:
                return 
            at_r2 = []
            at_p = get_matching( self.Prev, p_add_p )
            at_t = get_matching( self, p_add_t )
            at_r1 = get_replacement( self.Prev, p_rep_pp, p_con_pp, )
            if self.Prev.Prev:
                at_r2 = get_rep_2( self.Prev, self.Prev.Prev, p_rep_pp_p, p_con_pp_p )
            for at in at_p + at_t + at_r1 + at_r2 :
                tmp_residue.add_atom( at )
        else:
            at_r1, at_r2 = [], []
            if concap or c:
                at_r1 = get_replacement( self, p_rep_tt, p_con_tt, )
                at_r2 = get_rep_2( self, self.Prev, p_rep_p_t, p_con_p_t, )
            else:
                at_r1 = get_replacement( self.Prev, p_rep_pp, p_con_pp, )

            at_r4, at_r5 = [], []
            at_p = get_matching( self.Prev, p_add_p )
            at_t = get_matching( self, p_add_t )
            at_n = get_matching( self.Next, p_add_n )
            at_r3 = get_replacement( self.Next, p_rep_nn, p_con_nn, )
            if self.Next.Next:
                at_r4 = get_rep_2( self.Next, self.Next.Next, p_rep_nn_n, p_con_nn_n )
            if self.Prev.Prev:
                at_r5 = get_rep_2( self.Prev, self.Prev.Prev, p_rep_pp_p, p_con_pp_p )

            for at in at_p + at_t + at_n + at_r1 + at_r2 + at_r3 + at_r4 + at_r5:
                tmp_residue.add_atom( at )
        for at in tmp_residue:
            at.Molecule = tmp_residue


        if residue or r:
            self.ready = tmp_residue
        if concap or c:
            self.concap = tmp_residue
        if bridge or b:
            self.bri = tmp_residue

class Residue( molecules.Molecule ):
    """Class designed to calculate MFCC procedure for proteins"""

    def __init__(self, *args, **kwargs):
        self.label_dict = {}
        super( Residue, self ).__init__( *args, **kwargs )
        self.Chain = None
        self.c_term = False
        self.n_term = False
        self.concap = False
        self.in_qm_region = False
        self.in_mm_region = True
        self._res_id = None
        self._res_name = None
        self.Prev = None
        self.Next = None
        self.Bridge = None
        self._level = None

#properties
        self._label = None
# Ready points to a residue which is fully saturated with side atoms
        self.ready = None
# con points to a residue which is fully saturated with side atoms as concap
        self.con = None
        #self.linewidth = {"X":25,"H":25, "N": 30, "C": 30, "O":40, "P" : 40,
        #                        'S' : 45 }
        #self.style = { "X": 'ko' ,"H":'wo', "N":'bo',"C":'go',"P":'ko', "O":'ro',

                                #'S' : 'yo'}
        self.color = { "X": 'black' ,"H":'white', "N":'blue',"C":'green',"P":'black', "O":'red', 'S' : 'yellow'}

    @property
    def N(self):
        return self.get_atom_by_pdbname( 'N' )
    @property
    def O(self):
        return self.get_atom_by_pdbname( 'O' )
    @property
    def OG1(self):
        return self.get_atom_by_pdbname( 'OG1' )
    @property
    def OG2(self):
        return self.get_atom_by_pdbname( 'OG2' )
    @property
    def C(self):
        return self.get_atom_by_pdbname( 'C' )
    @property
    def CA(self):
        return self.get_atom_by_pdbname( 'CA' )
    @property
    def CB(self):
        return self.get_atom_by_pdbname( 'CB' )
    @property
    def CB1(self):
        return self.get_atom_by_pdbname( 'CB1' )
    @property
    def CB2(self):
        return self.get_atom_by_pdbname( 'CB2' )

    @property
    def order_nr(self):
        if self.Cluster:
            return self.Cluster.index( self )
        else:
            return self.res_id
    def center( self, atom = None ):
        if atom is None:
            atom = self[0]
        vec = np.zeros(3) - atom.r
        for at in self:
            at.x, at.y, at.z = vec + at.r

    def rotate(self, r1, r2, t1, t2, t3):
        self.inv_rotate( r1, r2 )
        for at in self:
            r1 = molecules.Rotator.get_Rz( t1 )
            r2 = molecules.Rotator.get_Ry_inv( t2 )
            r3 = molecules.Rotator.get_Rz( t3 )
            at.x, at.y, at.z = reduce(lambda a, x : np.einsum('ij,j',x,a),[r1,r2,r3], at.r)

    def gather_ready( self, 
            residue = False, r = False,
            concap = False, c = False,
            bridge = False, b = False,
            level = 1 ):
        p = Pattern()
        if residue or r:
            tmp_residue = self.copy_info()
            res_type = "res"
        elif concap or c:
            tmp_residue = self.copy_info()
            tmp_residue = self.copy_info() 
            tmp_residue.concap = True
            res_type = "con"
        elif bridge or b:
            tmp_residue = self.copy_info() 
            tmp_residue = Residue()
            tmp_residue.AA = self.AA
            tmp_residue.concap = True
            tmp_residue.bridge = True
            res_type = "bri"
        else:
            return
        tmp_residue._level = level

        p_rep_pp_p = p.get()
        
#To set defaults for pattern for this residue
        if self.res_name in proline_dict:
            res_name = "pro"
        elif self.res_name in custom_dict:
            res_name = "cus"
        else:
            res_name = "reg"

        p_con_pp_p = {}
        p_rep_pp_p = re.compile( r'(?!x)x' )
# Patterns to just add atoms, previos depends on type of previous
        if self.Prev:
            if self.Prev.res_name in proline_dict:
                p_add_p = p.get( res_type = res_type, 
                        res_name = "pro",
                        level = level,
                        what_pos ='p', 
                        what_todo = 'add')
                p_rep_pp = p.get( res_type = res_type, 
                        res_name = "pro",
                        level = level,
                        what_pos ='pp', 
                        what_todo = 'rep')
                p_con_pp = p.get( res_type = res_type, 
                        res_name = "pro",
                        level = level,
                        what_pos ='pp', 
                        what_todo = 'con')
                p_con_p_t = p.get( res_type = res_type, 
                        res_name = "pro",
                        level = level,
                        what_pos ='p_t', 
                        what_todo = 'con')
                p_rep_p_t = p.get( res_type = res_type, 
                        res_name = "pro",
                        level = level,
                        what_pos ='p_t', 
                        what_todo = 'rep')

                if self.Prev.Prev and not c:
                    p_rep_pp_p = p.get( res_type = res_type, 
                        res_name = "pro",
                        level = level,
                        what_pos ='pp_p', 
                        what_todo = 'rep')
                    p_con_pp_p = p.get( res_type = res_type, 
                        res_name = "pro",
                        level = level,
                        what_pos ='pp_p', 
                        what_todo = 'con')

            elif self.Prev.res_name in custom_dict:
                p_add_p = p.get( res_type = res_type, 
                        res_name = "cus",
                        level = level,
                        what_pos ='p', 
                        what_todo = 'add')
                p_rep_pp = p.get( res_type = res_type, 
                        res_name = "cus",
                        level = level,
                        what_pos ='pp', 
                        what_todo = 'rep')
                p_con_pp = p.get( res_type = res_type, 
                        res_name = "cus",
                        level = level,
                        what_pos ='pp', 
                        what_todo = 'con')
            else:
                p_add_p = p.get( res_type = res_type, 
                        res_name = "reg",
                        level = level,
                        what_pos ='p', 
                        what_todo = 'add')
                p_rep_pp = p.get( res_type = res_type, 
                        res_name = "reg",
                        level = level,
                        what_pos ='pp', 
                        what_todo = 'rep')
                p_con_pp = p.get( res_type = res_type, 
                        res_name = "reg",
                        level = level,
                        what_pos ='pp', 
                        what_todo = 'con')
                p_con_p_t = p.get( res_type = res_type, 
                        res_name = "reg",
                        level = level,
                        what_pos ='p_t', 
                        what_todo = 'con')

                p_rep_p_t = p.get( res_type = res_type, 
                        res_name = "reg",
                        level = level,
                        what_pos ='p_t', 
                        what_todo = 'rep')

                if self.Prev.Prev and not c:
                    p_rep_pp_p = p.get( res_type = res_type, 
                        res_name = "reg",
                        level = level,
                        what_pos ='pp_p', 
                        what_todo = 'rep')
                    p_con_pp_p = p.get( res_type = res_type, 
                        res_name = "reg",
                        level = level,
                        what_pos ='pp_p', 
                        what_todo = 'con')



# Set adding for this residue, should depend on self
        p_add_t = p.get( res_type = res_type, 
                res_name = res_name,
                level = level,
                what_pos ='t', 
                what_todo = 'add')
        if self.Next:
            if self.Next.res_name in proline_dict:
                p_add_n = p.get( res_type = res_type, 
                        res_name = "pro",
                        level = level,
                        what_pos ='n', 
                        what_todo = 'add')
                p_rep_nn = p.get( res_type = res_type, 
                        res_name = "pro",
                        level = level,
                        what_pos ='nn', 
                        what_todo = 'rep')
                p_con_nn = p.get( res_type = res_type, 
                        res_name = "pro",
                        level = level,
                        what_pos ='nn', 
                        what_todo = 'con')

                if self.Next.Next:
                    p_con_nn_n = p.get( res_type = res_type, 
                        res_name = "pro",
                        level = level,
                        what_pos ='nn_n', 
                        what_todo = 'con')
                    p_rep_nn_n = p.get( res_type = res_type, 
                        res_name = "pro",
                        level = level,
                        what_pos ='nn_n', 
                        what_todo = 'rep')

            elif self.Next.res_name in custom_dict:
                p_add_n = p.get( res_type = res_type, 
                        res_name = "cus",
                        level = level,
                        what_pos ='n', 
                        what_todo = 'add')
                p_rep_nn = p.get( res_type = res_type, 
                        res_name = "cus",
                        level = level,
                        what_pos ='nn', 
                        what_todo = 'rep')
                p_con_nn = p.get( res_type = res_type, 
                        res_name = "cus",
                        level = level,
                        what_pos ='nn', 
                        what_todo = 'con')
            else:
                p_add_n = p.get( res_type = res_type, 
                        res_name = "reg",
                        level = level,
                        what_pos ='n', 
                        what_todo = 'add')
                p_rep_nn = p.get( res_type = res_type, 
                        res_name = "reg",
                        level = level,
                        what_pos ='nn', 
                        what_todo = 'rep')
                p_con_nn = p.get( res_type = res_type, 
                        res_name = "reg",
                        level = level,
                        what_pos ='nn', 
                        what_todo = 'con')
                if self.Next.Next:
                    p_con_nn_n = p.get( res_type = res_type, 
                        res_name = "reg",
                        level = level,
                        what_pos ='nn_n', 
                        what_todo = 'con')
# Patterns to replace these atoms, e.g. the CA from previous residue at level = 1
                    p_rep_nn_n = p.get( res_type = res_type, 
                        res_name = "reg",
                        level = level,
                        what_pos ='nn_n', 
                        what_todo = 'rep')
# Patterns to replace these atoms, e.g. the CA from previous residue at level = 1
#
        p_rep_tt = p.get( res_type = res_type, 
                res_name = res_name,
                level = level,
                what_pos ='tt', 
                what_todo = 'rep')
# Patterns to connect, governs which atoms to connect the H that replaced
# The atom given in pattern replace above
#
# e.g. the H that replaced the CA at level 1 should connect to N, and scaled in 
# distance
        p_con_tt = p.get( res_type = res_type, 
                res_name = res_name,
                level = level,
                what_pos ='tt', 
                what_todo = 'con')

        if self.n_term:
            at_r3 = []
            at_t = get_matching( self, p_add_t )
            at_n = get_matching( self.Next, p_add_n )
            at_r1 = get_replacement( self, p_rep_tt, p_con_tt, )
            at_r2 = get_replacement( self.Next, p_rep_nn, p_con_nn, )
            if self.Next.Next:
                at_r3 = get_rep_2( self.Next, self.Next.Next, p_rep_nn_n, p_con_nn_n )
            for at in at_t + at_n + at_r1 + at_r2 + at_r3:
                tmp_residue.add_atom( at )

        elif self.c_term:
            if concap or c:
                return 
            at_r2 = []
            at_p = get_matching( self.Prev, p_add_p )
            at_t = get_matching( self, p_add_t )
            at_r1 = get_replacement( self.Prev, p_rep_pp, p_con_pp, )
            if self.Prev.Prev:
                at_r2 = get_rep_2( self.Prev, self.Prev.Prev, p_rep_pp_p, p_con_pp_p )
            for at in at_p + at_t + at_r1 + at_r2 :
                tmp_residue.add_atom( at )
        else:
            at_r1, at_r2 = [], []
            if concap or c:
                at_r1 = get_replacement( self, p_rep_tt, p_con_tt, )
                at_r2 = get_rep_2( self, self.Prev, p_rep_p_t, p_con_p_t, )
            else:
                at_r1 = get_replacement( self.Prev, p_rep_pp, p_con_pp, )

            at_r4, at_r5 = [], []
            at_p = get_matching( self.Prev, p_add_p )
            at_t = get_matching( self, p_add_t )
            at_n = get_matching( self.Next, p_add_n )
            at_r3 = get_replacement( self.Next, p_rep_nn, p_con_nn, )
            if self.Next.Next:
                at_r4 = get_rep_2( self.Next, self.Next.Next, p_rep_nn_n, p_con_nn_n )
            if self.Prev.Prev:
                at_r5 = get_rep_2( self.Prev, self.Prev.Prev, p_rep_pp_p, p_con_pp_p )

            for at in at_p + at_t + at_n + at_r1 + at_r2 + at_r3 + at_r4 + at_r5:
                tmp_residue.add_atom( at )

        if residue or r:
            self.ready = tmp_residue
        if concap or c:
            self.con = tmp_residue
        if bridge or b:
            self.bri = tmp_residue

    def dist_to_res(self, other):
        return np.linalg.norm(other.get_com() - self.get_com())

    def get_com(self):
        com = np.zeros( 3 )
        m = 0
        for i in self:
            m += mass_dict[ i.element ]
            com += np.array([ float(i.x) , float(i.y) , float(i.z) ]) * mass_dict[ i .element ]
        return com / m

    def net_charge(self):
        tmp_props = Props()
        for at in self:
            tmp_props += at.Props
        return tmp_props["charge"]

    def hasNextNext(self):
        if self.Next:
            if self.Next.Next:
                return True
        return False
    def hasPrevPrev(self):
        if self.Prev:
            if self.Prev.Prev:
                return True
        return False

#Residue
    def __str__(self):
        base = "-".join( [self.Chain.chain_id, self.res_name + str(self.res_id)] )
        if self.concap:
            base += '-con'
        return base

    def get_relevant_residues(self):
        if self.n_term:
            return [self.ready.copy(), self.Next.ready.copy()]
        elif self.c_term:
            return [self.ready.copy(), self.Prev.ready.copy()]
        else:
            return [self.ready.copy(), self.Next.ready.copy(), self.Prev.ready.copy()]

    def get_relevant_concaps(self):
        if self.n_term:
            return [self.con.copy()]
        elif self.c_term:
            return [self.Prev.con.copy()]
        else:
            return [self.con.copy(), self.Prev.con.copy()]
    def mfcc_props(self):
        """After this functions all atoms here have final properties.
        """
        for at in self:
            p = molecules.Property()
            print "\nStarting with %s" %at.label
            print "Property at %.2f" % p['beta'][9]
            for R in self.get_relevant_residues():
                R.populate_bonds()
                for hx in R.get_dummy_h():
                    assert len ( hx.bonds ) == 1
                    hx.bonds.values()[0].Property += hx.Property
                    R.remove( hx )
                for R_at in R:
                    if R_at.label == at.label:
                        p += R_at.Property
                        print "Adding charge from res: %s, label: %s" %(R.res_name+str(R.res_id), R_at.label)
                        print p.q
            for C in self.get_relevant_concaps():
                C.populate_bonds()
                for hx in C.get_dummy_h():
                    assert len ( hx.bonds ) == 1
                    hx.bonds.values()[0].Property += hx.Property
                    C.remove( hx )
                for C_at in C:
                    if C_at.label == at.label:
                        p -= C_at.Property
                        print "Subtracting charge from %s" %C_at.label
                        print p.q
            print "Finished with %s" %at.label
            print p.q
            at.Property = p

    def add_atom( self, atom):
        self.append( atom )
        self.label_dict[ atom.label ] = atom
        atom.Molecule = self

    def get_dummy_h(self):
        """Return a list of dummy hydrogens that replaced heavy ones"""
        ats = []
        for i in self:
            if i.label.split('-')[-1] == 'XH':
                ats.append(i)
        if len(ats) == 0:
            print "Warning, no dummy hydrogens found in %s" % (self)
            return
        return ats
    def get_atom_by_label(self, label):
        try:
            for i in self:
                return self.label_dict[ label ]
        except KeyError:
            pass
    def existsPdbAtom(self, pdb_name ):
        for atom in self:
            if atom.pdb_name == pdb_name:
                return True
        return False

    def getAtom(self, pdb_name):
        #try:
        for i in self:
            if i.pdb_name == pdb_name:
                return i.copy()

    def copy_self(self):
        return self.copy()

    def copy( self):
        """Copy residue method"""
        new = Residue()

        [ new.add_atom( i.copy() ) for i in self ]
        new.c_term = self.c_term
        new.n_term = self.n_term
        new._res_id = self._res_id
        new._res_name = self.res_name
        new.AA = self.AA
#Keep information if this is concap
        new.concap = self.concap
        new._level = self._level


        new.Chain = self.Chain

        new.in_qm_region = self.in_qm_region
        new.in_mm_region = self.in_mm_region

        new.Next = self.Next
        new.Prev = self.Prev
        new.Bridge = self.Bridge

        return  new

    def copy_info(self):
        new = Residue()
        new.c_term = self.c_term
        new.n_term = self.n_term
        new.AA = self.AA

        new.in_qm_region = self.in_qm_region
        new.in_mm_region = self.in_mm_region

        new._res_id = self._res_id
        new._res_name = self.res_name
        new.Next = self.Next
        new.Prev = self.Prev
        new.Bridge = self.Bridge
        new.Chain = self.Chain
        return  new


    def add_props_to_atoms(self, prop_list):
#prop_list contains a list of atoms that have a Prop from the loprop MolFrag output
        for atm in self:
            for newatm in prop_list:
                if atm.label == newatm.label:
                    atm.Props = newatm.Props.copy()

    def get_xyz(self):
        """ Return string representation of all atoms coordinates for writing file"""
        string = "%d\n\n" %len(self)
        for i in self:
            string += "%s %s %s %s\n" %( i.element, i.x, i.y, i.z )
        return string.rstrip('\n')

    def get_mol_string(self, *args, **kwargs):
        return self.get_mol( *args, **kwargs)

    def get_mol(self, basis = ['ano-1 2','ano-1 4 3 1', 'ano-2 5 4 1' ] ):
        charge = 0

        elem = [x.element for x in self]
        atomtypes = len( uniq( elem ) )

        if "H" in elem:
            h_atoms = [ x for x in self if x.element == "H" ]
        if "C" in elem:
            c_atoms = [ x for x in self if x.element == "C" ]
        if "N" in elem:
            n_atoms = [ x for x in self if x.element == "N" ]
        if "O" in elem:
            o_atoms = [ x for x in self if x.element == "O" ]
        if "S" in elem:
            s_atoms = [ x for x in self if x.element == "S" ]

        if not self.is_concap:
            if self.n_term:
                charge += 1
            elif self.c_term:
                charge -= 1
            if res_dict[ self.res_name ] in chargeDict:
                charge += chargeDict[ res_dict[ self.res_name] ]
#Temporary fix to give level 3 concaps negative and positive even if they are not
#n or c terminal
        if self._level == 3:
            charge = 0
            for at in self:
                if at.pdb_name == "H1":
                    charge += 1
                    break
                if at.pdb_name == "OC1":
                    charge -= 1
                    break
                

        string = ""
        string += 'ATOMBASIS\n'
        string += 'The studied system is %s\n' % ( self ) 
        string += 'File generated by pdbreader 4.0 -- //By Ignat Harczuk --\n'
        string += 'Atomtypes=%s Charge=%s Angstrom Nosymm\n'%( str(atomtypes ) , str(charge) )
        
        if "H" in elem:
            string += 'Charge=1.0 Atoms=%d Basis=%s\n'%( len(h_atoms), basis[0])
            for k in h_atoms:
                string += '{0:15s}{1:10.5f}{2:10.5f}{3:10.5f}\n'.format( *tuple([ k.label, float(k.x), float(k.y), float(k.z) ]) )
        if "C" in elem:
            string += 'Charge=6.0 Atoms=%d Basis=%s\n'%( len(c_atoms), basis[1])
            for k in c_atoms:
                string += '{0:15s}{1:10.5f}{2:10.5f}{3:10.5f}\n'.format( *tuple([ k.label, float(k.x), float(k.y), float(k.z) ]) )
        if "N" in elem:
            string += 'Charge=7.0 Atoms=%d Basis=%s\n'%( len(n_atoms), basis[1])
            for k in n_atoms:
                string += '{0:15s}{1:10.5f}{2:10.5f}{3:10.5f}\n'.format( *tuple([ k.label, float(k.x), float(k.y), float(k.z) ]) )
        if "O" in elem:
            string += 'Charge=8.0 Atoms=%d Basis=%s\n'%( len(o_atoms), basis[1])
            for k in o_atoms:
                string += '{0:15s}{1:10.5f}{2:10.5f}{3:10.5f}\n'.format( *tuple([ k.label, float(k.x), float(k.y), float(k.z) ]) )
        if "S" in elem:
            string += 'Charge=16.0 Atoms=%d Basis=%s\n'%( len(s_atoms), basis[2])
            for k in s_atoms:
                string += '{0:15s}{1:10.5f}{2:10.5f}{3:10.5f}\n'.format( *tuple([ k.label, float(k.x), float(k.y), float(k.z) ]) )
        return string.rstrip( '\n' )

class NewChain( molecules.Cluster):
    """Will behaive like Cluster and rely everything on getters and setters to avoid bugs in overwriting properties"""
    def __init__(self, *args, **kwargs):
        self._snapshot = None
        self._time = None
        self._freq = None
        self._System = None
        super( NewChain, self ).__init__( *args, **kwargs )

    @property
    def System(self):
        if self._System:
            return self._System
        return None
    @System.setter
    def System(self, val):
        self._System = val

    @property
    def time(self):
        if self.System:
            return self.System.time
        return self._time

    @time.setter
    def time(self, val):
        self._time = val

    @property
    def snapshot(self):
        if self._System is not None:
            return self._System.snapshot
        return None

    @snapshot.setter
    def snapshot(self, val):
        self._snapshot = val

    @property
    def freq(self):
        if self._System is not None:
            return self._System.freq
        return None

    @freq.setter
    def freq(self, val):
        self._freq = val

    @property
    def chain_id(self):
        if self._chain_id is not None:
            return self._chain_id
        return None

    @chain_id.setter
    def chain_id(self, val):
        self._chain_id = val

        

class Chain( molecules.Cluster ):
    def __init__(self ):
        super( Chain, self).__init__()
        self.chain_id = 'A'
        self.res_count = 0

        self.ready_residues = []
        self.ready_concaps = []
        self.ready_bridges = []

        self.customDict = {'CRO':'X1', 'CRO2':'X2'}

        self.rootname = None

    def mpi_mfcc(self, procs = 4):
        """Undergoes full mfcc for all residues in this chain using MPI"""
        pass

    def copy( self):
        """Copy Chain method"""
        new = Chain()
        for mol in self:
            new.add_residue( mol.copy() )
        return new

    def center( self, atom = None ):
        """Translates the chain so that atom is in origo"""
        if atom is None:
            atom = self[0]
        vec = np.zeros(3) - atom.r
        for at in [at for mol in self for at in mol]:
            at.x, at.y, at.z = vec + at.r

    def plot_maya( self, p1, p2, p3, copy = True ):
        """Wrapper function to plot quivers of beta around whole chain
        p1 is the first point to put in origo to center around
        p2 is the point which will lie in the z axis 
        p3 will be projected to the x-axis
        """
        
        if copy:
            copy = self.copy()
        else:
            copy = self

        for mol in [mol for mol in copy if isinstance( mol, molecules.Molecule) ]:
            mol.populate_bonds()

        v, t1, t2, t3 = utilz.center_and_xz( p1, p2, p3 )

#Beta section
        b = utilz.ut2s( copy.p.b )

# plotting section
        for r in copy:
            r.translate_by_r( v )
            r.inv_rotate( t1, t2, t3 )

        x, y, z = utilz.E_at_sphere(r_max = 50, r_points = 10)
        bv = utilz.b_at_sphere( b, x, y, z )

        mlab.figure(figure=None, bgcolor=(1,1,1), fgcolor=None, engine=None, size=(400, 350))
        m = mlab.quiver3d( x, y, z, bv[...,0], bv[...,1], bv[...,2],
                colormap = 'BrBG' )

#Plot bonds
        for mol in copy:
            for each in mol.bond_dict:
                for key in mol.bond_dict[ each ]:
                    mlab.plot3d( [key.x, each.x],
                             [key.y, each.y],
                             [key.z, each.z], color = (0,0,0,) )
        for mol in copy:
            for at in mol:
                mlab.points3d([at.x], [at.y], [at.z], 
                        color = color_dict[ at.element ],
                        resolution = 50,
                        scale_factor = scale_factor_dict[ at.element] )

    def __str__(self):
        """docstring for __str__"""
        return "chain_" + self.chain_id
        
    def save(self, name = False):
        if name:
            fobj = open( "%s" %name, 'wb' )
        else:
            fobj = open( "%s.p" %self, 'wb' )
        pickle.dump( self, fobj, protocol = 2 )
        fobj.close()

    def write_pro( self, FILE ):
        f =  open( FILE , 'w' )
        for i in self:
            if i.in_mm_region:
                for atm in i:
                    if atm.in_qm_region:
                        continue
                    f.write( "{0:15s}{1:10s}{2:10s}{3:10s}\n".format(atm.label, atm.x, atm.y, atm.z ) )

    def write_mol_concaps(self, args ):
        for res in self.ready_concaps:
            name = args.pdbfile.rstrip(".pdb") 
            _file = open( "-".join( [ name, "concap", res.res_name, res.res_id ] ) + '.mol' ,'w')
            _file.write( res.get_mol() )

    def write_mol_bridges(self, args ):
        for res in self.ready_bridges:
            name = args.pdbfile.rstrip(".pdb") 
            _file = open( "-".join( [ name, "bridge", res.res_name, res.res_id ] ) + '.mol' ,'w')
            _file.write( res.get_mol() )

    def write_mol_residues( self, args ):
        for res in self.ready_residues:
            name = args.pdbfile.rstrip(".pdb") 
            _file = open( "-".join( [ name, "residue", res.res_name, res.res_id ] ) + '.mol' ,'w')
            _file.write( res.get_mol( ) )

    def add_residue ( self, residue ):
        if isinstance( residue, Residue ):
            self.append( residue )
            residue.Chain = self
            self.res_count += 1

    def connect_residues(self,):

        for i in range( len(self) ):

            if self[i].n_term:
                self[i].Next = self[i + 1]

            elif self[i].c_term:
                self[i].Prev = self[i - 1]

            else:
                self[i].Next = self[i + 1]
                self[i].Prev = self[i - 1]

    def build_bridges(self, level = 1):
        Pat = Pattern( level )
        for j in self:
            if j.Bridge:
                tmp_residue = j.copy() 
                tmp_residue = Residue()
                tmp_residue.AA = Residue.AA
                tmp_residue.concap = True
                tmp_residue.bridge = True
                for k in j:
                    if Pat.bri_b.match( k.pdb_name ):
                        tmp_residue.add_atom( k.copy() )
                    elif Pat.rep_bri_b_b.match( k.pdb_name ):
                        tmp_atom = k.copy()
                        tmp_atom.element = "H"
                        tmp_atom._label += "-XH"
                        H = tmp_atom.r
                        for others in (Pat.con_bri_b_b[ tmp_atom.pdb_name ] ):
                            C = j.getAtom( others )
                            dist = Pat.pdb_dist_dict[ C.pdb_name ]
                            C = C.r
                            tmp_atom.setXyz (C + (H - C) * dist / np.linalg.norm(H - C))
                            tmp_residue.add_atom( tmp_atom )
                for k in j.Bridge:
                    if Pat.bri_b.match( k.pdb_name ):
                        tmp_residue.add_atom( k.copy() )
                    elif Pat.rep_bri_b_b.match( k.pdb_name):
                        tmp_atom = k.copy()
                        tmp_atom.element = "H"
                        tmp_atom._label += "-XH"
                        H = tmp_atom.r
                        for others in (Pat.con_bri_b_b[ tmp_atom.pdb_name ] ):
                            C = j.Bridge.getAtom( others )
                            dist = Pat.pdb_dist_dict[ C.pdb_name ]
                            C = C.r
                            tmp_atom.setXyz (C + (H - C) * dist / np.linalg.norm(H - C))
                            tmp_residue.add_atom( tmp_atom )
                self.ready_bridges.append( tmp_residue )
class System( list ):
    def __init__(self):
        pdbfile = None
        pass

    @staticmethod
    def read_protein_from_file( FILE, in_AA = True ):
        pat = re.compile(r'^ATOM|^HETATM|^TER|^END')
        text = [f for f in open(FILE).readlines() if pat.match(f)]
        """ 
        return all Protein chains of the PDB file in a cla System
        """
        pat_element = re.compile( r'([A-Z])' )
#tmp is tracking array used to keep track of res_name changing
        tmp = []
        firstEntry = True

#Initiate a chain with chain_id corresponding to "A" usually in the PDB file
        tmp_chain = Chain()
        tmp_chain.chain_id = text[0][21:22].strip()

#Add the chain to System
        world = System()
        world.pdbfile = FILE
        world.basefile = re.compile(r'(.*).(pdb|gro)').search( FILE ).group(1)
        world.profile = world.basefile + ".pro"
        world.potfile = world.basefile + ".pot"

        world.append( tmp_chain )

        tmp_residue = Residue( AA = in_AA )

        for i in range(len( text )):
#Have to define this here so that first and last amino acid has correct name
            res_name = text[i][17:21].strip()
            try:
                res_id = int( text[i][22:26].strip() )
            except ValueError:
                pass
#Initiate values for Atom,
            chain_id = text[i][21:22].strip()
            x = text[i][30:38].strip()
            y = text[i][38:46].strip()
            z = text[i][46:54].strip()
            pdb_name = text[i][11:16].strip()

            try:
                element = pat_element.search( pdb_name ).group(1)
#This will cause IndexError for TER when no other entries are there
            except IndexError:
                element = None
            except AttributeError:
                element = None

#This is end of chain, need to put residue in chain and start at beginning for new chain
            if text[i][0:3] == "TER" or text[i][0:3] == "END":

                tmp_residue.c_term = True
                tmp_chain.add_residue ( tmp_residue )
                tmp_residue.Chain = tmp_chain
                tmp_residue = Residue(AA = in_AA)
                tmp = []
                firstEntry = True
                continue

            if ( chain_id != tmp_chain.chain_id):

                tmp_chain = Chain()
                tmp_chain.chain_id = chain_id
                world.append( tmp_chain )

                tmp_residue = Residue(AA = in_AA)
                tmp = []
                firstEntry = True


#Check to only make acids out of predefiend 4-letter keywords in global known_res dictionary
            if text[i][17:21].strip() not in res_dict:
                continue
# First residue added to tracking array
            if firstEntry:
                tmp_residue.n_term = True
                tmp.append( int( text[i][22:26].strip() ) )
                firstEntry = False
                
# If tracking arrays last element is not the same as current, means we have new residue
            if tmp[-1] != res_id :
                tmp_chain.add_residue( tmp_residue )
                tmp_residue.Chain = tmp_chain
                tmp_residue = Residue( AA = in_AA)
                tmp.append( res_id )


#Create an Atom here and add it to the residue
            tmp_atom = Atom()
            tmp_atom.x, tmp_atom.y, tmp_atom.z = map( float, [x, y, z] )
            tmp_atom.element = element
            tmp_atom.pdb_name = pdb_name
            tmp_atom._res_name = res_name
            tmp_atom._res_id = res_id
            tmp_atom._label = "%d-%s-%s" %( res_id, res_name, pdb_name )

            tmp_residue.add_atom( tmp_atom )
            tmp_residue._res_id = res_id
            tmp_residue._res_name = res_name
        return world

    def set_qm_and_mm_regions( self, qm_res_ids = [], cutoff = 5.0 ):
        for ch in self:
            for res in ch:
                res.in_mm_region = False

        for ch in self:
            for res in ch:
                if int(res.res_id) in qm_res_ids:
                    res.in_qm_region = True

        for ch in self:
            for res in ch:
                if res.in_qm_region:
                    for res2 in [r for r in ch if r != res ]:
                        if res.dist_to_res(res2) < cutoff:
                            res2.in_mm_region = True

    def update_atoms_in_qm_region( self ):
        """ This function takes the whole systems MM atoms, and if they are included in the QM
        set their label to in_qm_region, so that the write_pro command excludes them"""

#set all atoms in qm to in_qm_region
        qm_list = []
        for ch in self:
            for res in ch:
                if res.in_qm_region:
                    for atm in res.ready:
                        qm_list.append( atm.label )

#set atoms in MM region to QM
        for ch in self:
            for res in ch:
                if res.in_mm_region:
                    for atm in res:
                        if atm.label in qm_list:
                            atm.in_qm_region = True


    def get_targz_and_mol(self):
        """Routine to find mol and .tar.gz files from DALTON
        that lead to the target .POT file"""
        tmp = []
        tmpMols = self.get_mol_files()
        for i in self.get_targz_files():
            for j in tmpMols:
                if i.split( "_" )[1].rstrip(".tar.gz") == j.rstrip (".mol"):
                    tmp.append( j.rstrip ('.mol'))
        return tmp

    def get_net_charge(self, args):
        chainDict = {}
        for ch in self:
            chainDict[ ch.chain_id ] = 0.0
            for res in ch:
                if res_dict[ res.res_name ] in chargeDict:
                    chainDict[ ch.chain_id ] += chargeDict[ res_dict[ res.res_name] ]
        return chainDict

    def write_xyz(self, ind ):
        for ch in self:
            for co in ch.ready_bridges:
                if co.res_id == ind:
                    open( "-".join( [ "bridge", co.res_name, co.res_id ] ) + ".xyz", 'w').write( co.get_xyz() )
                    print "Wrote bridge-%s-%s"%(co.res_name, co.res_id) + ".xyz" 

            for co in ch.ready_concaps:
                if co.res_id == ind:
                    open( "-".join( [ "concap", co.res_name, co.res_id ] ) + ".xyz", 'w').write( co.get_xyz() )
                    print "Wrote concap-%s-%s"%(co.res_name, co.res_id) + ".xyz" 

            for co in ch.ready_residues:
                if co.res_id == ind:
                    open( "residue-%s-%s"%(co.res_name, co.res_id) + ".xyz", 'w').write( co.get_xyz() )
                    print "Wrote residue-%s-%s"%(co.res_name, co.res_id) + ".xyz" 

    def get_mol_files(self):
        molFiles = [ f for f in os.listdir( os.getcwd() ) if f.endswith( ".mol") ]
        return molFiles

    def get_targz_files(self):
        targzFiles = [ f for f in os.listdir( os.getcwd() ) if f.endswith( ".tar.gz") ]
        return targzFiles

    def find_sulfur_bridges(self):
        cys = []
        for ch in self:
            for res in ch:
                if res.res_name == "CYS":
                    if res.existsPdbAtom( "HG1" ):
                        #No bridge if has this atom
                        continue
                    else:
                        cys.append( res )

#For all cys that form sulfur bridge, find the partner residue 
        for i in range(len(cys)):
            for j in range( i , len(cys)):
                if i == j:
                    continue
                if cys[i].getAtom( "SG" ).dist_to_atom( cys[j].getAtom( "SG" ) ) < 2.8 :
                    cys[i].Bridge = cys[j]
                    cys[j].Bridge = cys[i]

    def paralell_molfrag( self, name, global_array, dal_prefix = "linear_" ):
        """ 
        Special ruotine that is called to parallellize write_pot.

        All name files shared by a pool of workers is obtained by matching all .mol
        and .tar.gz files in current directory.

        Name is name of molfile that has a corresponding .tar.gz file with dal_prefix in
        front of it.

        """
        mol = name + ".mol"
        targz = dal_prefix + name + ".tar.gz"

    def write_pot(self, FILE, qmmm_border = 1.0, level = 1, ncpu = 1, scratch = '/tmp'):
        """ 
        
        1) Opens all matching .tar.gz and .mol files in directory,
        and looping over all atoms in mm_region and sets their properties.

        2) For mm_region_atoms within qmmm_co of qm region,
        transfers properties to next neightbour.

        NOTE: Try to paralellize since MolFrag.output_potential_file takes most time.
        """

        f_warnings = open( "warnings.log", 'w' )

# Separate MM and QM atoms for easier property transfer
        mm_atoms = [ atm for ch in self for res in ch for atm in res
            if res.in_mm_region and not atm.in_qm_region] 

        qm_atoms = [ atm for ch in self for res in ch for atm in res.ready
            if res.in_qm_region] 

        for i in qm_atoms:
            for j in mm_atoms:
                r =  i.dist_to_atom( j )
                if r < qmmm_border:
                    j.in_qmmm_border = True

#First implement level 2 capping, define all connections for level 2 capping
# This is to know which atom to transfer the properties from the dummmy XH
# before summing over all common atoms begins
#
# Since the pdblabel of XH atom that replaced a CA will be CA, it means that
# at level 1, we will try to move properties from the XH to either C or N.
# 
# (at level 3, if multiple options are presented, move properties equally)
#

        ints = map( float, range(-1,2) )
        if level == 1:
            connect_dict = { "CA" : [ "C" , "N" ] , "CD" : ["N"] , "CA1" : [ "N1" ] , 
                    "CA3" : [ "C3" ], "CB" : ["SG"] }

        if level == 2:
            connect_dict = { "C1" : ["CA1"] , "CD": ["N"] , "CB" : ["CA"],
                    "N"   : ["CA"] , "N3" : ["CA3"] , "C" : ["CA"] ,
                    "CG"  : ["CA"] , "CA" : ["C", "N"], 
                    "CB1" : ["CA"] }

        pat_xyz = re.compile(r'^\s*(\w|-)+\s+(-*\d*.+\d+)\s+(-*\d*.+\d+)\s+(-*\d*.+\d+) *$')

# Paralellize over all residues in get_targz_and_mol later
#
        for each in self.get_targz_and_mol():
            # read all the atoms in the .mol file

            target_dir =  os.path.join( scratch , "loprop", each )
            tmp_res = Residue()

            for line in open( each + ".mol").readlines():
                if pat_xyz.match( line ):
                    f = line.split()
                    tmp_atom = Atom()
                    tmp_atom.AA = True
                    tmp_atom.x = float(f[1])
                    tmp_atom.y = float(f[2])
                    tmp_atom.z = float(f[3].strip())
                    tmp_atom.element = f[0].split('-')[-1][0]
                    tmp_atom.pdb_name = f[0].split('-')[-1]
                    tmp_atom._label = f[0]
                    res_name = each.split('-')[3]
                    res_id = each.split('-')[4]
                    tmp_res.add_atom ( tmp_atom )
                    tmp_res._res_name = res_name
                    tmp_res.res_id = res_id

            if not os.path.isdir( target_dir ):
                os.makedirs( target_dir )

            tarfile.open( "linear_" + each + ".tar.gz", 'r:gz').extractall(path=target_dir)
            pot_lines = MolFrag( tmpdir = target_dir,
                    pf = penalty_function( alpha = 2.0 ),
                    sf = shift_function).output_potential_file(
                            maxl = 0,
                            pol = 1
                            ).split("\n")
            pat_pot = re.compile( r'^1 ' )
            relevant = [line for line in pot_lines if pat_pot.match( line )]

            assert len(relevant) == len(tmp_res)

            for i in range(len(relevant)):
                tmp_res[ i ].Props["charge"] =  float(relevant[i].split()[4]) 
                tmp_res[ i ].Props["polar"] =   float(relevant[i].split()[5]) 

# Small code snippet to check so that the read residue/concap has integer charge
            o = 0
            for oo in [-1.0, 0.0, 1.0]:
                if almost_eq( oo, tmp_res.net_charge(), thr=0.01 ):
                    o += 1
            if o != 1:
                print "ERROR item %s does not have integer charge" %each
                print tmp_res.net_charge()
                raise SystemExit

            for ch in self:
                if "residue" in each:
                    for orig in [res.ready for res in ch if res.in_mm_region]:
                        if orig.res_id == each.split('-')[-1]:
                            orig.add_props_to_atoms( tmp_res )
                elif "concap" in each:
                    for orig in [res.con for res in ch if res.in_mm_region]:
                        if orig.res_id == each.split('-')[-1]:
                            orig.add_props_to_atoms( tmp_res )

# First transfer the properties of hydrogens that originally replaced a heavy atom
            for ch in self:
                for res in ch:
                    for atm in res.ready:
                        # -XH is hydrogen that replaced heavy atom
                        if atm.label.endswith("-XH"):
# atm.pdb_name still has the original heavy atoms pdb entry
# entry will be those atoms which are to have the hydrogens properties
# loop over all entries twice
                            for entry in connect_dict[ atm.pdb_name ]:
                                try:
                                    alabel = "-".join( atm.label.split('-')[0:2]+[ entry ] )
                                    res.get_atom_by_label( alabel ).Props = Props()
                                    atm.transfer_props( res.get_atom_by_label( alabel ) )
                                except AttributeError:
                                    f_warnings.write( "Could not get the connected heavy from dict, attrib error: \n" )
                                    f_warnings.write( "alabel: %s\n" % alabel )
                                    f_warnings.write("atm.label: %s\n" % atm.label)
                                    f_warnings.write( "res.res_name: %s res.res_id: %s\n\n" %( res.res_name, res.res_id ) )
            for ch in self:
                for con in [res.con for res in ch if res.con is not None]:
                    for atm in con:
                        if atm.label.endswith("-XH"):
                            for entry in connect_dict[atm.pdb_name]:
                                try:
                                    alabel = "-".join( atm.label.split('-')[0:2]+[ entry ] )
                                    con.get_atom_by_label( alabel ).Props = Props()
                                    atm.transfer_props( con.get_atom_by_label( alabel ) )
                                except AttributeError:
                                    f_warnings.write("Could not get the connected heavy from dict, attrib error: \n" )
                                    f_warnings.write("alabel: %s\n" %alabel )
                                    f_warnings.write("atm.label: %s\n" % atm.label)
                                    f_warnings.write("res.res_name: %s res.res_id: %s\n\n" %( res.res_name, res.res_id ) )

# Go through each label read from atoms in mm_region,
# add its properties from each encountered residue and 
# subtract it from each encountered concap
        final = []

        for atom in mm_atoms:
            if atom.in_qmmm_border:
                if atom.element == "N":
                    if len( atom.get_closest( cutoff =2.0 ) ) != 2:
                        print "ERROR: qmmm border %s had more than 2 neighbors" %atom
                atom.transfer_own_props_to_list( atom.get_closest( cutoff = 2.0) )
            label = atom.label
            tmp_props = Props()
            for ch in self:
                for res in [r.ready for r in ch if r.in_mm_region]:
                    atm = res.get_atom_by_label ( label )
                    if atm is not None:
                        tmp_props += atm.Props

                for res in [res.con for res in ch if res.con is not None]:
                    atm = res.get_atom_by_label ( label )
                    if atm is not None:
                        tmp_props -= atm.Props

            final.append( " ".join( [atom.label] +map(str,[atom.x, atom.y, atom.z]) + map(str, [tmp_props["charge"], tmp_props["polar"] ] ) ))

        string = "".join( ["AA\n", "%d 0 1 1\n" %len(final)] +\
                ["{0:15s}{1:10s}{2:10s}{3:10s}{4:10s}{5:10s}\n".format(\
                line.split()[0].split('-')[0], line.split()[1], line.split()[2], line.split()[3], line.split()[4], line.split()[5]) 
                for line in final ])

        
        f = open( FILE, 'w' )
        f.write( string )
        return string


    def write_qmmm( self, co = 5, qmmm_xyz = False ):
        """ """
        name = self.pdbfile.rstrip(".pdb") 
        for ch in self:
            for res in ch:
#For all residues defined by qm region
                if res.in_qm_region:
                    _file = open( "-".join( [ name, "qm", "residue", res.res_name, res.res_id ] ) + '.mol' ,'w')
                    _file.write( res.ready.get_mol() )

                elif res.in_mm_region:
                    _file = open( "-".join( [ name, "mm", "residue", res.res_name, res.res_id ] ) + '.mol' ,'w')
                    _file.write( res.ready.get_mol() )

#residues i, i+1 and i-1 in mm region
                    _file = open( "-".join( [ name, "mm", "residue", res.Next.res_name, res.Next.res_id ] ) + '.mol' ,'w')
                    _file.write( res.Next.ready.get_mol() )

                    _file = open( "-".join( [ name, "mm", "residue", res.Prev.res_name, res.Prev.res_id ] ) + '.mol' ,'w')
                    _file.write( res.Prev.ready.get_mol() )

#concaps i and i-1 in mm_region
            for res in [r.con for ch in self for r in ch if r.con is not None]:
                if res.in_mm_region:
                    _file = open( "-".join( [ name, "mm", "concap", res.res_name, res.res_id ] ) + '.mol' ,'w')
                    _file.write( res.get_mol() )
                    _file = open( "-".join( [ name, "mm", "concap", res.Prev.res_name, res.Prev.res_id ] ) + '.mol' ,'w')
                    _file.write( res.Prev.con.get_mol() )

        if qmmm_xyz:
            resid = None
            for ch in self:
                atoms = []
                for res in ch:
                    if res.in_qm_region:
                        resid = res.res_id
                    if res.in_qm_region:
                        for at in res.ready:
                            atoms.append(  at.xyz() )
                    elif res.in_mm_region:
                        for at in res.ready:
                            atoms.append(  at.xyz() )
            xyz = self.basefile + "_%sqm_%.2fco.xyz" %( resid, co )
            f_ = open( xyz, 'w' )
            f_.write("%d\n\n" %len(atoms))
            [f_.write( at + "\n" ) for at in atoms ]

class NewSystem( list ):
    """Can hold instances of Clusters, Molecules, Atoms in it
    used to separate trajectory configurations between different snapshots.

    Each trajectory snapshot thus holds a NewSystem class
    
    """
    def __init__(self, *args, **kwargs):
        self._snapshot = None
        self._time = None
        self._freq = None
        super(NewSystem, self).__init__()

        if type(args) == tuple:
            if len(args) == 1:
                if type(args[0]) == list:
                    for i in args[0]:
                        self.add( i )
                else:
                    self.add( args[0] )
            else:
                for at in args:
                    self.add( at )

    def connect_everything(self):
        for atom in self.atoms:
            atom.System = self
        for mol in self.molecules:
            mol.System = self
            mol.connect_everything()
        for cluster in [c for c in self if isinstance( c, molecules.Cluster) ]:
            cluster.System = self
            cluster.connect_everything()

# To define N- /C- terminals and set Next/ Prev attributes
    def connect_residues(self,):
        for ch in self:
            ch[0].n_term = True
            ch[-1].c_term = True
            for i, res in enumerate(ch):
                if ch[i].n_term:
                    ch[i]._Next = ch[i + 1]

                elif ch[i].c_term:
                    ch[i]._Prev = ch[i - 1]

                else:
                    ch[i]._Next = ch[i + 1]
                    ch[i]._Prev = ch[i - 1]
    @property
    def A(self):
        return self.get_chain_by_char ('A')
    @property
    def B(self):
        return self.get_chain_by_char ('B')
    @property
    def C(self):
        return self.get_chain_by_char ('C')
    @property
    def X(self):
        return self.get_chain_by_char ('X')
    
    def get_chain_by_char(self, char):
        for ch in self:
            if ch.chain_id == char:
                return ch
        print "No chain with identifier %s in %s" %(char, self )
        return

    @staticmethod
    def get_full_protein_from_string( _string, in_AA = True ):
        pat = re.compile(r'^ATOM|^HETATM|^TER|^END')
        text = _string.split('\n')
        """ 
        return all Protein chains of the PDB file in a cla System
        """
        pat_element = re.compile( r'([A-Z])' )
#tmp is tracking array used to keep track of res_name changing
        tmp = []
        firstEntry = True

#Initiate a chain with chain_id corresponding to "A" usually in the PDB file
        tmp_chain = Chain()
        tmp_chain.chain_id = text[0][21:22].strip()

#Add the chain to System
        world = NewSystem()

        world.append( tmp_chain )

        tmp_residue = Residue( AA = in_AA )

        for i in range(len( text )):
#Have to define this here so that first and last amino acid has correct name
            res_name = text[i][17:21].strip()
            try:
                res_id = int( text[i][22:26].strip() )
            except ValueError:
                pass
#Initiate values for Atom,
            chain_id = text[i][21:22].strip()
            x = text[i][30:38].strip()
            y = text[i][38:46].strip()
            z = text[i][46:54].strip()
            pdb_name = text[i][11:16].strip()

            try:
                element = pat_element.search( pdb_name ).group(1)
#This will cause IndexError for TER when no other entries are there
            except IndexError:
                element = None
            except AttributeError:
                element = None

#This is end of chain, need to put residue in chain and start at beginning for new chain
            if text[i][0:3] == "TER" or text[i][0:3] == "END":

                tmp_residue.c_term = True
                tmp_chain.add_residue ( tmp_residue )
                tmp_residue.Chain = tmp_chain
                tmp_residue = Residue(AA = in_AA)
                tmp = []
                firstEntry = True
                continue

            if ( chain_id != tmp_chain.chain_id):

                tmp_chain = Chain()
                tmp_chain.chain_id = chain_id
                world.append( tmp_chain )

                tmp_residue = Residue(AA = in_AA)
                tmp = []
                firstEntry = True


#Check to only make acids out of predefiend 4-letter keywords in global known_res dictionary
            if text[i][17:21].strip() not in res_dict:
                continue
# First residue added to tracking array
            if firstEntry:
                tmp_residue.n_term = True
                tmp.append( int( text[i][22:26].strip() ) )
                firstEntry = False
                
# If tracking arrays last element is not the same as current, means we have new residue
            if tmp[-1] != res_id :
                tmp_chain.add_residue( tmp_residue )
                tmp_residue.Chain = tmp_chain
                tmp_residue = Residue( AA = in_AA)
                tmp.append( res_id )


#Create an Atom here and add it to the residue
            tmp_atom = Atom()
            tmp_atom.x, tmp_atom.y, tmp_atom.z = map( float, [x, y, z] )
            tmp_atom.element = element
            tmp_atom.pdb_name = pdb_name
            tmp_atom._res_name = res_name
            tmp_atom._res_id = res_id
            tmp_atom._label = "%d-%s-%s" %( res_id, res_name, pdb_name )

            tmp_residue.add_atom( tmp_atom )
            tmp_residue._res_id = res_id
            tmp_residue._res_name = res_name
        world.connect_everything()
        return world




    def find_sulfur_bridges(self):
        """Note: identical to System.find_sulfur_bridges"""
        cys = []
        for ch in self:
            for res in ch:
                if res.res_name == "CYS":
                    if "HG1" in res:
                        #No bridge if has this atom
                        continue
                    else:
                        cys.append( res )

#For all cys that form sulfur bridge, find the partner residue 
        for i in range(len(cys)):
            for j in range( i , len(cys)):
                if i == j:
                    continue
                if cys[i].get_atom( "SG" ).dist_to_atom( cys[j].get_atom( "SG" ) ) < 2.8 :
                    cys[i]._Bridge = cys[j]
                    cys[j]._Bridge = cys[i]

    @property
    def snapshot(self):
        if self._snapshot:
            return self._snapshot
        return None
    @snapshot.setter
    def snapshot(self, val):
        self._snapshot = val

    @property
    def time(self):
        if self._time:
            return self._time
        return None
    @time.setter
    def time(self, val):
        self._time = val

    @property
    def freq(self):
        if self._freq:
            return self._freq
        return None
    @freq.setter
    def freq(self, val):
        self._freq = val


#NewSystem method of adding objects
    def add(self, item):
        if isinstance( item, molecules.Cluster):
            self.append( item )
            item.System = self
        elif isinstance( item, molecules.Molecule):
            c = NewChain( item )
            c.System = self
            self.append( c )

    @property
    def atoms(self):
        return [a for chain in self for mol in chain for a in mol if isinstance(a , molecules.Atom ) ]
    @property
    def molecules(self):
        return [m for chain in self for m in chain if isinstance(m , molecules.Molecule ) ]
    @classmethod
    def from_pdb_string( cls, _string ):
        """Assuming the string is a pdb format, read in all chains and stuff"""
        res, meta = all_residues_from_pdb_string( _string )

        chains = [ NewChain(ch) for ch in utilz.splitter( res, lambda x: x.chain_id )]

        system = cls( *chains )
        system.meta = meta

        for ch in chains:
            ch[0].n_term = True
            ch[-1].c_term = True
        for ch in chains:
            ch.connect_everything()
            ch.System = system
        for chain in chains:
            for res in chain:
                chain.chain_id = res.chain_id
        return system

    def save(self, fname = "system.p"):
        pickle.dump( self, open(fname, 'wb' ), protocol = 2 )
    @staticmethod
    def load(fname = 'system.p'):
        if not os.path.isfile( fname ):
            raise IOError
        pick = pickle.load( open(fname, 'rb' ) )
        if not isinstance( pick, NewSystem ):
            raise TypeError("Wrong pickled class")
        return pick
    
        
class World( list ):
    """Can hold instances of Systems, Clusters, Molecules, Atoms in it"""
    def __init__(self, *args, **kwargs):
        super(World, self).__init__()
        self._project = None

        if args is not None:
            for each in args:
                if isinstance( each, NewSystem):
                    self.append( each )
        
    def connect(self):
        for s in [s for s in self if isinstance(s, NewSystem) ]:
            s.connect_everything

    def add(self, item):
        if isinstance(item, NewSystem ):
            self.append( item )
 
    def save(self, fname = "world.p"):
        pickle.dump( self, open(fname, 'wb' ), protocol = 2 )

    
    @property
    def molecules(self):
        return [m for s in self for ch in s for m in ch if isinstance( m, molecules.Molecule ) ]

    @staticmethod
    def load(fname = 'world.p'):
        if not os.path.isfile( fname ):
            raise IOError
        pick = pickle.load( open(fname, 'rb' ) )
        if not isinstance( pick, World ):
            raise TypeError("Wrong pickled class")
        return pick
    



def main():

    """
    By default, this program stores and entire pdb file as a target POT file
    used in DALTON pe /qmmm.

    By specifying -qm <res_id>, by default stores a .POT file for the
    residues found at a center-of-mass distance within 
    -co --cutoff <dist> of residue with number <res_id> 

    """
    args = run_argparse( sys.argv )

#System is a representation of the pdbfile supplied, i.e. separate split chains, waters, etc
    Syst = System.read_protein_from_file( args.pdbfile )
    Syst.find_sulfur_bridges()


#Build a list of molecules containing atoms to be 
#printed to .mol files for the residues, concaps and bridges

    if args.qm:
        Syst.set_qm_and_mm_regions( qm_res_ids = args.qm, cutoff = args.cutoff )

    [ res.connect_residues() for res in Syst ]

    [ res.build_residues( level = args.level ) for res in Syst ]
    [ res.build_concaps(  level = args.level ) for res in Syst ]
    [ res.build_bridges(  level = args.level ) for res in Syst ]

    Syst.update_atoms_in_qm_region()


#Write all files if requested
    if args.write:
        if "mol" in args.write:
            [ res.write_mol_residues( args, ) for res in Syst ]

        if "con" in args.write:
            [ ch.write_mol_concaps(args,) for ch in Syst ]

        if "bri" in args.write:
            [ res.write_mol_bridges(args,) for res in Syst ]

        if "pro" in args.write:
            [ch.write_pro( Syst.profile ) for ch in Syst ]

        if "pot" in args.write:
            Syst.write_pot( FILE = Syst.potfile,
                    qmmm_border = 1.0, level = args.level, 
                    ncpu = 4 )

        if "qmmm" in args.write:
            Syst.write_qmmm( co = args.cutoff, qmmm_xyz = args.qmmm_xyz )

    if args.write_xyz:
        for i in args.write_xyz:
            Syst.write_xyz( i )




if __name__ == "__main__":
    main()
