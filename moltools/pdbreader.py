#!/usr/bin/env python

__all__ = [ 'Pattern', 'Residue', 'Chain', 'System', 'World',
        'all_residues_from_pdb_string',
        'all_residues_from_pdb_file',
        ]

import os, re, sys, argparse, tarfile, ctypes, multiprocessing, pickle, logging

import numpy as np
from copy import deepcopy
#from mayavi import mlab

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

from .molecules import Molecule, Cluster, Bond, charge_dict, a0
from .property import Property
from .molecules import Atom as MoleculesAtom
from .utilz import unique, splitter

try:
    from .loprop.loprop import MolFrag, penalty_function, shift_function
except ImportError:
    pass

try:
    from applequist import gaussian
except ImportError:
    pass


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


pdb_dist_dict = { "CA" : 1.1, "CD" : 1.1, 
        "SG" : 1.1, "N" : 1.1, "C": 1.1,
        "N1" : 1.1, "C3" : 1.1, "CA1" : 1.1, "CA3" : 1.1,
        "CB" : 1.1, "C1": 1.1 , "N3" : 1.1,
        "CB1" :1.1, "CB2" :1.1,
}




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
    to the same chain type
    
    Will set all ._label properties for all atoms from here on
    """
    pat = re.compile(r'^ATOM|^HETATM')
    pat_meta = re.compile(r'^TITLE|^REMARK|^CRYST1|^MODEL')
    text = [f for f in _string.split('\n') if pat.match(f)]
    meta_text = [f for f in _string.split('\n') if pat_meta.match(f)]

    atoms = []
    res_ids = []
    chain_ids = []
    chain_dict = {}
    pat_element = re.compile( r'([A-Z]|[a-z])' )

    for i, line in enumerate( text ):
        res_name = text[i][17:21].strip()
        res_id = int( text[i][22:26].strip() )
        pdb_name = text[i][11:16].strip()
#If some digits used in specific pdb naming like gromacs/avogadro
        element = pat_element.search( pdb_name ).groups()[0]
        chain_id = text[i][21:22].strip()
        if chain_id == "":
            chain_id = "X"
        x = text[i][30:38].strip()
        y = text[i][38:46].strip()
        z = text[i][46:54].strip()
        x, y, z = map( float, [x, y, z] )

        atoms.append( Atom( x = x, y = y, z = z,
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

    res_ids = unique( res_ids )
    chain_ids = unique( chain_ids )

    res = [Residue([a for a in atoms if (a.res_id == r and a.chain_id == c)], AA = in_AA) for c in chain_ids for r in res_ids if r in chain_dict[c] ]

    for each in res:
        each.res_id = each[0].res_id
        each.res_name = each[0].res_name
        each._chain_id = each[0].chain_id
        for at in each:
#Pick .label by default since it is connected to this residue
            at._label = at.label
    

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
            dist = pdb_dist_dict[ C.pdb_name ]
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
            dist = pdb_dist_dict[ C.pdb_name ]
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
        with distance according to the heavy atom types distance in pdb_dist_dict

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
        self[ ( 'reg', 'res', 'n', 'add', 2 ) ] = re.compile(r'N$|H$|HN$|2HN$|CA$|HA$|HA1$|HA2$')
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
        self[ ( 'pro', 'con', 'pp_p', 'con', 2 ) ] = {}

        self[ ( 'reg', 'con', 'nn_n', 'con', 2 ) ] = { }
        self[ ( 'reg', 'con', 'p_t', 'con', 2 ) ] = { }
        self[ ( 'reg', 'con', 'pp_p', 'con', 2 ) ] = {}

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
        self[ ( 'reg', 'res', 'nn', 'rep', 2 ) ] = re.compile(r'C$|CB$')
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
        self[ ( 'reg', 'con', 'nn', 'rep', 2 ) ] = re.compile(r'C$|CB$')
        self[ ( 'cus', 'con', 'nn', 'rep', 2 ) ] = re.compile(r'C2$|CB1')
        self[ ( 'pro', 'con', 'nn', 'rep', 2 ) ] = re.compile(r'CB$|CD$|CD2$|C$')
        self[ ( 'reg', 'con', 'bb', 'rep', 2 ) ] = re.compile(r'CA$')

        #Extra dummy patterns
        self[ ( 'pro', 'con', 'nn_n', 'rep', 2 ) ] = re.compile(r'DUMMY')
        self[ ( 'pro', 'con', 'p_t', 'rep', 2 ) ] = re.compile(r'DUMMY')
        self[ ( 'pro', 'con', 'pp_p', 'rep', 2 ) ] = re.compile(r'DUMMY')

        self[ ( 'reg', 'con', 'nn_n', 'rep', 2 ) ] = re.compile(r'DUMMY')
        self[ ( 'reg', 'con', 'p_t', 'rep', 2 ) ] = re.compile(r'DUMMY')
        self[ ( 'reg', 'con', 'pp_p', 'rep', 2 ) ] = re.compile(r'DUMMY')

        #Extra dummy patterns



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

class Atom( MoleculesAtom ):
    def __init__(self, *args, **kwargs):
        self._chain_id = None
        self._name = None
        self._res_name = None
        super( Atom, self ).__init__( *args, **kwargs )

        if kwargs != {}:
            setattr( self, "_chain_id",  kwargs.get( "chain_id", None ) )
            setattr( self, "_res_name",  kwargs.get( "res_name", None ) )

    def prop_from_dummy( self ):
        """After function this atom has properties taken from
        neighbouring hydrogens which replace heavy previous atoms,
        
        and also the bonds between them"""
        p = Property()
        for bond in self.bonds:
            if bond._Atom2.is_dummy():
                p += bond.p + bond._Atom2.p
        return p


    @property
    def res_name(self):
        if self._res_name:
            return self._res_name
        return "XXX"
    @res_name.setter
    def res_name(self,val):
        self._res_name = val

#Atom
    @property
    def chain_id(self):
        if self._chain_id:
            return self._chain_id
        if self.Molecule:
            if self.Molecule.Chain:
                return self.Molecule.Chain.chain_id
        return None

class Residue( Molecule ):
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

        self.label_dict = {}

        super( Residue, self ).__init__( *args, **kwargs )

        self.concap = None
        self.ready = None
        self.is_concap = None
        self.is_ready = None
        self.is_bridge = None

        self.c_term = False
        self.n_term = False

#NexResidue setters
    @property
    def Prev(self):
        return self._Prev
    @property
    def Next(self):
        return self._Next
    @Prev.setter
    def Prev(self, val):
        self._Prev = val
    @Next.setter
    def Next(self, val):
        self._Next = val



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

#Residue
    def __str__(self):
        base = "-".join( [self.chain_id, self.res_name + str(self.res_id)] )
        if self.is_concap:
            base += '-concap'
        if self.is_ready:
            base += '-ready'
        return base

#Residue
    def add_atom( self, atom):
        if isinstance( atom, MoleculesAtom ):
            if any( (atom.r == x.r).all() for x in self):
                return
            self.append( atom )
            self.label_dict[ atom.label ] = atom
            atom.Molecule = self

#Residue
    def get_relevant_residues(self):
        if self.n_term:
            return map(lambda x:x.copy(), [self.ready, self.Next.ready] )
        elif self.c_term:
            return map(lambda x:x.copy(), [self.ready, self.Prev.ready] )
        else:
            return map(lambda x:x.copy(), [self.ready, self.Next.ready, self.Prev.ready] )

    def get_relevant_concaps(self):
        if self.n_term:
            return map(lambda x:x.copy(), [self.concap] )
        elif self.c_term:
            return map(lambda x:x.copy(), [self.Prev.concap] )
        else:
            return map(lambda x:x.copy(), [self.concap, self.Prev.concap] )

#Residue
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

#Residue
    def copy( self ):
        """Copy Residue method, return new Residue where all atoms
        , their bonds and properties are copies
        

        Warning, all the bonds in this returned instance may
        have pointers to outside Molecules
        """
        new_res = Residue()

        for atom in self:
            new_atom = atom.copy()
            for b in atom.bonds:
                at1 = b._Atom1.copy()
                at2 = b._Atom2.copy()
                at1.Molecule = new_res
                at2.Molecule = new_res
                new_bond = Bond( at1, at2 )
                if at1.res_id == at2.res_id:
                    new_bond._Molecule = new_res
                else:
                    new_bond._Molecule = None

                new_bond.p = b.p.copy_property()
                new_atom.add_bond( new_bond )
            new_res.add_atom( new_atom )

        new_res.c_term = self.c_term
        new_res.n_term = self.n_term
        new_res._res_id = self._res_id
        new_res._res_name = self.res_name
        new_res.AA = self.AA
#Keep information if this is reay/concap/bridge
        new_res.is_ready = self.is_ready
        new_res.is_concap = self.is_concap
        new_res.is_bridge = self.is_bridge
        new_res._level = self._level


        new_res._chain_id = self.chain_id

        new_res.Next = self.Next
        new_res.Prev = self.Prev
        new_res.Bridge = self.Bridge

        return new_res



#Residue
    def mfcc_props(self):
        """After this functions all atoms here have final properties.
        New implementation using bond midpoint and also between residues
        """
        self.populate_bonds( cluster = 1 )
        for at in self:
            at.p = Property()

        relevant_centers = []
#First update carbons that are attached to fake hydrogens with new props
        for res in self.get_relevant_residues():
            for center in res.get_ats_and_bonds():
                if isinstance( center, MoleculesAtom ):
                    center.p += center.prop_from_dummy()
                relevant_centers.append( center )

        for con in self.get_relevant_concaps():
            for center in con.get_ats_and_bonds():
                if isinstance( center, MoleculesAtom ):
                    center.p += center.prop_from_dummy()
                relevant_centers.append( center )

#Loop over all atoms and se if they are in residue or concap 
#Save the points inbetween to put then om neighbors for later
        points_inbetween = []
        for center_1 in self.get_ats_and_bonds():
            tmp_p = Property()
            for center_2 in relevant_centers:
                if np.allclose( center_1.r, center_2.r):
#This bond between residues belongs to no specific molecule
                    if not center_2._Molecule:
                        if center_2._Atom1._Molecule.is_ready:
                            tmp_p += center_2.p
                        elif center_2._Atom1._Molecule.is_concap:
                            tmp_p -= center_2.p
                        points_inbetween.append( center_1 )
#These are atoms and bonds which exist in the final residue conf
                    else:
                        if center_2._Molecule.is_ready:
                            tmp_p += center_2.p
                        elif center_2._Molecule.is_concap:
                            tmp_p -= center_2.p
            center_1.p += tmp_p

#Hacky solution so far to only compare the x-coordinate,  WARNING
#Not the full vector point
        points = unique( points_inbetween, key = lambda x: x.r[0],
                get_original = True )
        for point in points:
            for at in [ point._Atom1, point._Atom2 ]:
                if at in self:
                    at.p = at.p + point.p/2.0
        for point in points:
            for at in [ point._Atom1, point._Atom2 ]:
                if at in self:
                    at.bonds.remove( point )
        self.LoProp = True

#Property of Residue
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
        if self[0]:
            tmp_id = self[0]._res_id
            for at in self:
                try:
                    assert tmp_id == at._res_id
                except AssertionError:
                    logging.error( "No _res_id in Residue and not all atoms in same residue")
            if tmp_id is not None:
                return tmp_id
        tmp_id = 0
        return tmp_id

    def get_atom(self, pdb_name):
        #try:
        for i in self:
            if i.pdb_name == pdb_name:
                return i.copy()

    @res_id.setter
    def res_id(self, val):
        self._res_id = val

#Residue
    @property
    def res_name(self):
        if self._res_name:
            return self._res_name
        return "YYY"
    @res_name.setter
    def res_name(self, val):
        self._res_name = val
#Residue
    @property
    def chain_id(self):
        if self._chain_id is not None:
            return self._chain_id
        if self.Cluster:
            return self.Cluster.chain_id
        if self[0]:
            tmp_ch = self[0].chain_id
            for at in self:
                try:
                    assert tmp_ch == at.chain_id
                except AssertionError:
                    logging.error( "No Chain object or _chain_id in Residue and not all atoms have same chain_id")
        else:
            tmp_ch = "X"
        return tmp_ch
    @chain_id.setter
    def chain_id(self, val ):
        self._chain_id = val

#Residue Method
    def copy_info(self):
        new = Residue()
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



#Residue Method
    def gather_ready( self, 
            residue = False, r = False,
            concap = False, c = False,
            bridge = False, b = False,
            level = 2 ):
        p = Pattern()
        if residue or r:
            tmp_residue = self.copy_info()
            tmp_residue.is_ready = True
            res_type = "res"
        elif concap or c:
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
class Chain( Cluster):
    """Will behaive like Cluster and rely everything on getters and setters to avoid bugs in overwriting properties"""
    def __init__(self, *args, **kwargs):
        self._snapshot = None
        self._time = None
        self._freq = None
        self._System = None


        super( Chain, self ).__init__( *args, **kwargs )
#Temporary hack to get the chain ID of residue or list of residues
        if args is not ():
            self.chain_id = args[0][0]._chain_id

# To define N- /C- terminals and set Next/ Prev attributes
    def connect_residues(self,):
        self[0].n_term = True
        self[-1].c_term = True
        for i, res in enumerate(self):
            if self[i].n_term:
                self[i]._Next = self[i + 1]

            elif self[i].c_term:
                self[i]._Prev = self[i - 1]

            else:
                self[i]._Next = self[i + 1]
                self[i]._Prev = self[i - 1]

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

#setter for Chain
    @property
    def chain_id(self):
        if self._chain_id is not None:
            return self._chain_id
        return None

    @chain_id.setter
    def chain_id(self, val):
        self._chain_id = val

class System( list ):
    """Can hold instances of Clusters, Molecules, Atoms in it
    used to separate trajectory configurations between different snapshots.

    Each trajectory snapshot thus holds a System class
    
    """
    def __init__(self, *args, **kwargs):
        self._snapshot = None
        self._time = None
        self._freq = None
        super(System, self).__init__()

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
    def gather_ready(self, level = 2):

        for each in self.residues:
            each.gather_ready( concap = True, level = level )
            each.gather_ready( residue = True, level = level )

    def connect_everything(self):
        for atom in self.atoms:
            atom.System = self
        for mol in self.molecules:
            mol.System = self
            mol.connect_everything()
        for cluster in [c for c in self if isinstance( c, Cluster) ]:
            cluster.System = self
            cluster.connect_everything()

    @property
    def residues(self):
        return [m for m in self.molecules if isinstance( m, Residue ) ]

    @property
    def coc(self):
        return sum( [at.r * charge_dict[at.element] for at in self.atoms])\
                /sum( map(float,[charge_dict[at.element] for at in self.atoms]) )
    @property
    def p(self):
        return self.sum_property
    @property
    def sum_property(self):
        """
Return the sum properties of all properties in System
        """
        conv = 1.0
        p = Property()
        coc = self.coc
        if self.AA:
            conv = 1/a0

        el_dip = np.array([ conv*(center.r- coc)*center.p.q for mol in self.molecules for center in mol.get_ats_and_bonds()  ])
        nuc_dip = np.array([ conv*(center.r-coc)*charge_dict[center.element] for mol in self.molecules for at in mol.get_ats_and_bonds() ])
        dip_lop = np.array([ center.p.d for mol in self.molecules for center in mol.get_ats_and_bonds() ])
        dip = el_dip + nuc_dip
        d = (dip + dip_lop).sum(axis=0)
        p = Property()

        for center in [c for mol in self.molecules for c in mol.get_ats_and_bonds()]:
            p.q += center.p.q
            p.a += center.p.a
            p.b += center.p.b
        return p

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

# Slice solution for System, will return System when slicing object
    def __add__(self, other):
        return System(list.__add__(self, other))
    def __getslice__(self, i, j):
        return self.__getitem__(slice(i, j))
    def __getitem__(self, item):
        if isinstance( item, slice ):
            result = list.__getitem__(self, item)
            try:
                return System(result)
            except TypeError:
                return result
        else:
            return super(System,self).__getitem__( item )
#System
    def get_xyz_string(self):
        _string = "%d\n\n" %( len(self.atoms) )
        for at in self.atoms:
            _string += at.xyz_string()
        return _string

    def min_dist_atoms_separate_res_chain(self, AA_cutoff = 1.5 ):
        """Return list of atoms which have an other atom closer than 1.5 AA to them
        and are not in the same residue and also different chains
        
        """
        tmp = []
        if not self.AA:
            AA_cutoff /= a0
        ats = self.atoms
        for i1, at1 in enumerate( ats[:-1] ):
            for i2, at2 in enumerate( ats[i1:]):
                if (at1.chain_id == at2.chain_id) and (at1.res_id == at2.res_id):
                    continue
                if at1.dist_to_atom ( at2 ) < AA_cutoff:
                    tmp.append( at1 )
                    tmp.append( at2 )
        return unique( tmp )


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

    @property
    def AA(self):
        AA = self[0].AA
        for chain in self:
            try:
                AA == chain.AA
            except AssertionError:
                logging.error("All chains in system %s are not of same unit" %self)
        return AA

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
        world = System()

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
                tmp_chain.add ( tmp_residue )
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
                tmp_chain.add( tmp_residue )
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
                    if "HG1" in [x.pdb_name for x in res]:
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


#System method of adding objects
    def add(self, item):
        if isinstance( item, Cluster):
            self.append( item )
            item.System = self
        elif isinstance( item, Molecule):
            c = Chain( item )
            c.System = self
            self.append( c )

#System method of adding objects
    @property
    def atoms(self):
        return [a for chain in self for mol in chain for a in mol if isinstance(a , MoleculesAtom ) ]
    @property
    def molecules(self):
        return [m for chain in self for m in chain if isinstance(m , Molecule ) ]
    @classmethod
    def from_pdb_string( cls, _string,
            connect = False,
            detect_term = False ):
        """Assuming the string is a pdb format, read in all chains and stuff"""
        res, meta = all_residues_from_pdb_string( _string )

        nc =  Chain( splitter(res, lambda x: x.chain_id)[0] )
        chains = [ Chain(ch) for ch in splitter( res, lambda x: x.chain_id )]

        system = cls( *chains )
        system.meta = meta

        if connect:
            for ch in chains:
                ch.connect_residues()
        if detect_term:
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
        if not isinstance( pick, System ):
            raise TypeError("Wrong pickled class")
        return pick
    
        
class World( list ):
    """Can hold instances of Systems, Clusters, Molecules, Atoms in it"""
    def __init__(self, *args, **kwargs):
        super(World, self).__init__()
        self._project = None

        if args is not None:
            for each in args:
                if isinstance( each, System):
                    self.append( each )
        
    def connect(self):
        for s in [s for s in self if isinstance(s, System) ]:
            s.connect_everything

    def add(self, item):
        if isinstance(item, System ):
            self.append( item )
 
    def save(self, fname = "world.p"):
        pickle.dump( self, open(fname, 'wb' ), protocol = 2 )

    
    @property
    def molecules(self):
        return [m for s in self for ch in s for m in ch if isinstance( m, Molecule ) ]

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
    Syst = System.from_pdb_string( open(args.pdbfile).read() )
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
