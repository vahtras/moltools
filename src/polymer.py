import pdbreader
import molecules, re, copy
import copy as copy
import numpy as np
import utilz
import logging

charge_dict = {"H": 1.0, "C": 6.0, "N": 7.0, "O": 8.0, "S": 16.0,
        "P" : 15, "X" : 0.0 }

class NewPattern( pdbreader.Pattern ):

    def __init__(self, *args, **kwargs):
        super( NewPattern, self).__init__(*args, **kwargs)
#Ready level 3 to add
        self[ ( 'reg', 'res', 'p', 'add', 3 ) ] = re.compile(r'H1$|H2$|C$|CA$|HN$')
        self[ ( 'reg', 'res', 't', 'add', 3 ) ] = re.compile(r'.*')
        self[ ( 'reg', 'res', 'n', 'add', 3 ) ] = re.compile(r'H1$|H2$|C$|CA$|HC$')

#Concaps level 3 to add
        self[ ( 'reg', 'con', 'p', 'add', 3 ) ] = re.compile(r'DUMMY')
        self[ ( 'reg', 'con', 't', 'add', 3 ) ] = re.compile(r'CA$|C$|H1$|H2$|HN$')
        self[ ( 'reg', 'con', 'n', 'add', 3 ) ] = re.compile(r'C$|CA$|H1$|H2$|HC$')

#Residues level 3 to connect
        self[ ( 'reg', 'res', 'pp', 'con', 3 ) ] = { "CB1" :["CA"], "CB2" : ["CA"] }
        self[ ( 'reg', 'res', 'tt', 'con', 3 ) ] = { }
        self[ ( 'reg', 'res', 'nn', 'con', 3 ) ] = { "CB1" : ["CA"], "CB2" : ["CA"] }
        self[ ( 'reg', 'res', 'pp_p', 'con', 3 ) ] = { "CA" : ["C"] }
        self[ ( 'reg', 'res', 'nn_n', 'con', 3 ) ] = { "C" : ["CA"] }

#Concap level 3 to connect
        self[ ( 'reg', 'con', 'pp', 'con', 3 ) ] = { }
        self[ ( 'reg', 'con', 'p_t', 'con', 3 ) ] = { "CA": ["C"] }
        self[ ( 'reg', 'con', 'tt', 'con', 3 ) ] = { "CB1" : ["CA"], "CB2" :["CA"] }

        self[ ( 'reg', 'con', 'nn', 'con', 3 ) ] = { "CB1" : ["CA"], "CB2" : ["CA"]  }
        self[ ( 'reg', 'con', 'nn_n', 'con', 3 ) ] = { "C" : ["CA"] }

#Custom heavy to transform to hydrogen and scale them in distance
#Residues level 3
        self[ ( 'reg', 'res', 'pp', 'rep', 3 ) ]  = re.compile(r'CB2$|CB1$')
        self[ ( 'reg', 'res', 'tt', 'rep', 3 ) ] = re.compile(r'BALLONY')
        self[ ( 'reg', 'res', 'nn', 'rep', 3 ) ]  = re.compile(r'CB2$|CB1$')
        self[ ( 'reg', 'res', 'pp_p', 'rep', 3 ) ]  = re.compile(r'CA$')
        self[ ( 'reg', 'res', 'nn_n', 'rep', 3 ) ]  = re.compile(r'C$')

#Concap level 3
        self[ ( 'reg', 'con', 'p_t', 'rep', 3 ) ]  = re.compile(r'CA$')
        self[ ( 'reg', 'con', 'tt', 'rep', 3 ) ] = re.compile(r'CB1$|CB2$')
        self[ ( 'reg', 'con', 'nn', 'rep', 3 ) ]  = re.compile(r'CB2$|CB1$')
        self[ ( 'reg', 'con', 'nn_n', 'rep', 3 ) ]  = re.compile(r'C$')



class Atom( pdbreader.Atom ):

    def __init__(self, *args, **kwargs):
        super(Atom, self).__init__( *args, **kwargs)
        self._label = None

    #@property
    #def label(self):
    #    if self._label is not None:
    #        return self._label
    #    return "-".join( [ str(self.Molecule.res_id), self.Molecule.res_name, self.pdb_name ] )

class Monomer( pdbreader.Residue ):
    """Right now specific for pmma but late can be made general to any
    monomer type"""



    def __init__(self, *args, **kwargs ):
        super( Monomer, self ).__init__( *args, **kwargs )
        self.c_term = False
        self.n_term = False
        self._label = None
        self.is_ready = False
        self.Cluster = Polymer()
        """Standard connectivities"""
        self._r = 1.0
        self._angle = 104.5
        self._dihedral = 0.0
        self._rot_angle = 0.0
        self.hidden = []


    def __str__(self):
        st =  "-".join( [str(self.res_id), self.res_name] )
        if self.concap:
            st += '-con'
        elif self.is_ready:
            st += '-ready'

        return st

    def copy_self(self):
        return copy.deepcopy(self)

    def copy_info(self):
        new = Monomer()
        new.c_term = self.c_term
        new.n_term = self.n_term
        new.AA = self.AA
        new.in_qm_region = self.in_qm_region
        new.in_mm_region = self.in_mm_region
        new._res_id = self._res_id
        new._res_name = self.res_name
        new.Next = self.Next
        new.Prev = self.Prev
        new.Chain = self.Chain
        return new

    def gather_ready( self, 
            residue = False, r = False,
            concap = False, c = False,
            bridge = False, b = False,
            level = 1 ):
        """Override the pdbreader.Residue.gather_ready method due to new pattern"""
        p = NewPattern()
        if residue or r:
            tmp_residue = self.copy_info()
            res_type = "res"
            tmp_residue.is_ready = True
        elif concap or c:
            tmp_residue = self.copy_info()
            tmp_residue = self.copy_info() 
            tmp_residue.concap = True
            res_type = "con"
        else:
            print "Give some option to what to make"
            return 

        p_rep_pp_p = p.get()
        
#To set defaults for pattern for this residue
        res_name = "reg"

        p_con_pp_p = {}
        p_rep_pp_p = re.compile( r'(?!x)x' )
# Patterns to just add atoms, previos depends on type of previous
        if self.Prev:
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
            at_t = pdbreader.get_matching( self, p_add_t )
            at_n = pdbreader.get_matching( self.Next, p_add_n )
            at_r1 = pdbreader.get_replacement( self, p_rep_tt, p_con_tt, )
            at_r2 = pdbreader.get_replacement( self.Next, p_rep_nn, p_con_nn, )
            if self.Next.Next:
                at_r3 = pdbreader.get_rep_2( self.Next, self.Next.Next, p_rep_nn_n, p_con_nn_n )
            for at in at_t + at_n + at_r1 + at_r2 + at_r3:
                tmp_residue.add_atom( at )

        elif self.c_term:
            if concap or c:
                return 
            at_r2 = []
            at_p = pdbreader.get_matching( self.Prev, p_add_p )
            at_t = pdbreader.get_matching( self, p_add_t )
            at_r1 = pdbreader.get_replacement( self.Prev, p_rep_pp, p_con_pp, )
            if self.Prev.Prev:
                at_r2 = pdbreader.get_rep_2( self.Prev, self.Prev.Prev, p_rep_pp_p, p_con_pp_p )
            for at in at_p + at_t + at_r1 + at_r2 :
                tmp_residue.add_atom( at )
        else:
            at_r1, at_r2 = [], []
            if concap or c:
                at_r1 = pdbreader.get_replacement( self, p_rep_tt, p_con_tt, )
                at_r2 = pdbreader.get_rep_2( self, self.Prev, p_rep_p_t, p_con_p_t, )
            else:
                at_r1 = pdbreader.get_replacement( self.Prev, p_rep_pp, p_con_pp, )

            at_r4, at_r5 = [], []
            at_p = pdbreader.get_matching( self.Prev, p_add_p )
            at_t = pdbreader.get_matching( self, p_add_t )
            at_n = pdbreader.get_matching( self.Next, p_add_n )
            at_r3 = pdbreader.get_replacement( self.Next, p_rep_nn, p_con_nn, )
            if self.Next.Next:
                at_r4 = pdbreader.get_rep_2( self.Next, self.Next.Next, p_rep_nn_n, p_con_nn_n )
            if self.Prev.Prev:
                at_r5 = pdbreader.get_rep_2( self.Prev, self.Prev.Prev, p_rep_pp_p, p_con_pp_p )

            for at in at_p + at_t + at_n + at_r1 + at_r2 + at_r3 + at_r4 + at_r5:
                tmp_residue.add_atom( at )

        if residue or r:
            self.ready = tmp_residue
        if concap or c:
            self.con = tmp_residue
        if bridge or b:
            self.bri = tmp_residue

    @property
    def HN(self):
        return self.first_h
    @property
    def HC(self):
        return self.last_h
    @property
    def CX(self):
        return self.get_atom_by_pdb_name( 'CX' )

    @property
    def first_h(self):
        for at in self.hidden:
            if at.pdb_name == 'HN':
                return at
        return self.get_atom_by_pdb_name( 'HN' )

    @property
    def last_h(self):
        for at in self.hidden:
            if at.pdb_name == 'HC':
                return at
        return self.get_atom_by_pdb_name( 'HC' )


    @property
    def last_heavy(self):
        return self.get_atom_by_pdb_name( 'CA' )

    @property
    def first_heavy(self):
        return self.get_atom_by_pdb_name( 'C' )

    def get_atom_by_pdb_name(self, label, dup = False):
        at = []
        for i in self:
            if i.pdb_name == label:
                at.append(i)
        if len(at) > 1 and not dup:
            print "Warning: Duplicate with pdb name %s, returning %s" %(label, at[0].label )
            return at[0]
        elif len(at) > 1 and dup:
            return at

        elif len(at) == 0:
            print "No %s in %s" %(label, self)
            return
        return at[0]

    def get_mol_string(self, basis = ("ano-1 2", "ano-1 4 3 1",
        "ano-2 5 4 1" ) ):
        if len( basis ) > 1:
            el_to_rowind = {"H" : 0, "C" : 1, "O" : 1, "N" : 1,
                    "S" : 2, "P" : 2}
        else:
            el_to_rowind = {"H" : 0, "C" : 0, "O" : 0, "N" : 0, "S" : 0 }
        st = ""
        s_ = ""
        if self.AA: s_ += " Angstrom"
        uni = utilz.unique([ at.element for at in self])
        st += "ATOMBASIS\n\n\nAtomtypes=%d Charge=0 Nosymm%s\n" %(len(uni), s_)
        for el in uni:
            st += "Charge=%s Atoms=%d Basis=%s\n" %( str(charge_dict[el]),
                    len( [all_el for all_el in self if (all_el.element == el)] ),
                    basis[ el_to_rowind[el] ])
            for i in [all_el for all_el in self if (all_el.element == el) ]:
                st += "%s %.5f %.5f %.5f\n" %(i.element, i.x, i.y, i.z ) 
        return st

    def get_mol(self, basis = ['ano-1 2','ano-1 4 3 1', 'ano-2 5 4 1' ],
            level = 0):
        """get the Residue.get_mol, but set charge to 0"""
        charge = 0

        elem = [x.element for x in self]
        atomtypes = len( utilz.unique( elem ) )

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



    @staticmethod
    def from_pdb( f_, in_AA= True, out_AA=False):
#Create an Atom here and add it to the residue
        pat = re.compile(r'^HETATM|^ATOM')
        pat_element = re.compile( r'([A-Z])' )
        
        tmp_mono = Monomer( AA = in_AA )
        for i, line in enumerate(open(f_).readlines()):
            if not pat.match( line ):continue
            res_name = line[17:21].strip()
            try:
                res_id = int( line[22:26].strip() )
            except ValueError:
                pass
#Initiate values for Atom,
            chain_id = line[21:22].strip()
            x = line[30:38].strip()
            y = line[38:46].strip()
            z = line[46:54].strip()
            pdb_name = line[11:16].strip()
            try:
                element = pat_element.search( pdb_name ).group(1)
#This will cause IndexError for TER when no other entries are there
            except IndexError:
                element = None
            except AttributeError:
                element = None
            tmp_atom = Atom()
            tmp_atom.Molecule = tmp_mono
            tmp_atom.x, tmp_atom.y, tmp_atom.z = map( float, [x, y, z] )
            tmp_atom.element = element
            tmp_atom.pdb_name = pdb_name
            tmp_atom.res_name = res_name
            tmp_atom._res_id = res_id
            #tmp_atom.label = "%d-%s-%s" %( res_id, res_name, pdb_name )
            tmp_mono._res_name = res_name
            tmp_mono._chain_id = chain_id
            tmp_mono.add_atom( tmp_atom )
        if in_AA and not out_AA:
            tmp_mono.to_AU()
        return tmp_mono

    def hide(self, atom):
        self.hidden.append( atom )
        self.remove( atom )

    def unhide(self, atom):
        self.append( atom )
        self.hidden.remove( atom )

    def unhide_all(self):
        while len( self.hidden ) > 0:
            at = self.hidden[0]
            self.unhide( at )

    @property
    def with_hidden(self):
        cpy = copy.deepcopy(self)
        for at in self.hidden:
            cpy.append(at)
        return cpy

    def surface(self):
        for at in self.hidden:
            self.append( at )
        self.hidden = []
        return self

    @property
    def res_id(self):
        if self._res_id is not None:
            return self._res_id
        else:
            return 1

    def rotate_around(self, p1, p2, theta = 0.0):
        """Rotate clockwise around line formed from point p1 to point p2 by theta"""
        for at in self + self.hidden:
            at.x, at.y, at.z = utilz.rotate_point_by_two_points( at.r, p1, p2, theta)

    def reflect_around(self, p1, p2, p3):
        """Rotate clockwise around line formed from point p1 to point p2 by theta"""
        for at in self + self.hidden:
            at.x, at.y, at.z = utilz.reflect_point_by_three_points( at.r, p1, p2, p3 )

    def translate_by_r(self, r):
        """Need to override class.Molecules and also include hidden hydrogens"""
        super( Monomer, self ).translate_by_r( r )
        for at in self.hidden:
            at.x += r[0]
            at.y += r[1]
            at.z += r[2]

    def inv_rotate(self, t1, t2, t3):
        """Need to override class.Molecules and also include hidden hydrogens"""
        r1 = utilz.Rz_inv(t1)
        r2 = utilz.Ry(t2)
        r3 = utilz.Rz_inv(t3)
        for at in self + self.hidden:
            at.x, at.y, at.z = np.einsum('ab,bc,cd,d', r3, r2, r1, at.r )
            at.p = at.p.inv_rotate( t1, t2, t3 )

    
    def connect_monomer(self, other, last = True):
        """By default connects to the last point of the Monomer"""
        if last:
            r_ca = self.last_heavy.r
            e_old = (self.last_h.r - r_ca)/ np.linalg.norm( (self.last_h.r - r_ca) )
            r_new = e_old * other._rn
            other.translate_by_r( r_ca + r_new - other.first_heavy.r )

            self.c_term = False
            other.c_term = True
            other._res_id = self._res_id + 1

            other.hide( other.first_h )
            self.hide( self.last_h )
        else:
            r_ca = self.last_heavy.r
            e_old = (self.first_heavy.r - first_h.r)/ np.linalg.norm( (self.first_heavy.r - self.first_h.r) )
            r_new = e_old * other._rn
            other.translate_by_r( r_ca + r_new - other.first_heavy.r )

            self.c_term = False
            other.c_term = True
            other._res_id = self._res_id + 1

            other.hide( other.first_h )
            self.hide( self.last_h )

class Polymer( molecules.Cluster ):
    @staticmethod
    def from_monomer( mono ):
        """Initiate polymer by calling it with a monomer as argument"""
        P = Polymer()
        P.chain_id = "A"
        mono.c_term = True
        mono.n_term = True
        mono._res_id = 1
        mono.Cluster = P
        mono.Polymer = P
        P.append( mono )
        return P

    def add_monomer(self, mono, last = True):
        mono = copy.deepcopy( mono )
        last = self[-1]
        last.Next = mono
        mono.Prev = last
        mono.Cluster = self
        mono.Polymer = self
        last.connect_monomer( mono, last = last )
        self.append( mono )

    @property
    def with_hidden(self):
        """Return self but with all hidden hydrogens which are on N and C
        terminal ends included"""

        p = Polymer.from_monomer( copy.deepcopy( self[0].with_hidden ) )

        for r in self[1:-1]:
            p.add_monomer( copy.deepcopy( r.with_hidden ) )
        p.add_monomer( copy.deepcopy( self[-1].with_hidden) )

        return p

    def get_qm_mol_string(self, basis = ("ano-1 2 1", "ano-1 3 2 1", "ano-2 5 4 1" ), ):
        """Override Cluster method to include hidden hydrogens"""

        if len( basis ) > 1:
            # Set row index number to periodic table one
            el_to_rowind = {"H" : 0, "C" : 1, "O" : 1, "N" : 1 , "S" : 2  }
        else:
            # Keep all 0, since basis is only len 1
            el_to_rowind = {"H" : 0, "C" : 0, "O" : 0, "N" : 0, "S" : 0 }

        st = ""
        comm1 = "QM: " + " ".join( [ str(m) for m in self if m.in_qm] )[:72]
        comm2 = "MM: " + " ".join( [ str(m) for m in self if m.in_mm] )[:73]
        uni = utilz.unique([ at.element for mol in self for at in mol if mol.in_qm])
        s_ = ""
        if self.AA: s_ += "Angstrom"

        st += "ATOMBASIS\n%s\n%s\nAtomtypes=%d Charge=0 Nosymm %s\n" %( \
                comm1,
                comm2,
                len(uni),
                s_)
        for el in uni:
            st += "Charge=%s Atoms=%d Basis=%s\n" %( str(charge_dict[el]),
                    len( [all_el for mol in self for all_el in mol if ((all_el.element == el) and mol.in_qm )] ),
                     basis[ el_to_rowind[el] ] )
            for i in [all_el for mol in self for all_el in mol if ((all_el.element == el) and mol.in_qm) ]:
                st += "{0:5s}{1:10.5f}{2:10.5f}{3:10.5f}\n".format( i.element, i.x, i.y, i.z )
        return st




