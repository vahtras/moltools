import pdbreader
import molecules, re, copy
import numpy as np
import utilz

class NewPattern( pdbreader.Pattern ):

    def __init__(self, *args, **kwargs):
        super( NewPattern, self).__init__(*args, **kwargs)
#ready level 3 to add
        self[ ( 'reg', 'res', 'p', 'add', 3 ) ] = re.compile(r'CA$|CB\d{1}$|HB\d{1}$')
        self[ ( 'reg', 'res', 't', 'add', 3 ) ] = re.compile(r'.*')
        self[ ( 'reg', 'res', 'n', 'add', 3 ) ] = re.compile(r'C$|CA$|CB\d{1}$|HB\d{1}$')

#Concaps level 3 to add
        self[ ( 'reg', 'con', 'p', 'add', 3 ) ] = re.compile(r'DUMMY')
        self[ ( 'reg', 'con', 't', 'add', 3 ) ] = re.compile(r'CA$|CB\d{1}$|HB\d{1}$')
        self[ ( 'reg', 'con', 'n', 'add', 3 ) ] = re.compile(r'C$|CA$|CB\d{1}$|HB\d{1}$')

#Residues level 3 to connect
        self[ ( 'reg', 'res', 'pp', 'con', 3 ) ] = { "OG1" : ["CB2"], 
                "OG2" : ["CB2"], "C" : ["CA"] }
        self[ ( 'reg', 'res', 'tt', 'con', 3 ) ] = { }
        self[ ( 'reg', 'res', 'nn', 'con', 3 ) ] = { "OG1" : ["CB2"], 
                "OG2" : ["CB2"], "C" : ["CA"] }
#Concap level 3 to connect
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
        self[ ( 'reg', 'res', 'pp', 'rep', 3 ) ]  = re.compile(r'OG\d{1}$|C$')
        self[ ( 'reg', 'res', 'tt', 'rep', 3 ) ] = re.compile(r'BALLONY')
        self[ ( 'reg', 'res', 'nn', 'rep', 3 ) ]  = re.compile(r'OG\d{1}$|C$')


class Monomer( pdbreader.Residue ):
    """Right now specific for pmma but late can be made general to any
    monomer type"""


    def __init__(self, *args, **kwargs ):
        super( Monomer, self ).__init__( *args, **kwargs )
        self.c_term = False
        self.n_term = False
        """Standard connectivities"""
        self._r = 1.0
        self._angle = 104.5
        self._dihedral = 180.0
        self._rot_angle = np.pi
        self.hidden = []

    @property
    def first_h(self):
        return self.get_atom_by_pdb_name( 'HN' )
    @property
    def last_h(self):
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
            tmp_atom = pdbreader.Atom()
            tmp_atom.Molecule = tmp_mono
            tmp_atom.x, tmp_atom.y, tmp_atom.z = map( float, [x, y, z] )
            tmp_atom.element = element
            tmp_atom.pdb_name = pdb_name
            tmp_atom.res_name = res_name
            tmp_atom._res_id = res_id
            tmp_atom.label = "%d-%s-%s" %( res_id, res_name, pdb_name )
            tmp_mono.add_atom( tmp_atom )
        if in_AA and not out_AA:
            tmp_mono.to_AU()
        return tmp_mono

    def hide(self, atom):
        self.hidden.append( atom )
        self.remove( atom )

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

    def connect_monomer(self, other):
        r_ca = self.last_heavy.r
        e_old = (self.last_h.r - r_ca)/ np.linalg.norm( (self.last_h.r - r_ca) )
        r_new = e_old * other._r
        other.translate_by_r( r_ca + r_new - other.first_heavy.r )

        self.c_term = False
        other.c_term = True
        other._res_id = self._res_id + 1

        other.hide( other.first_h )
        self.hide( self.last_h )

        p1 = self.last_heavy.r
        p2 = other.first_heavy.r
        other.rotate_around( p1, p2, other.res_id * other._rot_angle )



class Polymer( molecules.Cluster ):
    @staticmethod
    def from_monomer( mono ):
        """Initiate polymer by calling it with a monomer as argument"""
        P = Polymer()
        P.chain_id = "A"
        mono.c_term = True
        mono.n_term = True
        mono._res_id = 1
        P.append( mono )
        return P

    def add_monomer(self,mono):
        mono = copy.deepcopy( mono )
        last = self[-1]
        last.Next = mono
        mono.Prev = last
        last.connect_monomer( mono )
        self.append( mono )


