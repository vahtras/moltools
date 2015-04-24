import molecules, re

class Monomer( molecules.Molecule ):

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
            tmp_atom = molecules.Atom()
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
