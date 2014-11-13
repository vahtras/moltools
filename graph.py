#!/usr/bin/env python

#from functions import print_f, fun
from sys import argv
import re
import numpy as np

from molecules import Cluster



def run_argparse():
    import argparse

    A1 = argparse.ArgumentParser()

    A1.add_argument('-d', action="store_true", default=False)
    A1.add_argument('-qm', dest='qm',
            default= 'qm.xyz',
            help = 'file to read qm region, default [qm.xyz]')
    A1.add_argument('-mm', dest='mm',
            default= 'mm.xyz',
            help = 'file to read mm region, default [mm.xyz]')
    A1.add_argument('-outfile', 
            default= "new_mm",
            help = "File to write the new mm_region to, default new_mmX.pro, X is cutoff")
    A1.add_argument('-region_cutoff', dest='region_cutoff',
            type = float,
            default= 5.0,
            help = 'Cutoff for which atoms to include in qm/mm border, default 15 A')
    A1.add_argument('-cutoff', dest='cutoff',
            type = float,
            default= 1.2,
            help = 'Cutoff for which atoms to remove from qm/mm border, default 1.2 A')


    A1.add_argument('-xyz_files', action = 'store_true',
            default = False )
    A1.add_argument('-groups', action = 'store_true',
            default = False )
    A1.add_argument('-waters', type = float)
##############How many columns to read from mm_region depends on these numbers

    A1.add_argument('-angular', type = str,
            choices = ['0', '1', '2', '3'],
            default = '0',
            help = 'multipole moment number, default 0')
    A1.add_argument('-polar', type = str,
            choices = ['0', '1', '2'],
            default = '1',
            help = 'polarizability number, default 1')
    A1.add_argument('-hyper', type = str,
            choices = ['0', '1', '2'],
            default = '1',
            help = 'polarizability number, default 1')
    A1.add_argument('--charge', type = str,
            default = '0',
            choices = ['0', '1'],)
    A1.add_argument('--angular', type = str,
            choices = ['0', '1', '2', '3'],
            default = '0',
            help = 'multipole moment number, default 0')
    A1.add_argument('--polar', type = str,
            choices = ['0', '1', '2'],
            default = '0')
    args = A1.parse_args( argv[1:])



    
    return args




class Props(object):
    def __init__(self):
        """
        class that holds loprop properties, update later for full compatibility
        Right now only point charge, point polarizability
        """
        charge = 0.0
        angular = 0.0
        polar = 0.0

    def __str__(self):
        return " ".join(map( str , [self.charge, self.polar ] ))

class Atom:
    def __init__(self, line, args, qm = False):

        """
        Input is a line from an atom entry that has a label, x, y, z, and optional
        charge, dipoles, and polarizabilites

        args holds definitions of quantum numbers

-------------------------------------------------------------

        label is a string of which atom it is, unique:
        >>> atype
        '63H-V-HA'

        x, y, z are floats represtenting coordinates
        >>> x
        1.0

        props is a class with mm float point properties such as charge, multipole etc.
        >>> props.Charge
        -0.213
        >>> props.Polar
        8.213

        CloseAtoms is a list of atoms that are close to current atom.
        at first is empty but will be filled later on by class Cluster.
        """
        self.CloseAtoms = []
        self.Props = Props()

        l = int(real_l(args.angular))
        a = int(real_a(args.polar))

        line = line.split()

#Exit if given q, l and a don't agree with supplied MM file
        if (4 + l + a) != len(line):
            print "Wrong supplied combination of --charge, --angular and --polar"

        if qm:
            return

#If q = l = a = 0, read only xyz
        n = False

#This is xyz 
        if ( l+a ) == 0:
            self.label = line[0]
            self.x, self.y, self.z = map(float, [line[1], line[2], line[3]])

#This is the  charges
        if l == 0:
            self.label = line[0]
            self.x, self.y, self.z = map(float, [line[1], line[2], line[3]])
            self.Props.Charge = float( line[4] )

#This is charges and polarizabilities
        if a == 0:
            self.label = line[0]
            self.x, self.y, self.z = map(float, [line[1], line[2], line[3]])
            self.Props.Charge = float( line[4] )
            self.Props.Polar = float( line[5] )
#
##Charges and multipoles
#        elif (a) == 0:
#            n = Atom( atom = i[0],
#                    xyz = [i[1], i[2], i[3]],
#                    charge = i[4],
#                    angular = i[5:5+l+1]
#                    )
##All included
#        else:
#            n = Atom( atom = i[0],
#                    xyz = [i[1], i[2], i[3]],
#                    charge = i[4],
#                    angular = i[5:5+l+1],
#                    polar = i[5+l+1:5+l+a+1]
#                    )
        try:
            self.atype = re.compile(r'^\d+([A-Z]{1})').match( line[0] ).group(1)
        except AttributeError:
            self.atype = line[0]


    def __eq__(self, other):
        if (self.x == other.x) and (self.y == other.y) and (self.z == other.z):
            return True
        else:
            return False

    def __str__(self):
        return self.label +' '+ str(self.x) + str(self.y) + str(self.z) +' '+str(self.Props.Charge) + ' CloseAtoms: %s'%" ".join(map(lambda x: \
                x.label, self.CloseAtoms))

    def Same(self,other):
        """
        Returns true if coordinates of self and other Matches

        Threshhold of 0.1 Angstrom
        """
        if self.Dist(other) < 0.1:
            return True
        else:
            return False

    def HasCloseAtom(self,other):
        """
        returns true if self and other are CloseAtoms.
        """
        if other in self.CloseAtoms:
            return True
        else: 
            return False

    def AddCloseAtom(self, other):
        """
        Check if neighbour is not the same
        """
        if self != other:
            self.CloseAtoms.append(other)

    def RemoveCloseAtom(self, other):
        """
        Check if neighbour is there, then remove 
        """
        if other in self.CloseAtoms:
            self.CloseAtoms.remove(other)

    def Dist(self,other):
        r=np.sqrt((self.x - other.x)**2 +
                (self.y - other.y )**2 +
                (self.z - other.z )**2 )
        return r

class Cluster:

    atom_index = {'H':0, 'C':1, 'N':2, 'O':3, 'S':4}
    Dist_matrix = [
            [ 0.0, 1.2 , 1.2 , 1.2, 1.2 ],
            [ 1.2, 1.6 , 1.75, 1.7, 1.7 ],
            [ 1.2, 1.75, 0.0 , 1.7, 2.5 ],
            [ 1.2, 1.7 , 1.7 , 0.0, 2.5 ],
            [ 1.2, 1.7 , 2.5 , 2.5, 0.0 ]]
            

    def __init__(self, label):
        self.label = label
        self.Atomlist = []

    def __len__(self):
        return len(self.Atomlist)

    def __getitem__(self, val):
        return self.Atomlist[val]

    def AddAtom(self, atom):
        """
        If atom doesn't exist, then add it
        """
        if atom in self:
            return
        self.Atomlist.append( atom )

    def ExistsAtom(self, atom):
        """
        Search Atomlist for atom, if it exists, return True
        """
        for i in self.Atomlist:
            if i.Same(atom):
                return True
        return False

    def UpdateCloseAtoms(self, MM_C ):
        for i in self.Atomlist:
            for j in self.Atomlist:
                if j in i.CloseAtoms:
                    continue

                r = self.Dist_matrix[ self.atom_index[ i.atype ]][ self.atom_index[ j.atype ]]
                if i.Dist(j) < r and j in MM_C:
                    i.AddCloseAtom(j)

    def RemoveAtom(self, atom, args):
        """
        Removes the Atom that has coordinates Matching x, y and z,
        but before removing, passes on it's properties evenly to
        CloseAtoms, and removes itself from their neighbour lists.

        If this is connected to a hydrogen, first take hydrogen properties to itself, and remove the hydrogen, that way, group properties are concerted isntead of atom properties
        """
#Check if Atom exists in Atomlist
        if not self.ExistsAtom( atom ):
            return

        if len( atom.CloseAtoms) == 0:
            self.Atomlist.remove(atom)
            return
        try:
            q = atom.Props.Charge / len(atom.CloseAtoms)
            a = atom.Props.Polar / len(atom.CloseAtoms)
        except ZeroDivisionError:
            print "ERROR Tried to remove atom with no CloseAtoms, try increasing region cutoff!"
            exit(10)

        if args.groups:
            dic = {}
            d = []
            for i in atom.CloseAtoms:
                r = atom.Dist(i)
                dic[r] = i
                d.append( r )
            keep = dic[ max(d) ]
            atom.RemoveCloseAtom( keep )

            tmp = []
            for i in atom.CloseAtoms:
                atom.Props.Charge += i.Props.Charge
                atom.Props.Polar += i.Props.Polar
                tmp.append( i )
            for i in tmp:
                self.Atomlist.remove( i )
                atom.RemoveCloseAtom( i )
            keep.Props.Charge += atom.Props.Charge
            keep.Props.Polar += atom.Props.Polar
            self.Atomlist.remove( atom )

            
        else:
            tmp = []
            for i in atom.CloseAtoms:
                i.Props.Charge += q
                i.Props.Polar += a
                tmp.append(i)
            for i in tmp:
                atom.RemoveCloseAtom( i )
            self.Atomlist.remove( atom )


def real_l(l):
    """
    Returns number of angular momentum numbers to read depending on the number
    >>>real_l(2)    %Quadrupoles
    3
    """
    l = str(l)
    if l == '0':
        return '0'
    elif l == '1':
        return '3'
    elif l == '2':
        return '6'
    elif l == '3':
        return '9'
    print 'BUG discovered, angular momentum number is out of bound'
    raise ValueError
def real_a(a):
    """
    Returns number of polarizability numbers to read depending on the number
    >>>real_a(1)    %isotropic
    1
    """
    a = str(a)
    if a == '0':
        return '0'
    elif a == '1':
        return '1'
    elif a == '2':
        return '6'
    print 'BUG discovered, polarizability number is out of bound'
    raise ValueError

def main():
    """
    Main function for given input of MM and QM file, along with cutoff threshhold,
    reorders the MM file to only include mm sites close to QM within a given distance.

    Purpose is to study the automatic variation of the QM size given a fixed pdb
    crystal structure. For covalently bonded systems in QM/MM.
    """

#Read options from terminal, list of options available with graph -h
    args = run_argparse()

#Define a cluster class for QM, MM and QMMM regions
    QM_C = Cluster( "QM" )
    MM_C = Cluster( "MM")
    QMMM_C = Cluster( "QMMM")

#Search pattern that matches a coordinate entry that has an arbitrary label for atom
    pat = re.compile(r'^\s*\w+\S*\s+-*\d+\.\d*\s+-*\d+\.\d*\s+-*\d+\.\d*')

#Store all atoms from QM region
    for i in open(args.qm).readlines():
        if pat.match(i):
            QM_C.AddAtom( Atom( i, args, qm= True ))

#Store all atoms in mm region that do not have identical coordinates to any atom in the QM region
    for i in open(args.mm).readlines():
        if pat.match(i):
#Add all atoms to MM region by default, unless if its water and water cutoff
            if args.waters:
                if i.split()[0].split('-')[1] == 'T3':
                    if QM_C.Atomlist[7].Dist( Atom( i, args) ) > args.waters:
                        continue
            mm_atom = Atom( i, args)

            MM_C.AddAtom( mm_atom )
            for j in QM_C:
                r = j.Dist( mm_atom )
                if  r < 0.1:
#Remove MM atoms without invoking property transfer methods for obvious QM atoms
                    MM_C.Atomlist.remove( mm_atom )
                if  r < args.region_cutoff:
#Atoms in a region within bond capping region placed in QM/MM region
                    QMMM_C.AddAtom( mm_atom )


    print "Done storing MM and QMMM clusters"
    print "QM atoms: %d\nMM atoms: %d\nQMMM atoms: %d\n" %(len(QM_C), len(MM_C), len(QMMM_C))
    print "---------------------------\nUpdating close atoms in QMMM:"

    QMMM_C.UpdateCloseAtoms( MM_C )
    print "Done updating"
    print "Information on QMMM cluster:"
    for i in QMMM_C:
        if i.label == "62C-X1-C":
            print i

#Now time for elimination of atoms which are poisonous.
#Calculate distance between each atom in the boarder region and remove it
#if it is too cloose to qm as defined by cutoff
    #print_f.pre_elimination(qm_region, qm_mm_region, args.cutoff)


    print "---------------------------\nBegin elimination\n\n"
    for i in QMMM_C:
        for j in QM_C:
            r = i.Dist(j)
            if r < float(args.cutoff):
#If this atom still exists in mm region, it should be removed, and its properties transfered defiend by args implemented in RemoveAtom method
                MM_C.RemoveAtom( i, args )
                continue

    print "-----------------------------\nPost elimination:"
    print "Lenght of MM: %d\nlength of QMMM: %d" %(len(MM_C), len(QMMM_C))


    f = open(args.outfile + "_" + str(args.cutoff) + '.pro', 'w')
    f.write('AA\n')
    f.write('%d %d %d %d\n'%tuple([len(MM_C),int(args.angular),int(args.charge),len(args.polar) ] ))
    for i in MM_C:
        f.write("{0:10s}  {1:10f}  {2:10f}  {3:10f}  {4:10f}  {5:10f}\n".format( i.label, i.x, i.y, i.z, i.Props.Charge, i.Props.Polar))
    f.close()

    pat = re.compile(r'(^\d+)')
    f = open(args.outfile + "_" + str(args.cutoff) + '.pot', 'w')
    f.write('AA\n')
    f.write('%d %d %d %d\n'%tuple([len(MM_C),int(args.angular),int(args.charge),len(args.polar) ] ))
    for i in MM_C:
        f.write("{0:10s}  {1:10f}  {2:10f}  {3:10f}  {4:10f}  {5:10f}\n".format( pat.match( i.label).group(1) , i.x, i.y, i.z, i.Props.Charge, i.Props.Polar))
    f.close()


    if args.xyz_files:
        f = open( args.outfile + "_QM_MM.xyz", 'w')
        f.write( str(len( QM_C) + len(MM_C)) + "\n\n" )
        for i in QM_C:
            f.write("{0:15s}{1:10f}{2:10f}{3:10f}\n".format( i.atype, i.x, i.y, i.z))
        for i in MM_C:
            f.write("{0:15s}{1:10f}{2:10f}{3:10f}\n".format( i.atype, i.x, i.y, i.z))
        f.close()

        f = open( args.outfile + "_QM_QMMM.xyz", 'w')
        f.write( str(len( QMMM_C)+ len(QM_C)) + "\n\n" )

        for i in QMMM_C:
            f.write("{0:15s}{1:10f}{2:10f}{3:10f}\n".format( i.atype, i.x, i.y, i.z))
        for i in QM_C:
            f.write("{0:15s}{1:10f}{2:10f}{3:10f}\n".format( i.atype, i.x, i.y, i.z))
        f.close()

if __name__ == '__main__':
    main()
