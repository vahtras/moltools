#!/usr/bin/env python
#-*- coding: utf-8 -*-

import re, os, math, operator, numpy
from sys import argv
import argparse

def uniq(input):
    output = []
    rest = []
    for x in input:
        if x not in output:
            output.append(x)
        else:
            rest.append(x)
    return output# , rest
def main():
    A = argparse.ArgumentParser( add_help= True)
    A.add_argument( "-b", default = "Turbomole-TZVP" ,help = "Name of basis")
    A.add_argument( "-c", type = int, default = 0)
    A.add_argument( "-f", type = str )
    A.add_argument( "--ano", type = str , help='toggle ano basis, supply size: [small, middle, heavy, large ]')
    A.add_argument( "--anobasis", nargs = '*' , help='set ano size to use, default [H: 2s, C/N/O: 3s2p1d, S:4s3p2d1f]')
    A.add_argument( "--center", type = int, default = 5 )
    A.add_argument( "--cmethod", type = str, default ='Centerofmass' )
    A.add_argument( "--method", type = str , default = 'B3LYP' )

    a = A.parse_args()

    if not a.f:
        print "Supply pdb file using -f <file>"
        raise SystemExit
    name = a.f.split('.')[0]

    basis = a.b
    charge = a.c

# Dictionaries, small is 6-31G, middle is 6-31+G*, heavy is ?? large is 6-311++G**
    adict = { '1' : 'H', '6': 'C', '7': 'N', '8': 'O', '16':'S' }
    small = { "H" : "ano-1 1", "C" : "ano-1 3 2" ,
            "N" : "ano-1 3 2", "O" : "ano-1 3 2",
            "S" : "ano-2 4 3" }
    middle = { "H" : "ano-1 2", "C" : "ano-1 4 3 1" ,
            "N" : "ano-1 4 3 1", "O" : "ano-1 4 3 1",
            "S" : "ano-2 5 4 1" }
    heavy = { "H" : "ano-1 3 1", "C" : "ano-1 5 3 2" ,
            "N" : "ano-1 5 3 2", "O" : "ano-1 5 3 2" }
    large = { "H" : "ano-1 4 1", "C" : "ano-1 5 4 1" ,
            "N" : "ano-1 5 4 1", "O" : "ano-1 5 4 1",
            "S" : "ano-2 7 6 1" }

    ano_dict = { "small": small, "middle": middle, "heavy": heavy, "large": large}

    ftext = open(a.f).readlines()

#Get atoms
    pat1 = re.compile(r'^(ATOM|HETATM)')
    alist = []
    hlist = []
    olist = []
    for i in ftext:
        if pat1.search(i):
            if ( i[11:16].strip() == "SW") or (i[11:16] == "DW"):
                continue
            if i[11:16].strip()[0] == "H":
                hlist.append( PdbEntry(i) )
            if i[11:16].strip()[0] == "O":
                olist.append( PdbEntry(i) )
            alist.append( PdbEntry(i) )
#
    xmin = 10000.0; ymin = 10000.0; zmin = 10000.0; 
    xmax = -10000.0; ymax = -10000.0; zmax = -10000.0; 
    for i in alist:
        if float(i.x) < xmin:
            xmin = float(i.x)
        if float(i.y) < ymin:
            ymin = float(i.y)
        if float(i.z) < zmin:
            zmin = float(i.z)
        if float(i.x) > xmax:
            xmax = float(i.x)
        if float(i.y) > ymax:
            ymax = float(i.y)
        if float(i.z) > zmax:
            zmax = float(i.z)
    center = [ xmax - xmin, ymax -ymin, zmax- zmin]
    wlist = []
    for i in olist:
        tmp = Water()
        i.inWater= True
        tmp.AddAtom(i)
        for j in hlist:
            if j.inWater:
                continue
            if i.Dist(j) < 1.0:
                tmp.AddAtom( j )
        tmp.SetCenter()
        tmp.SetCenterOfMass()
        tmp.GetDistFromCenter( center )
        wlist.append( tmp )
    wlist.sort( key=operator.attrgetter('distFromCenter'))


    wlist = wlist[0:a.center]
    alist = []

    ave = numpy.array([0.0 for i in range(len(wlist) -1) ])
    tmin =  1000.0
    tmax = -1000.0
    for i in wlist:
        tave = []
        for j in wlist:
            if i == j:
                continue
            if a.cmethod == "Regular":
                tave.append( i.Dist(j) )
            elif a.cmethod == "Centerofmass":
                tave.append( i.DistCenterOfMass(j) )
        srt = sorted(tave)
        if srt[0] < tmin:
            tmin = srt[0]
        if srt[-1] > tmax:
            tmax = srt[-1]
        ave += numpy.array( srt )
        for j in i.atomlist:
            alist.append( j )
    ave /= len(wlist)

    name = a.f.rstrip('.pdb')

#Write mol file from remaining atoms

    h = []; c = []; n = []; o = []; s = []; atoms=[]
    pat = re.compile(r'^ *([A-Z])')
    for j in alist:
        if j.Name == 'H': h.append( j ) ; atoms.append('H')
        if j.Name == 'C': c.append( j ) ; atoms.append('C')
        if j.Name == 'N': n.append( j ) ; atoms.append('N')
        if j.Name == 'O': o.append( j ) ; atoms.append('O')
        if j.Name == 'S': s.append( j ) ; atoms.append('S')
    atomtypes = len(uniq(atoms))  

    name = a.f.rstrip('.pdb')
    comment = ""

    if a.ano:
        name += '_%s' %a.ano
    if a.method == 'B3LYP':
        name += '_%s' %a.method
    try:
        comment += "Average : %.2f %.2f %.2f Method: %s\n"\
                %(ave[0], ave[1], ave[2], a.cmethod)
        comment +=  "Min: %.2f Max: %.2f\n" %(tmin, tmax)
        name += '_%d' %a.center
    except IndexError:
        comment += "Average : %.2f Method: %s\n"\
                %(ave[0], a.cmethod)
        comment +=  "Min: %.2f Max: %.2f\n" %(tmin, tmax)
        name += '_%d' %a.center

    name += '.mol'
    if a.ano:
        if a.ano not in ano_dict:
            print "Supply correct ano type please, [small, middle, heavy, large]"
            raise SystemExit
        filen = open( name , 'w' )
        filen.write('ATOMBASIS\n')
        filen.write(comment)
        filen.write('Atomtypes=%s Charge=%s Angstrom Nosymm\n'%(str(atomtypes),str(charge)))
        if h:
            filen.write('Charge=1.0 Atoms=%d Basis=%s\n'%(len(h),ano_dict[a.ano]['H'] ))
            for k in h:
                filen.write('{0:15s}{1:15s}{2:15s}{3:15s}\n'.format(
                    *tuple([ k.Name, k.x, k.y, k.z ])))
        if c:
            filen.write('Charge=6.0 Atoms=%d Basis=%s\n'%(len(c),ano_dict[a.ano]['C']))
            for k in c:
                filen.write('{0:15s}{1:15s}{2:15s}{3:15s}\n'.format(
                    *tuple([ k.Name, k.x, k.y, k.z ])))
        if n:
            filen.write('Charge=7.0 Atoms=%d Basis=%s\n'%(len(n),ano_dict[a.ano]['N']))
            for k in n:
                filen.write('{0:15s}{1:15s}{2:15s}{3:15s}\n'.format(
                    *tuple([ k.Name, k.x, k.y, k.z ])))
        if o:
            filen.write('Charge=8.0 Atoms=%d Basis=%s\n'%(len(o),ano_dict[a.ano]['O']))
            for k in o:
                filen.write('{0:15s}{1:15s}{2:15s}{3:15s}\n'.format(
                    *tuple([ k.Name, k.x, k.y, k.z ])))
        if s:
            filen.write('Charge=16.0 Atoms=%d Basis=%s\n'%(len(s),ano_dict[a.ano]['S']))
            for k in s:
                filen.write('{0:15s}{1:15s}{2:15s}{3:15s}\n'.format(
                    *tuple([ k.Name, k.x, k.y, k.z ])))
        filen.close()
    else:
        filen = open( name , 'w' )
        filen.write('BASIS\n')
        filen.write('%s\n' % basis)
        filen.write(comment)
        filen.write('Atomtypes=%s Charge=%s Angstrom Nosymm\n'%(str(atomtypes),str(charge)))
        if h:
            filen.write('Charge=1.0 Atoms=%d \n'%len(h))
            for k in h:
                filen.write('{0:15s}{1:15s}{2:15s}{3:15s}\n'.format(
                    *tuple([ k.Name, k.x, k.y, k.z ])))
        if c:
            filen.write('Charge=6.0 Atoms=%d \n'%len(c))
            for k in c:
                filen.write('{0:15s}{1:15s}{2:15s}{3:15s}\n'.format(
                    *tuple([ k.Name, k.x, k.y, k.z ])))
        if n:
            filen.write('Charge=7.0 Atoms=%d \n'%len(n))
            for k in n:
                filen.write('{0:15s}{1:15s}{2:15s}{3:15s}\n'.format(
                    *tuple([ k.Name, k.x, k.y, k.z ])))
        if o:
            filen.write('Charge=8.0 Atoms=%d \n'%len(o))
            for k in o:
                filen.write('{0:15s}{1:15s}{2:15s}{3:15s}\n'.format(
                    *tuple([ k.Name, k.x, k.y, k.z ])))
        if s:
            filen.write('Charge=16.0 Atoms=%d \n'%len(s))
            for k in s:
                filen.write('{0:15s}{1:15s}{2:15s}{3:15s}\n'.format(
                    *tuple([ k.Name, k.x, k.y, k.z ])))
        filen.close()

class Atom:
    def __init__(self, line):
        self.Name = line.split()[0].strip()
        self.x = line.split()[1].strip()
        self.y = line.split()[2].strip()
        self.z = line.split()[3].strip()
        self.inWater = False

class Water:
    def __init__(self):
        self.atomlist = []
    def AddAtom(self, other):
        self.atomlist.append( other )
    def SetCenter(self):
        assert ( len(self.atomlist) == 3)
        x = 0; y= 0; z = 0;
        for i in self.atomlist:
            x+= float(i.x); y+= float(i.y); z+= float(i.z); 
        self.center = [ x/3, y/3, z/3]
    def SetCenterOfMass(self):
        assert ( len(self.atomlist) == 3)
        x = 0; y= 0; z = 0;
        for i in self.atomlist:
            if i.Name == 'O':
                x+= 16.0*float(i.x); y+= 16.0*float(i.y); z+= 16.0*float(i.z); 
            else:
                x+= float(i.x); y+= float(i.y); z+= float(i.z); 
        self.centerOfMass = numpy.array([ x, y, z]) / 18.0

    def GetDistFromCenter( self, coord):
        xyz = self.center
        self.distFromCenter = math.sqrt( (xyz[0] -coord[0])**2 +\
                (xyz[1] -coord[1])**2  + (xyz[2] -coord[2])**2 )
    def Dist(self, other):
        xyz1 = self.center
        xyz2 = other.center
        return math.sqrt( (xyz1[0] - xyz2[0])**2 + \
            (xyz1[1] - xyz2[1])**2 + (xyz1[2] - xyz2[2])**2 )
    def DistCenterOfMass(self, other):
        xyz1 = self.centerOfMass
        xyz2 = other.centerOfMass
        return math.sqrt( (xyz1[0] - xyz2[0])**2 + \
            (xyz1[1] - xyz2[1])**2 + (xyz1[2] - xyz2[2])**2 )


class PdbEntry:
    def __init__(self, line ):
        self.PDBName = line[11:16].strip()
        self.Name = line[11:16].strip()[0]
        self.ResName = line[17:20]
        self.ResNum =line[22:26].strip()
        self.x = line[30:38].strip()
        self.y = line[38:46].strip()
        self.z = line[46:54].strip()
        self.inWater = False
    def __str__(self):
        return (' ' + self.Name + '(PDBName=' + self.PDBName 
                + ',ResName=' + self.ResName + ',ResNum=' 
                + self.ResNum + ')' + '\t' + self.x +
                '  ' + self.y + '  '+ self.z +"\n" )
    def Dist(self, other):
        return math.sqrt((float(self.x)-float(other.x))**2 + \
            (float(self.y)-float(other.y))**2 + (float(self.z)-float(other.z))**2)

if __name__ == '__main__':
    main()
