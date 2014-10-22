#!/usr/bin/env python
#-*- coding: utf-8 -*-

import argparse, re, fractions

import numpy as np
import math as m

from molecules import Water, Atom, Property, Methanol

a0 = 5.2917721092e-11

class Generator( dict ):
    """
    class to create molecules, write dalton .mol files 
    using -params for study with use_calculator.py

    water currently implemented

    """
    def __init__(self, *args, **kwargs):

        self[ ("water","a_hoh", "degree") ] = 104.5
        self[ ("water","r_oh", "AA") ] = 1.83681

    @staticmethod
    def build_molecule( molecule ):
        """ molecule is a string to build, return class """

    def gen_mols_params(self, mol, AA = False,):

        r = np.r_[ self.optionsR[ "min" ] : self.optionsR[ "max" ] : \
                complex( "%sj"%self.optionsR[ "points" ] ) ]

        tau = np.r_[ self.optionsTau[ "min" ] : self.optionsTau[ "max" ] : \
                complex( "%sj"%self.optionsTau[ "points" ] ) ]

        theta = np.r_[ self.optionsTheta[ "min" ] : self.optionsTheta[ "max" ] : \
                complex( "%sj"%self.optionsTheta[ "points" ] ) ]

        rho1 = np.r_[ self.optionsRho1[ "min" ] : self.optionsRho1[ "max" ] : \
                complex( "%sj"%self.optionsRho1[ "points" ] ) ]

        rho2 = np.r_[ self.optionsRho2[ "min" ] : self.optionsRho2[ "max" ] : \
                complex( "%sj"%self.optionsRho2[ "points" ] ) ]

        rho3 = np.r_[ self.optionsRho3[ "min" ] : self.optionsRho3[ "max" ] : \
                complex( "%sj"%self.optionsRho3[ "points" ] ) ]

        r_oh = self[ ("water","r_oh", "AA") ]
        a_hoh = np.pi * self[ ("water", "a_hoh", "degree" )] / 180.0

        for i in r:
            for j in tau:
                for k in theta:
                    for l in rho1:
                        for m in rho2:
                            for n in rho3:
                                w1 = self.get_water_params( [0, 0, 0], r_oh, a_hoh, AA = AA)
                                x, y, z = self.polar_to_cartesian( i, j, k )
                                w2 = self.get_water_params( [x,y,z], r_oh, a_hoh)
                                w2.rotate( l, m, n )
                                name = ""
                                name += "-".join( map( str, ["%3.2f"%i, "%3.2f"%j, "%3.2f"%k, "%3.2f"%l, "%3.2f"%m, "%3.2f"%n] ) )
                                name += ".mol"
                                self.write_mol( [w1, w2], name )

    def vary_parameters( self, *args ):
        """Given two parameters, for e.g. r and theta, keeps all other static"""
        if args:
            for j in args:
                for i in j:
                    if i == "r":
                        self.varyR = True
                        self.optionsR = j[i]
                    if i == "tau":
                        self.varyTau = True
                        self.optionsTau = j[i]
                    if i == "theta":
                        self.varyTheta = True
                        self.optionsTheta = j[i]
                    if i == "rho1":
                        self.varyRho1 = True
                        self.optionsRho1 = j[i]
                    if i == "rho2":
                        self.varyRho2 = True
                        self.optionsRho2 = j[i]
                    if i == "rho3":
                        self.varyRho3 = True
                        self.optionsRho3 = j[i]

    def get_water_params(self, origin, r, theta, AA = True ):
        h1 = Atom() ; h2 = Atom() ; o = Atom()
        d = (m.pi/2 - theta/2)
        o.element = "O" ; h1.element = "H" ; h2.element = "H"
        o.x = origin[0] ; o.y = origin[1] ; o.z = origin[2] 
        h1.x = (origin[0] + r * m.cos(d)) ; h1.y = origin[1] ; h1.z = (origin[2] + r*m.sin(d))
        h2.x = (origin[0] - r * m.cos(d)) ; h2.y = origin[1] ; h2.z = (origin[2] + r*m.sin(d))
        w = Water(); w.append( o) ;w.append( h2 ) ;w.append( h1 ) 
        w.a_hoh = theta
        w.r_oh = r
        w.center = origin
        w.euler1 = 0.00
        w.euler2 = 0.00
        w.euler3 = 0.00

        if AA:
            w.AA = True
            w.h1.AA = True
            w.h2.AA = True
            w.o.AA  = True
        else:
            w.AA = False
            w.h1.AA = False
            w.h2.AA = False
            w.o.AA  = False
        return w

    def get_mol( self, center, mol, AA = True ):
        """return molecule in origo, all molecules have different definition
        of euler

        for water place O in origo
        for methanol place C=O bond in origo
        
        """

        if mol == "water":
            h1 = Atom( { "AA" : AA} )
            h2 = Atom( { "AA" : AA} )
            o =  Atom( { "AA" : AA} )

            r_oh = self[ ("water","r_oh", "AA") ]
            if not AA:
                r_oh = r_oh / a0

            a_hoh = self[ ("water","a_hoh","degree") ]

            d = (m.pi/2 - a_hoh/2)
            origin = np.array( [ 0, 0, 0] )
            o.element = "O" ; h1.element = "H" ; h2.element = "H"
            o.x = origin[0]
            o.y = origin[1]
            o.z = origin[2] 

            h1.x = (origin[0] + r_oh * m.cos(d))
            h1.y = origin[1] 
            h1.z = (origin[2] + r_oh*m.sin(d))

            h2.x = (origin[0] - r_oh * m.cos(d)) 
            h2.y = origin[1] 
            h2.z = (origin[2] + r_oh*m.sin(d))

            w = Water()
            w.append( o )
            w.append( h1 )
            w.append( h2 )
            if AA:
                w.AA = True
                w.h1.AA = True
                w.h2.AA = True
                w.o.AA  = True
            else:
                w.AA = False
                w.h1.AA = False
                w.h2.AA = False
                w.o.AA  = False
            return w
        elif mol == "methanol":
            r_co = self[ ("methanol", "r_co", "AA" )]["r_co"]
            r_oh = self[ ("methanol", "r_oh", "AA" )]["r_oh"]
            r_ch = self[ ("methanol", "r_ch", "AA" )]["r_ch"]
            a_coh = self[ ("methanol", "a_coh", "AA" )]["a_coh"]
            a_hco = self[ ("methanol", "a_hco", "AA" )]["a_hco"]
            a_hch = self[ ("methanol", "a_hch", "AA" )]["a_hch"]

            c = Atom( {"x":0,"y":-r_co/2,"z":0,"AA" : AA , "element":"C"} )
            o = Atom( {"x":0,"y":r_co/2,"z":0,"AA" : AA , "element":"C"} )

    def readWaters(self, fname):
        """From file with name fname, return a list of all waters encountered"""
#If the file is plain xyz file

        atoms = []
        if fname.endswith( ".xyz" ) or fname.endswith(".mol"):
            pat_xyz = re.compile(r'^\s*(\w+)\s+(-*\d*.+\d+)\s+(-*\d*.+\d+)\s+(-*\d*.+\d+) *$')
            for i in open( fname ).readlines():
                if pat_xyz.match(i):
                    f = pat_xyz.match(i).groups()
                    tmpAtom = Atom()
                    tmpAtom.AA = True
                    tmpAtom.x = float(f[1])
                    tmpAtom.y = float(f[2])
                    tmpAtom.z = float(f[3])
                    tmpAtom.element = f[0][0]
                    atoms.append( tmpAtom )

        elif fname.endswith( ".pdb" ):
            pat1 = re.compile(r'^(ATOM|HETATM)')
            for i in open( fname ).readlines():
                if pat1.search(i):
                    #Ignore charge centers for polarizable water models
                    if ( i[11:16].strip() == "SW") or (i[11:16] == "DW"):
                        continue
                    tmpAtom = Atom(i[11:16].strip()[0], \
                            float(i[30:38].strip()), \
                            float(i[38:46].strip()), \
                            float(i[46:54].strip()), \
                            int(i[22:26].strip()) )

                    if fnameAAorAU == "AU":
                        if args.opAAorAU == "AA":
                            tmpAtom.toAA()
                    elif fnameAAorAU == "AA":
                        if args.opAAorAU == "AU":
                            tmpAtom.toAU()
                    atoms.append( tmpAtom )
        elif fname.endswith( ".out" ):
            pat_xyz = re.compile(r'^(\w+)\s+(-*\d*.+\d+)\s+(-*\d*.+\d+)\s+(-*\d*.+\d+) *$')
            for i in open( fname ).readlines():
                if pat_xyz.match(i):
                    f = pat_xyz.match(i).groups()
                    tmpAtom = Atom(f[0][0], float(f[1]), float(f[2]), float(f[3]), 0)
                    if fnameAAorAU == "AU":
                        if args.opAAorAU == "AA":
                            tmpAtom.toAA()
                    elif fnameAAorAU == "AA":
                        if args.opAAorAU == "AU":
                            tmpAtom.toAU()
                    atoms.append( tmpAtom )
#loop over oxygen and hydrogen and if they are closer than 1 A add them to a water
        waters = []
        cnt = 1

        if fname.endswith( ".xyz" ) or fname.endswith(".mol"):
            for i in atoms:
                if i.element == "H":
                    continue
                if i.inWater:
                    continue
                tmp = Water(  )
                i.inWater = True
                tmp.append( i )
                for j in atoms:
                    if j.element == "O":
                        continue
                    if j.inWater:
                        continue
#If in cartesian:
                    if j.AA:
                        if i.distToAtom(j) < 1.1:
                            tmp.append ( j )
                            j.inWater = True
                    else:
                        if i.distToAtom(j) < 1.1/a0:
                            tmp.append ( j )
                            j.inWater = True
                tmp.number = cnt
                cnt += 1
                waters.append( tmp )
        elif fname.endswith( ".pdb" ):
#Find out the size of the box encompassing all atoms
            xmin = 10000.0; ymin = 10000.0; zmin = 10000.0; 
            xmax = -10000.0; ymax = -10000.0; zmax = -10000.0; 
            for i in atoms:
                if i.x < xmin:
                    xmin = i.x
                if i.y < ymin:
                    ymin = i.y
                if i.z < zmin:
                    zmin = i.z
                if i.x > xmax:
                    xmax = i.x
                if i.y > ymax:
                    ymax = i.y
                if i.z > zmax:
                    zmax = i.z
            center = np.array([ xmax - xmin, ymax -ymin, zmax- zmin]) /2.0
            wlist = []
            for i in atoms:
                if i.element != "O":
                    continue
                tmp = Water()
                i.inWater= True
#__Water__.append() method will update the waters residue number and center coordinate
#When all atoms are there
#Right now NOT center-of-mass
                tmp.append(i)
                for j in atoms:
                    if j.element != "H":
                        continue
                    if j.inWater:
                        continue
#1.05 because sometimes spc water lengths can be over 1.01
                        
                    if args.opAAorAU == "AA":
                        if i.dist(j) <= 1.05:
                            j.inWater = True
                            tmp.append( j )
                            if len(tmp.atomlist) == 3:
                                break
                    elif args.opAAorAU == "AU":
                        if i.dist(j) <= 1.05/a0:
                            j.inWater = True
                            tmp.append( j )
                            if len(tmp.atomlist) == 3:
                                break
                wlist.append( tmp )
            wlist.sort( key = lambda x: x.distToPoint( center ))
            center_water = wlist[0]
            cent_wlist = wlist[1:]
            cent_wlist.sort( key= lambda x: x.distToWater( center_water) )
            waters = [center_water] + cent_wlist[ 0:args.waters - 1 ]
        elif fname.endswith( ".out" ):
            for i in atoms:
                if i.element == "H":
                    continue
                if i.inWater:
                    continue
                tmp = Water(  )
                i.inWater = True
                tmp.append( i )
                for j in atoms:
                    if j.element == "O":
                        continue
                    if j.inWater:
                        continue
#If in cartesian:
                    if i.AA:
                        if i.dist(j) < 1.0:
                            tmp.append ( j )
                            j.inWater = True
                    else:
                        if i.dist(j) < 1.0/a0:
                            tmp.append ( j )
                            j.inWater = True
                tmp.number = cnt
                cnt += 1
                waters.append( tmp )
        return waters

    def write_mol(self, wlist, name = "tmp.mol" ):
        f_ = open (name, 'w')
        f_.write("ATOMBASIS\n\n\nAtomtypes=2 Charge=0 Nosymm\n")
        f_.write("Charge=1.0 Atoms=4 Basis=cc-pVTZ\n")
        for i in wlist:
            for j in i:
                if j.element != "H":
                    continue
                f_.write( "%s %.5f %.5f %.5f\n" %(j.element, j.x, j.y, j.z ) )
        f_.write("Charge=8.0 Atoms=2 Basis=cc-pVTZ\n")
        for i in wlist:
            for j in i:
                if j.element != "O":
                    continue
                f_.write( "%s %.5f %.5f %.5f\n" %(j.element, j.x, j.y, j.z ) )
        f_.close()

    def polar_to_cartesian(self, r, tau, theta):
        x, y, z = r* m.sin( theta )*m.cos( tau ) \
               , r* m.sin(  theta )*m.sin( tau )  \
               , r* m.cos(  theta ) 

        return x , y , z

if __name__ == '__main__':

    A = argparse.ArgumentParser( add_help= True)
#Related to generating two water molecules with specified 6 parameters
    A.add_argument( "-params", type = str, 
            default = "water", choices = ["water",
                "methanol" , "ethane"] )

    A.add_argument( "-r"     ,   type = float , default = 3.00  ) 
    A.add_argument( "-theta" ,   type = float , default = 0.00  ) 
    A.add_argument( "-tau"   ,   type = float , default = 0.00  ) 
    A.add_argument( "-rho1"  ,   type = float , default = 0.00  ) 
    A.add_argument( "-rho2"  ,   type = float , default = 0.00  ) 
    A.add_argument( "-rho3"  ,   type = float , default = 0.00  ) 

    A.add_argument( "-r_max"     ,   type = float , default = 10.0 ) 
    A.add_argument( "-theta_max" ,   type = float , default = np.pi/2    ) 
    A.add_argument( "-tau_max"   ,   type = float , default = np.pi/2    ) 
    A.add_argument( "-rho1_max"  ,   type = float , default = np.pi/2    ) 
    A.add_argument( "-rho2_max"  ,   type = float , default = np.pi      ) 
    A.add_argument( "-rho3_max"  ,   type = float , default = np.pi/2    ) 

    A.add_argument( "-r_points"     ,   type = int , default = 1) 
    A.add_argument( "-theta_points" ,   type = int , default = 1) 
    A.add_argument( "-tau_points"   ,   type = int , default = 1) 
    A.add_argument( "-rho1_points"  ,   type = int , default = 1) 
    A.add_argument( "-rho2_points"  ,   type = int , default = 1) 
    A.add_argument( "-rho3_points"  ,   type = int , default = 1) 

#Related to generating/getting one water
    A.add_argument( "-get_mol", 
            type = str, default = "water", choices = ["water",
                "methanol" , "ethane"] )
    A.add_argument( "-get_r_oh" , type = float,  default = 1.83681 )
    A.add_argument( "-get_a_hoh", type = float,  default = 104.5 )
    A.add_argument( "-get_rho1" , type = float,  default = 0.0 )
    A.add_argument( "-get_rho2" , type = float,  default = 0.0 )
    A.add_argument( "-get_rho3" , type = float,  default = 0.0 )
    A.add_argument( "-AA" ,  default = False, action = 'store_true' )
    args = A.parse_args()

    g = Generator()

    opts =  {
       "r"    :{"min":args.r,    "max": args.r_max,    "points":args.r_points,    },
       "tau"  :{"min":args.tau,  "max": args.tau_max,  "points":args.tau_points,  },
       "theta":{"min":args.theta,"max": args.theta_max,"points":args.theta_points,},
       "rho1" :{"min":args.rho1, "max": args.rho1_max, "points":args.rho1_points, },
       "rho2" :{"min":args.rho2, "max": args.rho2_max, "points":args.rho2_points, },
       "rho3" :{"min":args.rho3, "max": args.rho3_max, "points":args.rho3_points, },
             }


    if args.params:
        g.vary_parameters( opts )
        g.gen_mols_params( args.params , AA = False )

    if args.get_mol:
        mol = g.get_mol( center = [0, 0, 0 ], mol = args.get_mol , AA = False )
