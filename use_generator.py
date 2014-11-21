#!/usr/bin/env python
#-*- coding: utf-8 -*-

import argparse, re, fractions

import numpy as np
import math as m

from template import Template
from molecules import Molecule, Water, Atom, Property, Cluster

a0 = 0.52917721092

class Generator( dict ):
    """
    class to create molecules, write dalton .mol files 
    using -params for study with use_calculator.py

    water currently implemented

    """
    def __init__(self, *args, **kwargs):

#This waater is TIP3P model, generalzie later
        self[ ("water","a_hoh", "degree") ] = 104.5
        self[ ("water","r_oh", "AA") ] = 0.9572

        self[ ("methanol", "r_oh", "AA" ) ] = 0.967
        self[ ("methanol", "r_co", "AA" ) ] = 1.428
        self[ ("methanol", "r_ch", "AA" ) ] = 1.098

        self[ ("methanol", "a_coh", "degree" ) ] = 107.16
        self[ ("methanol", "a_hch", "degree" ) ] = 109.6
        self[ ("methanol", "a_hco", "degree" ) ] = 109.342

        self[ ("methanol", "d_hcoh", "h4", "degree" ) ] =  60.0
        self[ ("methanol", "d_hcoh", "h5", "degree" ) ] = -60.0
        self[ ("methanol", "d_hcoh", "h6", "degree" ) ] =  180.0

        if kwargs is not {}:
            self.options = kwargs.get( 'r', {"max": 10 , "min" : 3, "points":1} )

    @staticmethod
    def build_molecule( molecule ):
        """ molecule is a string to build, return class """

    @staticmethod
    def get_hfqua_dal():
        return """**DALTON INPUT
.RUN RESPONSE
.DIRECT
.PARALLELL
**WAVE FUNCTION
.HF
**RESPONSE
*QUADRATIC
.DIPLEN
**END OF DALTON INPUT"""

    def gen_mols_param(self, mol = "water", 
            basis = ["ano-1 2 1", "ano-1 3 2 1"],
            AA = True):

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
                                c= Cluster()
                                w1 = self.get_mol( [0, 0, 0], mol , AA = AA)
                                c.append( w1, in_qm = True )
                                x, y, z = self.polar_to_cartesian( i, j, k )
                                w2 = self.get_mol( [x,y,z], mol, AA = AA)
                                w2.rotate( l, m, n )

                                c.append( w2, in_qm = True )
                                name = ""
                                name += "-".join( map( str, ["%3.2f"%i, "%3.2f"%j, "%3.2f"%k, "%3.2f"%l, "%3.2f"%m, "%3.2f"%n] ) )
                                name += ".mol"

                                m = c.get_qm_mol_string( AA = AA,
                                        basis = tuple(basis),
                                        )
                                f_ = open(name, 'w')
                                f_.write( m)

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

    def get_mol( self, center, mol, AA = True ):
        """return molecule in center, all molecules have different definition
        of euler angles

        for water place O in origo
        for methanol place C=O bond in origo
        
        """

        if mol == "water":
#Geometrical parameters
            r_oh = self[ ("water","r_oh", "AA") ]
            a_hoh = self[ ("water","a_hoh","degree") ]
            if not AA:
                r_oh = r_oh / a0

            d = (90 - a_hoh/2 ) * np.pi / 180


            xo = center[0]
            yo = center[1]
            zo = center[2] 

            xh1 = (center[0] + r_oh * np.cos(d))
            yh1 =  center[1] 
            zh1 = (center[2] + r_oh* np.sin(d))

            xh2 = (center[0] - r_oh * np.cos(d)) 
            yh2 = center[1] 
            zh2 = (center[2] + r_oh* np.sin(d))

            h1 = Atom( **{ "AA" : AA,
                "x" : xh1,
                "y" : yh1,
                "z" : zh1,
                "element" : "H"} )
            h2 = Atom( **{ "AA" : AA,
                "x" : xh2,
                "y" : yh2,
                "z" : zh2,
                "element" : "H"} )
            o = Atom( **{ "AA" : AA,
                "x" : xo,
                "y" : yo,
                "z" : zo,
                "element" : "O"} )

            w = Water()
            w.append( o )
            w.append( h1 )
            w.append( h2 )
            
            return w

        elif mol == "methanol":

            r_co = self[ ("methanol", "r_co", "AA" )]
            r_oh = self[ ("methanol", "r_oh", "AA" )]
            r_ch = self[ ("methanol", "r_ch", "AA" )]

            a_coh = self[ ("methanol", "a_coh", "degree" ) ]
            #a_hch = self[ ("methanol", "a_hch", "degree" ) ]
            a_hco = self[ ("methanol", "a_hco", "degree" ) ]

            a_coh *= np.pi / 180
            a_hco *= np.pi / 180

            d_hcoh_4 = self[ ("methanol", "d_hcoh", "h4", "degree" ) ]
            d_hcoh_4 *= np.pi / 180
            d_hcoh_5 = self[ ("methanol", "d_hcoh", "h5", "degree" ) ]
            d_hcoh_5 *= np.pi / 180
            d_hcoh_6 = self[ ("methanol", "d_hcoh", "h6", "degree" ) ]
            d_hcoh_6 *= np.pi / 180

            if not AA:
                r_co, r_oh, r_ch = r_co/a0, r_oh/a0, r_ch/a0

            c1 = Atom( **{"x":0, "y":0, "z":-r_co/2, "AA": AA, "element":"C" } )
            o2 = Atom( **{"x":0, "y":0, "z": r_co/2, "AA": AA, "element":"O" } )

            h3 = Atom( **{"x":r_oh*np.cos( a_coh-np.pi/2),
                "y":0,
                "z":r_oh*np.sin( a_coh-np.pi/2) + r_co/2,
                "AA": AA, "element":"H" } )

            h4 = Atom( **{"x": r_ch*np.sin( a_hco ) * np.cos( d_hcoh_4 ),
                "y": r_ch*np.sin( a_hco) * np.sin( d_hcoh_4 ),
                "z": r_ch*np.cos( a_hco) - r_co/2 ,
                "AA": AA, "element":"H" } )
            h5 = Atom( **{"x": r_ch*np.sin( a_hco ) * np.cos( d_hcoh_5 ),
                "y": r_ch*np.sin( a_hco) * np.sin( d_hcoh_5 ),
                "z": r_ch*np.cos( a_hco) - r_co/2 ,
                "AA": AA, "element":"H" } )
            h6 = Atom( **{"x": r_ch*np.sin( a_hco ) * np.cos( d_hcoh_6 ),
                "y": r_ch*np.sin( a_hco) * np.sin( d_hcoh_6 ),
                "z": r_ch*np.cos( a_hco) - r_co/2 ,
                "AA": AA, "element":"H" } )

            m = Methanol()
            m.append(c1)
            m.append(o2)
            m.append(h3)
            m.append(h4)
            m.append(h5)
            m.append(h6)

            return m

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

    def build_pna( self,  xyz = "tmp.xyz", waters = 0, 
            min_r = 2.0,
            mult_r = 10):
        pna = Molecule.from_xyz( xyz )
        freqs = [ "0.0", "0.0238927", "0.0428227", "0.0773571" ] 

        np.random.seed(args.seed)

        c = Cluster()
        c.append(pna, in_qm = True)
        cnt = 0
        while cnt < waters:
# Random rotation angles
            t1 = np.random.uniform( 0, np.pi/2 )
            t2 = np.random.uniform( 0, np.pi   )
            t3 = np.random.uniform( 0, np.pi/2 )

# random length, rho and tau 
            r =  np.random.uniform( min_r , min_r * mult_r)
            tau =  np.random.uniform( 0, np.pi*2)
            theta =  np.random.uniform( 0,np.pi)

            center = self.polar_to_cartesian( r, tau, theta )

            wat = self.get_mol( center = pna.com + center,
                    mol = "water")

            wat.rotate( t1, t2, t3 )
            wat.res_id = cnt
            for at in wat:
                at.res_id = wat.res_id

            if c.mol_too_close( wat ):
                continue

#We are satisfied with this position, add properties to the water, and rotate them according to t1, t2, t3 so they match the water orientation
            c.append( wat, in_mm = True )
            cnt += 1

        for f_mm in freqs:
            for dist in ["nodist", "dist"]:
                for wat in [ m for m in c if m.in_mm ]:
                    t1, t2, t3 =  wat.get_euler()
                    kwargs_dict = Template().get( *("TIP3P", "HF", "ANOPVDZ",
                        dist == "dist",f_mm ) )
                    for at in wat:
                        Property.add_prop_from_template( at, kwargs_dict )
                    Property.transform_ut_properties( wat.h1.Property, t1,t2,t3 )
                    Property.transform_ut_properties( wat.h2.Property, t1,t2,t3 )
                    Property.transform_ut_properties( wat.o.Property,  t1,t2,t3 )
#Write out QM and MM region separately with properties
                open("pna.mol" ,'w').write(c.get_qm_mol_string(
                    basis= ("ano-1 2 1", "ano-1 3 2 1"),
                    AA = True))
                open("%dmm_%s_%s.pot" %(waters, f_mm, dist ),'w').write(c.get_qmmm_pot_string( AA = True))
                open("tmp.xyz", 'w').write( c.get_xyz_string() )

if __name__ == '__main__':

    A = argparse.ArgumentParser( add_help= True)

#Related to generating molecules with specified parameters
# Right now 2 waters implemented with 6 parameters

    A.add_argument( "-param", action = 'store_true', default = False )
    A.add_argument( "-param_mol", type = str, default = 'water' )



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
            type = str, choices = ["water",
                "methanol" , "ethane"] )
    A.add_argument( "-basis", type = str, nargs = "*",
            default =["ano-1 2 1", "ano-1 3 2 1" ] )
    A.add_argument( "-AA" ,  default = False, action = 'store_true' )

#########################
#      PNA RELATED
#########################

    A.add_argument( '-seed', type = int, default = 111 )

    A.add_argument( '-pna', action = 'store_true', default = False )
    A.add_argument( '-pna_waters', type = int, default = 10 )
    A.add_argument( '-pna_min_r', type = float, default = 0 )
    A.add_argument( '-pna_mult_r', type = float, default = 10 )

########################################################################

    args = A.parse_args()

    opts =  {
       "r"    :{"min":args.r,    "max": args.r_max,    "points":args.r_points,    },
       "tau"  :{"min":args.tau,  "max": args.tau_max,  "points":args.tau_points,  },
       "theta":{"min":args.theta,"max": args.theta_max,"points":args.theta_points,},
       "rho1" :{"min":args.rho1, "max": args.rho1_max, "points":args.rho1_points, },
       "rho2" :{"min":args.rho2, "max": args.rho2_max, "points":args.rho2_points, },
       "rho3" :{"min":args.rho3, "max": args.rho3_max, "points":args.rho3_points, },
             }

    g = Generator( **opts )

    if args.param:
        g.vary_parameters( opts )
        g.gen_mols_param( 
                mol = args.param_mol ,
                basis = args.basis,
                AA = False )
        raise SystemExit

    if args.get_mol:
        mol = g.get_mol( center = [0, 0, 0 ], mol = args.get_mol , AA = args.AA )

    if args.pna:
        g.build_pna( xyz = "pna_opt.xyz", 
                waters = args.pna_waters,
                min_r = args.pna_min_r,
                mult_r = args.pna_mult_r )
                

