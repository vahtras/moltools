__all__ = [ 'Generator' ]



import numpy as np
import molecules 

a0 = 0.52917721092
class Generator( dict ):
    """
    Used to create molecules, write dalton .mol files 
    using -param for study with use_calculator.py

    water currently implemented only

    plans to implement methanol

    """
    def __init__(self, *args, **kwargs):

#This waater is TIP3P model,
        self[ ("water", "tip3p", "a_hoh", "degree") ] = 104.52
        self[ ("water", "tip3p", "r_oh", "AA") ] = 0.9572

#This waater is SPC model,
        self[ ("water", "spc", "a_hoh", "degree") ] = 109.47
        self[ ("water", "spc", "r_oh", "AA") ] = 1.0

        self[ ("methanol", "gas_opt", "r_oh", "AA" ) ] = 0.967
        self[ ("methanol", "gas_opt", "r_co", "AA" ) ] = 1.428
        self[ ("methanol", "gas_opt", "r_ch", "AA" ) ] = 1.098

        self[ ("methanol", "gas_opt", "a_coh", "degree" ) ] = 107.16
        self[ ("methanol", "gas_opt", "a_hch", "degree" ) ] = 109.6
        self[ ("methanol", "gas_opt", "a_hco", "degree" ) ] = 109.342

        self[ ("methanol", "gas_opt", "d_hcoh", "h4", "degree" ) ] =  60.0
        self[ ("methanol", "gas_opt", "d_hcoh", "h5", "degree" ) ] = -60.0
        self[ ("methanol", "gas_opt", "d_hcoh", "h6", "degree" ) ] =  180.0

        
#Default options for water
        for val in ["r", "tau", "theta", "rho1", "rho2", "rho3", ]:
            self[ ( val, 'min') ]    = 0.0
            self[ ( val, 'max') ]    = 0.0
            self[ ( val, 'points') ] = 1
        self[ ( 'r', 'min') ]    = 5.0
        self[ ( 'r', 'max') ]    = 10.0
        self[ ( 'r', 'points') ] = 1

# Set by default all parameters to False
        for val in ["r", "tau", "theta", "rho1", "rho2", "rho3", ]:
            self[ ( val, "active" ) ]  = False

    @staticmethod
    def get_pe_b3lyp_dal( co = 1.0, AA = True, max_l = 2, sites = 3):
        r_order = max_l + 1
        if AA:
            aa = "AA"
        else:
            aa = "AU"
            co /= a0
        st = """**DALTON INPUT
.RUN WAVE FUNCTION
.DIRECT
.PARALLELL
.PEQM
*PEQM
.BORDER
REDIST -%d %.1f %s %d
**WAVE FUNCTION
.DFT
B3LYP
**END OF DALTON INPUT""" % (max_l+1, co, aa, sites)
        return st

    @staticmethod
    def get_qmmm_b3lyp_dal( damp = False):
        if damp:
            damp = "\n.DAMP"
        else:
            damp = ""
        st = """**DALTON INPUT
.RUN RESPONSE
.DIRECT
.PARALLELL
*QMMM
.QMMM%s
**WAVE FUNCTION
.DFT
B3LYP
**END OF DALTON INPUT""" % damp
        return st

    @staticmethod
    def get_lin_dal( _string ):
        if _string == 'hflin':
            return Generator.get_hflin_dal()
        elif _string == 'b3lyplin':
            return Generator.get_b3lyplin_dal()
        elif _string == 'camb3lyplin':
            return Generator.get_camb3lyplin_dal()
        elif _string == 'ccsdlin':
            return Generator.get_ccsdlin_dal()

    @staticmethod
    def get_hf_imag_dal( freq = ("0.0",), functional = 'B3PW91' ):
        _string = """**DALTON INPUT
.RUN RESPONSE
.DIRECT
.PARALLELL
**WAVE FUNCTION
.HF
**RESPONSE
*ABSORP
.ALPHA
.IMAG F
.FREQUE
"""
        _string += str(len(freq))  + '\n'
        freqs = " ".join( map( str, freq ) )
        _string += freqs
        _string += '\n'
        _string += "**END OF DALTON INPUT\n"
        return _string



    @staticmethod
    def get_b3lyplin_dal():
        return """**DALTON INPUT
.RUN RESPONSE
.DIRECT
.PARALLELL
**WAVE FUNCTION
.DFT
B3LYP
.INTERFACE
**INTEGRAL
.DIPLEN
.SECMOM
**RESPONSE
.PROPAV
XDIPLEN
.PROPAV
YDIPLEN
.PROPAV
ZDIPLEN
*LINEAR
.DIPLEN
**END OF DALTON INPUT""" 

    @staticmethod
    def get_b3lypqua_dal( ):
        return """**DALTON INPUT
.RUN RESPONSE
.DIRECT
.PARALLELL
**WAVE FUNCTION
.DFT
B3LYP
.INTERFACE
**INTEGRAL
.DIPLEN
.SECMOM
**RESPONSE
.PROPAV
XDIPLEN
.PROPAV
YDIPLEN
.PROPAV
ZDIPLEN
*QUADRATIC
.QLOP
.DIPLEN
**END OF DALTON INPUT""" 

    @staticmethod
    def get_hflin_freq_dal( freq = "0.0", au = True, nm = False  ):
        _string = """**DALTON INPUT
.RUN RESPONSE
.DIRECT
.PARALLELL
**WAVE FUNCTION
.HF
**INTEGRAL
.DIPLEN
.SECMOM
**RESPONSE
.PROPAV
XDIPLEN
.PROPAV
YDIPLEN
.PROPAV
ZDIPLEN
*LINEAR
.DIPLEN
.FREQUE
 1
 %s
""" %( freq )
        _string += "**END OF DALTON INPUT\n"
        return _string
    @staticmethod
    def get_hfqua_dal( ):
        return """**DALTON INPUT
.RUN RESPONSE
.DIRECT
.PARALLELL
**WAVE FUNCTION
.HF
.INTERFACE
**INTEGRAL
.DIPLEN
.SECMOM
**RESPONSE
.PROPAV
XDIPLEN
.PROPAV
YDIPLEN
.PROPAV
ZDIPLEN
*QUADRATIC
.QLOP
.DIPLEN
**END OF DALTON INPUT""" 

    @staticmethod
    def get_hflin_dal( ):
        return """**DALTON INPUT
.RUN RESPONSE
.DIRECT
.PARALLELL
**WAVE FUNCTION
.HF
.INTERFACE
**INTEGRAL
.DIPLEN
.SECMOM
**RESPONSE
.PROPAV
XDIPLEN
.PROPAV
YDIPLEN
.PROPAV
ZDIPLEN
*LINEAR
.DIPLEN
**END OF DALTON INPUT""" 

    @staticmethod
    def get_b3lyplin_freq_dal( freq = "0.0", au = True, nm = False  ):
        _string = """**DALTON INPUT
.RUN RESPONSE
.DIRECT
.PARALLELL
**WAVE FUNCTION
.DFT
B3LYP
.INTERFACE
**INTEGRAL
.DIPLEN
.SECMOM
**RESPONSE
.PROPAV
XDIPLEN
.PROPAV
YDIPLEN
.PROPAV
ZDIPLEN
*LINEAR
.DIPLEN
.FREQUE
 1
 %s
""" %( freq )
        _string += "**END OF DALTON INPUT\n"
        return _string

    @staticmethod
    def get_camb3lyplin_dal():
        return """**DALTON INPUT
.RUN RESPONSE
.DIRECT
.PARALLELL
**WAVE FUNCTION
.DFT
CAMB3LYP
.INTERFACE
**INTEGRAL
.DIPLEN
.SECMOM
**RESPONSE
.PROPAV
XDIPLEN
.PROPAV
YDIPLEN
.PROPAV
ZDIPLEN
*LINEAR
.DIPLEN
**END OF DALTON INPUT""" 



    @staticmethod
    def get_camb3lypqua_dal( ):
        return """**DALTON INPUT
.RUN RESPONSE
.DIRECT
.PARALLELL
**WAVE FUNCTION
.DFT
CAMB3LYP
.INTERFACE
**INTEGRAL
.DIPLEN
.SECMOM
**RESPONSE
.PROPAV
XDIPLEN
.PROPAV
YDIPLEN
.PROPAV
ZDIPLEN
*QUADRATIC
.QLOP
.DIPLEN
**END OF DALTON INPUT""" 



    @staticmethod
    def get_camb3lyplin_freq_dal( freq = "0.0", au = True, nm = False  ):
        _string = """**DALTON INPUT
.RUN RESPONSE
.DIRECT
.PARALLELL
**WAVE FUNCTION
.DFT
CAMB3LYP
.INTERFACE
**INTEGRAL
.DIPLEN
.SECMOM
**RESPONSE
.PROPAV
XDIPLEN
.PROPAV
YDIPLEN
.PROPAV
ZDIPLEN
*LINEAR
.DIPLEN
.FREQUE
 1
 %s
""" %( freq )
        _string += "**END OF DALTON INPUT\n"
        return _string



    @staticmethod
    def get_ccsdlin_dal():
        return """**DALTON INPUT
.RUN RESPONSE
**INTEGRALS
.DIPLEN
.SECMOM
**WAVE FUNCTION
.CC
*CC INPUT
.CCSD
*CCFOP
.DIPMOM
*CCLR
.DIPOLE
**END OF DALTON INPUT
"""
    @staticmethod
    def get_ccsdqua_dal():
        return """
**DALTON INPUT
.RUN RESPONSE
**INTEGRALS
.DIPLEN
.SECMOM
**WAVE FUNCTION
.CC
*CC INPUT
.CCSD
*CCFOP
.DIPMOM
*CCLR
.DIPOLE
*CCQR
.DIPOLE
**END OF DALTON INPUT
"""


    def get_mol( self, 
            center = [0,0,0], 
            mol = "water", 
            model = "tip3p",
            AA = False ):
        """return molecule in center, all molecules have different definition
        of euler angles

        for water place O in origo
        for methanol place C=O bond in origo
        
        """

        if mol == "water":
#Geometrical parameters, dependent om model
            if model == "tip3p":
                r_oh = self[ ("water", "tip3p", "r_oh", "AA") ]
                a_hoh = self[ ("water", "tip3p", "a_hoh","degree") ]

            if model == "spc":
                r_oh = self[ ("water", "spc", "r_oh", "AA") ]
                a_hoh = self[ ("water", "spc", "a_hoh","degree") ]

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

            h1 = molecules.Atom( **{ "AA" : AA,
                "x" : xh1,
                "y" : yh1,
                "z" : zh1,
                "element" : "H"} )
            h2 = molecules.Atom( **{ "AA" : AA,
                "x" : xh2,
                "y" : yh2,
                "z" : zh2,
                "element" : "H"} )
            o = molecules.Atom( **{ "AA" : AA,
                "x" : xo,
                "y" : yo,
                "z" : zo,
                "element" : "O"} )

            w = molecules.Water( AA = AA)
            w.append( o )
            w.append( h1 )
            w.append( h2 )
            
            return w

        elif mol == "methanol":

            r_co = self[ ("methanol", "gas_opt", "r_co", "AA" )]
            r_oh = self[ ("methanol", "gas_opt", "r_oh", "AA" )]
            r_ch = self[ ("methanol", "gas_opt", "r_ch", "AA" )]

            a_coh = self[ ("methanol", "gas_opt", "a_coh", "degree" ) ]
            #a_hch = self[ ("methanol","gas_opt",  "a_hch", "degree" ) ]
            a_hco = self[ ("methanol", "gas_opt", "a_hco", "degree" ) ]

            a_coh *= np.pi / 180
            a_hco *= np.pi / 180

            d_hcoh_4 = self[ ("methanol","gas_opt",  "d_hcoh", "h4", "degree" ) ]
            d_hcoh_4 *= np.pi / 180
            d_hcoh_5 = self[ ("methanol","gas_opt",  "d_hcoh", "h5", "degree" ) ]
            d_hcoh_5 *= np.pi / 180
            d_hcoh_6 = self[ ("methanol","gas_opt",  "d_hcoh", "h6", "degree" ) ]
            d_hcoh_6 *= np.pi / 180

            if not AA:
                r_co, r_oh, r_ch = r_co/a0, r_oh/a0, r_ch/a0

            c1 = molecules.Atom( **{"x":0, "y":0, "z":-r_co/2, "AA": AA, "element":"C" } )
            o2 = molecules.Atom( **{"x":0, "y":0, "z": r_co/2, "AA": AA, "element":"O" } )

            h3 = molecules.Atom( **{"x":r_oh*np.cos( a_coh-np.pi/2),
                "y":0,
                "z":r_oh*np.sin( a_coh-np.pi/2) + r_co/2,
                "AA": AA, "element":"H" } )

            h4 = molecules.Atom( **{"x": r_ch*np.sin( a_hco ) * np.cos( d_hcoh_4 ),
                "y": r_ch*np.sin( a_hco) * np.sin( d_hcoh_4 ),
                "z": r_ch*np.cos( a_hco) - r_co/2 ,
                "AA": AA, "element":"H" } )
            h5 = molecules.Atom( **{"x": r_ch*np.sin( a_hco ) * np.cos( d_hcoh_5 ),
                "y": r_ch*np.sin( a_hco) * np.sin( d_hcoh_5 ),
                "z": r_ch*np.cos( a_hco) - r_co/2 ,
                "AA": AA, "element":"H" } )
            h6 = molecules.Atom( **{"x": r_ch*np.sin( a_hco ) * np.cos( d_hcoh_6 ),
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



    def gen_mols_param(self, mol = "water", 
            model = 'tip3p',
            basis = ["ano-1 2", "ano-1 4 3 1"],
            AA = True,
            worst = False):
        r = np.linspace( self[ ('r', 'min')] , self[ ('r', 'max')] ,
            self[ ('r', 'points' ) ]  )
        tau = np.linspace( self[ ('tau', 'min')] , self[ ('tau', 'max')] ,
            self[ ('tau', 'points' ) ] )
        theta = np.linspace( self[ ('theta', 'min')] , self[ ('theta', 'max')] ,
            self[ ('theta', 'points' )  ] )
        rho1 = np.linspace( self[ ('rho1', 'min')], self[ ('rho1', 'max')],
            self[ ('rho1', 'points' )  ] )
        rho2 = np.linspace( self[ ('rho2', 'min')], self[ ('rho2', 'max')],
            self[ ('rho2', 'points' )  ] )
        rho3 = np.linspace( self[ ('rho3', 'min')], self[ ('rho3', 'max')],
            self[ ('rho3', 'points' )  ] )

        
        if model == 'tip3p':
            r_oh = self[ ("water", 'tip3p', "r_oh", "AA") ]
            a_hoh = np.pi * self[ ("water", 'tip3p', "a_hoh", "degree" )] / 180.0
        else:
            r_oh = self[ ("water", 'tip3p', "r_oh", "AA") ]
            a_hoh = np.pi * self[ ("water", 'tip3p', "a_hoh", "degree" )] / 180.0

        for i in r:
            for j in tau:
                for k in theta:
                    for l in rho1:
                        for m in rho2:
                            for n in rho3:
                                w1 = molecules.Water.get_standard( AA = AA)
                                w1.t( -w1.o.r )
                                if worst:
                                    w1 = self.get_mol( [0, 0, 0], 
                                            mol = mol,
                                            model = model, AA = AA)
                                    w1.populate_bonds()
                                    w1.populate_angles()
                                    w1.h1.scale_angle( 0.988 )
                                    w1.h1.scale_bond( 0.985 )
                                    w1.h2.scale_bond( 1.015 )
                                    w1.inv_rotate()
                                w2 = molecules.Water.get_standard( AA = AA )
                                w2.t( -w2.o.r )
                                x, y, z = self.polar_to_cartesian( i, j, k )

                                w2.rotate( l, m, n )
                                w2.t( np.array( [x, y, z]) )

                                name = ""
                                name += "-".join( map( str, ["%3.2f"%i, "%3.2f"%j, "%3.2f"%k, "%3.2f"%l, "%3.2f"%m, "%3.2f"%n] ) )
                                name += ".mol"

                                c = molecules.Cluster( w1, w2 )
                                tmp_mol = c.get_mol_string( basis = tuple(basis))
                                f_ = open(name, 'w')
                                f_.write( tmp_mol )
        return 0

    def vary_parameters( self, opts ):
        """Given two parameters, e.g. r and theta, keeps all other static
        param_list should be list of strings of parameters
        ["r":{"min": 2, "max":5, "points": 10}, "rho1" , ... ]

        Has sane defaults, but can be overrided by passing arguments to 
        main program as:

        -r_min 5
        -r_max 10
        -r_points 10

        Which overrides defaults 

        """
        for val in opts:
            self[ (val, 'active') ] = True
            self[ (val, 'min') ] = opts[val][ "min" ]
            self[ (val, 'max') ] = opts[val][ "max" ]
            self[ (val, 'points') ] = opts[val][ "points" ]

    def polar_to_cartesian(self, r, tau, theta):
        x, y, z = r* np.sin( theta )*np.cos( tau ) \
               , r* np.sin(  theta )*np.sin( tau )  \
               , r* np.cos(  theta ) 

        return x , y , z

    def one_mol_gen(self, mol = 'water', model = 'tip3p',):
        """
        Only implemented for water so far"""


        if mol == "water":
            d = self[ ("r_oh_dev", "max") ]
            p = self[ ("r_oh_dev", "points") ]
            r_d =  0.01*np.linspace( -d, d, p )

            d = self[ ("theta_hoh_dev", "max") ]
            p = self[ ("theta_hoh_dev", "points") ]
            theta_d =  0.01*np.linspace( -d, d, p )

            #a_hoh = self[ ( mol, model, "a_hoh", "degree" ) ] *np.pi/180
            #r_oh = self[ ( mol, model, "r_oh", "AA" ) ]

            for i in r_d:
                for j in r_d:
                    for k in theta_d:
                        scale_bond1 = 1 + i
                        scale_bond2 = 1 + j
                        scale_angle = 1 + k
                        names = map( lambda x:"%.3f"%x, [i, j, k] )
                        w = self.get_mol( mol = mol, model = model)
                        w.populate_bonds() ; w.populate_angles()
                        w.h1.scale_bond( scale_bond1 )
                        w.h2.scale_bond( scale_bond2 )
                        w.h1.scale_angle( scale_angle )
                        w.inv_rotate()
                        open( "_".join([model]+names) + ".mol",'w').write(w.get_mol_string())
    def build_pna( self,  xyz = "tmp.xyz", waters = 0,
            min_r = 2.0,
            mult_r = 10,
            seed = 111 ):
        pna = Molecule.from_xyz( xyz )
        freqs = [ "0.0", "0.0238927", "0.0428227", "0.0773571" ] 

        np.random.seed( seed )

        c = molecules.Cluster()
        c.add_mol(pna, in_qm = True)
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
            wat._res_id = cnt

            if c.mol_too_close( wat ):
                continue

#We are satisfied with this position, add properties to the water, and rotate them according to t1, t2, t3 so they match the water orientation
            c.add_mol( wat, in_mm = True )
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
                open("%dmm_%s_%s.pot" %(waters, f_mm, dist ),'w').write(c.get_qmmm_pot_string( in_AA = True ))
                open("tmp.xyz", 'w').write( c.get_xyz_string() )



