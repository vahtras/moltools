#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from .molecules import Molecule 
from .polymer import Monomer

def pmma_monomer( t = 0 ):
    """Return pmmma monomer building block as defined by SMILES
    format obtained in avogadro

    4 differen versions exist

    0 is standard from smiles( gives bas CH-3 overlap between neighbours
    1 is where CB1 is rotated pi/4.5 degrees around CA-CB1 bond
    2 is where also hydrogens in CD methyl group are rotated by np/2.5
    3 is the geometry optimized monomer from B3LYP/6-31+G* with the C-C bond beeing the average of the previous included and next
    included monomer in the optimization
    
    """
    builddir = 'build'
    molfile = 'pmma_monomer%d.pdb' %t
    FILE = os.path.join( os.path.dirname( os.path.realpath( __file__) ) , os.path.join( builddir, molfile ))

    m = Monomer.from_pdb( FILE, in_AA = True, out_AA = True )
#Set so that atoms are connected to this molecule
    for at in m:
        at.Molecule = m
    m._mono_name = "pmma"
    m._r = 1.46
    m._rn = 1.46

#This is for geometry optimized monomer, which as slighlty larger CA-HN bond
    if t == 3 :
        m._r = 1.577
        m._rn = 1.570
    m._angle = 104.5
    m._dihedral = 180.0
    return m


def sulfuric_acid():
    """Return geo opt. 
    molecule with sulfur in origo, one oxygen in xz plane"""
    builddir = "build"
    molfile = "sulfur_opt.xyz"
    FILE = os.path.join( os.path.dirname( os.path.realpath( __file__) ) , os.path.join( builddir, molfile ))
    m = Molecule.from_xyz( FILE, in_AA = True, out_AA = False )
    return m


def paranitro_aniline():
    """Return geo opt. 
    molecule with sulfur in origo, one oxygen in xz plane"""
    builddir = "build"
    molfile = "pna_opt.xyz"
    FILE = os.path.join( os.path.dirname( os.path.realpath( __file__) ) , os.path.join( builddir, molfile ))
    m = Molecule.from_xyz( FILE, in_AA= True, out_AA = False )
    return m

def tip3p():
    """Return geo opt. 
    molecule with sulfur in origo, one oxygen in xz plane"""
    builddir = "build"
    molfile = "tip3p.xyz"
    FILE = os.path.join( os.path.dirname( os.path.realpath( __file__) ) , os.path.join( builddir, molfile ))
    m = Water.get_standard( AA = False )
    for ind, at in enumerate( m ):
        at.order = ind + 1

    return m


if __name__ == '__main__':
    main()
