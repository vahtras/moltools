#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import molecules, os

def sulfuric_acid():
    """Return geo opt. 
    molecule with sulfur in origo, one oxygen in xz plane"""
    builddir = "build"
    molfile = "sulfur_opt.xyz"
    FILE = os.path.join( os.path.dirname( os.path.realpath( __file__) ) , os.path.join( builddir, molfile ))
    m = molecules.Molecule.from_xyz( FILE, in_AA = True, out_AA = False )
    return m


def paranitro_aniline():
    """Return geo opt. 
    molecule with sulfur in origo, one oxygen in xz plane"""
    builddir = "build"
    molfile = "pna_opt.xyz"
    FILE = os.path.join( os.path.dirname( os.path.realpath( __file__) ) , os.path.join( builddir, molfile ))
    m = molecules.Molecule.from_xyz( FILE, in_AA= True, out_AA = False )
    return m

def tip3p():
    """Return geo opt. 
    molecule with sulfur in origo, one oxygen in xz plane"""
    builddir = "build"
    molfile = "tip3p.xyz"
    FILE = os.path.join( os.path.dirname( os.path.realpath( __file__) ) , os.path.join( builddir, molfile ))
    m = molecules.Water.get_standard( AA = False )
    for ind, at in enumerate( m ):
        at.order = ind + 1

    return m


if __name__ == '__main__':
    main()
