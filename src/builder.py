#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import molecules, os


def build_sulfuric_acid():
    """Return molecule with sulfur in origo, one hydrogen in xz plane"""
    FILE = os.path.join( os.path.dirname( os.path.realpath( __file__) ) , "build/sulfur_opt.xyz" )
    m = molecules.Molecule.from_xyz( FILE )
    return m

        

if __name__ == '__main__':
    main()

