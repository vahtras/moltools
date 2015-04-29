[![Documentation Status](https://readthedocs.org/projects/dalton-tools/badge/?version=latest)](https://readthedocs.org/projects/dalton-tools/?badge=latest)
[![Build Status](https://travis-ci.org/fishstamp82/dalton_tools.svg?branch=master)](https://travis-ci.org/fishstamp82/dalton_tools)
[![Coverage Status](https://img.shields.io/coveralls/fishstamp82/dalton_tools.svg)](https://coveralls.io/r/fishstamp82/dalton_tools?branch=master)

#Welcome to moltools!

Code intended to aid in the analysis of calculations of water molecules.

    Right now porting from dalton_tools.git

## Installation:

`git clone git@github.com:fishstamp82/moltools.git moltools`

`export PYTHONPATH=$(pwd)/moltools:$PYTHONPATH`

> Tip: Add the path to the moltools directory to your initrc file of choice in order to have it automatically load.

_______

## Quick-start:


Run:

* \>>>`ipython`
* `in [1]: from molecules import Water, Cluster`

##### Create a water molecule with oxygen in origo, in atomic units
\>>> `w1 = Water().get_standard( AA = False )

##### Create an additional water molecule 
\>>> `w2 = Water().get_standard( AA = False )

##### Translate water 2 by 2.5 AU in the z-axis
\>>> `w2.translate_by_r( [0, 0, 2.5] )


##### Add them together into a Cluster

\>>> `c = Cluster( w1, w2)`

##### Set up which waters are in the QM region and which in the MM region

\>>> `open( "two_waters_input.mol", 'w').write( c.get_qm_mol_string( AA = False ))`


**********

######Visit [the documentation](http://dalton-tools.readthedocs.org/en/latest) for the API and more tutorials on the source code.
