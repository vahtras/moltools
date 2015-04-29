[![Documentation Status](https://readthedocs.org/projects/dalton-tools/badge/?version=latest)](https://readthedocs.org/projects/dalton-tools/?badge=latest)
[![Build Status](https://travis-ci.org/fishstamp82/dalton_tools.svg?branch=master)](https://travis-ci.org/fishstamp82/dalton_tools)
[![Coverage Status](https://img.shields.io/coveralls/fishstamp82/dalton_tools.svg)](https://coveralls.io/r/fishstamp82/dalton_tools?branch=master)

#Welcome to moltools!

Code intended to aid in the analysis of calculations of molecules.

    Right now porting from dalton_tools.git

## Current features:

MM properties for QMMM available for some molecules as templates.

Generation of a polymer from a monomer unit.

## Installation:

`git clone git@github.com:fishstamp82/moltools.git`

`export PYTHONPATH=$(pwd)/moltools/src:$PYTHONPATH`

Execute the following script so that private dependencies are uncommented.

`/src/scripts/make_master_pass_test.sh`

> Tip: Export the pythonpath variable in your initrc file of choice in order to have it automatically load.

_______

## Quick-start:


Run:

* `ipython`
* `in [1]: from molecules import Water, Cluster`

##### Create a water molecule with oxygen in origo, in atomic units
`in [2]: w1 = Water().get_standard( AA = False )`

##### Create an additional water molecule 
`in [3]: w2 = Water().get_standard( AA = False )`

##### Translate water 2 by 2.5 AU in the z-axis
`in [4]: w2.translate_by_r( [0, 0, 2.5] )`

##### Add them together into a Cluster
`in [5]: c = Cluster( w1, w2)`

##### You can always make a quick visualization of a Molecule / Cluster
`in [6]: c.plot( copy = True, center = True)`

##### Print a dalton input file using standard basis sets

`in [7]: open( "two_waters_input.mol", 'w').write( c.get_qm_mol_string())`


**********

######Visit [the documentation](http://dalton-tools.readthedocs.org/en/latest) for the API and more tutorials on the source code.
