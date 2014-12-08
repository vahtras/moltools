[![Documentation Status](https://readthedocs.org/projects/dalton-tools/badge/?version=latest)](https://readthedocs.org/projects/dalton-tools/?badge=latest)
[![Build Status](https://travis-ci.org/fishstamp82/dalton_tools.svg?branch=master)](https://travis-ci.org/fishstamp82/dalton_tools)
[![Coverage Status](https://img.shields.io/coveralls/fishstamp82/dalton_tools.svg)](https://coveralls.io/r/fishstamp82/dalton_tools?branch=master)

#Welcome to Dalton Tools!

Code intended to aid in the analysis of calculations of water molecules.
___

## Installation:

`git clone git@github.com:fishstamp82/dalton_tools.git dalton_tools`

`export PYTHONPATH=$(pwd)/dalton_tools:$PYTHONPATH`

> Tip: Add the path to the dalton_tools directory to your initrc file of choice in order to have it automatically load.

_______

## Quick-start:


Run:

* \>>>`ipython`
* \>>>`from use_generator import Generator`
* \>>>`from molecules import Cluster`

#### Create a water molecule with oxygen in origo, in atomic units
\>>> `w1 = Generator().get_mol( center = [0, 0, 0], mol = "water", AA = False )`

#### Create an additional water molecule with the oxygen located at z = 2.5, in atomic units
\>>> `w2 = Generator().get_mol( center = [0, 0, 2.5], mol = "water", AA = False )`

#### Add them together into a Cluster

\>>> `Cluster.add_mol( w1 )`

\>>> `Cluster.add_mol( w2 )`

#### Write the resulting configuration to a DALTON .mol file

\>>> `open( "two_waters_input.mol", 'w').write( Cluster.get_mol_string( AA = False ))`


**********

######Visit [the documentation](http://dalton-tools.readthedocs.org/en/latest) for the API and tutorials on the source code.
