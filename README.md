[![Documentation Status](https://readthedocs.org/projects/dalton-tools/badge/?version=latest)](https://readthedocs.org/projects/dalton-tools/?badge=latest)
[![Build Status](https://travis-ci.org/fishstamp82/dalton_tools.svg?branch=master)](https://travis-ci.org/fishstamp82/dalton_tools)
[![Coverage Status](https://img.shields.io/coveralls/fishstamp82/dalton_tools.svg)](https://coveralls.io/r/fishstamp82/dalton_tools?branch=master)

#Welcome to moltools!

	Code purpose: Wrap DALTON LoProp calculation into convinient functions callable by the Molecule instance using IPython or python scripts.

## Current features:

	Features include obtaining LoProp properties for solvent molecules/ligands or for proteins and polymers that are covalently bonded via the MFCC procedure.

	By integrating with the particles module, applequist equations are directly solvable for a system of classical molecules using damped charges/ dipole-moments directly from QM-obtainable properties.

	For localized Beta, this requires the latest development source of DALTON installed.

## Installation:

`git clone --recursive git@github.com:fishstamp82/moltools.git`

`export PYTHONPATH=$(pwd)/moltools/src:$PYTHONPATH`

> Tip: Export the pythonpath variable in your initrc file of choice in order to have it automatically load.


Execute the following script if you want to run DALTON computations in parallel using HPC clusters.

For the Linköping HPC triolith, execute:

`src/scripts/dalton_run_on_triolith.sh`

For Umeå HPC akka, execute:

`src/scripts/dalton_run_on_akka.sh`


## A quick-start:

Run:

* `ipython`
* `in [1]: from molecules import Water, Cluster`

##### Create a water molecule with oxygen in origo, in atomic units by default
`in [2]: w1 = Water().get_standard()`

##### Create an additional water molecule (atomic units by default)
`in [3]: w2 = Water()`

##### Translate water 2 by 2.5 AU in the z-axis
`in [4]: w2.translate_by_r( [0, 0, 2.5] )`

##### Add them together into a Cluster
`in [5]: c = Cluster( w1, w2 )`

##### You can always make a quick visualization of a Molecule / Cluster
`in [6]: c.plot( copy = True, center = True)`

##### Print a dalton input file using standard basis sets

`in [7]: open( "two_waters_input.mol", 'w').write( c.get_qm_mol_string())`

##### Attach some properties to the waters (The rotation of properties will be taken care of )

`in [8]: c.attach_properties( model = 'tip3p', method = 'HF', basis ='ANOPVDZ' )`

##### output the atomic/ molecular/ cluster propertiy via the .Property keyword, or via the quick-wrapper .p (.d for dipole, .a alpha .etc )

`in [9]: print c.p.a`
[]

******

######Visit [the documentation](http://dalton-tools.readthedocs.org/en/latest) for the API and more tutorials on the source code.
