[![Documentation Status](https://readthedocs.org/projects/dalton-tools/badge/?version=latest)](https://readthedocs.org/projects/dalton-tools/?badge=latest)
[![Build Status](https://travis-ci.org/fishstamp82/dalton_tools.svg?branch=master)](https://travis-ci.org/fishstamp82/dalton_tools)
[![Coverage Status](https://img.shields.io/coveralls/fishstamp82/dalton_tools.svg)](https://coveralls.io/r/fishstamp82/dalton_tools?branch=master)

#Welcome to moltools!

	Code intended to aid in the analysis of calculations of molecules.

	Right now porting from dalton_tools.git

## Current features:

	MM properties for QMMM available for some molecules as templates.

	MM properties for molecules and proteins calculatable directly from script* or interactive python. Also works for High-performance computers (tested on triolith and hpc in Umeå).

	Generation of a polymer from a monomer unit.

	*Needs vahtras/loprop.git to be in PYTHONPATH

## Installation:

`git clone git@github.com:fishstamp82/moltools.git`

`export PYTHONPATH=$(pwd)/moltools/src:$PYTHONPATH`

Execute the following script under the scripts subdir so that private dependencies are uncommented.

`src/scripts/make_master_pass_test.sh`

For Linköping HPC, execute:

`src/scripts/dalton_run_on_triolith.sh`

For Umeå HPC, execute:

`src/scripts/dalton_run_on_akka.sh`

> Tip: Export the pythonpath variable in your initrc file of choice in order to have it automatically load.


## Quick-start:

Run:

* `ipython`
* `in [1]: from molecules import Water, Cluster`

##### Create a water molecule with oxygen in origo, in atomic units
`in [2]: w1 = Water().get_standard( AA = False )`

##### Create an additional water molecule (atomic units by default)
`in [3]: w2 = Water()

##### Translate water 2 by 2.5 AU in the z-axis
`in [4]: w2.translate_by_r( [0, 0, 2.5] )`

##### Add them together into a Cluster
`in [5]: c = Cluster( w1, w2 )`

##### You can always make a quick visualization of a Molecule / Cluster
`in [6]: c.plot( copy = True, center = True)`

##### Print a dalton input file using standard basis sets

`in [7]: open( "two_waters_input.mol", 'w').write( c.get_qm_mol_string())`

##### Attach some properties to the waters (The rotation of properties will be taken care of )

`in [8]: c.attach_properties( model = 'tip3p', method = 'HF', basis ='ANOPVDZ' )

##### visualize atomic/ molecular/ cluster properties via the .Property keyword

`in [9]: c.p

**********

######Visit [the documentation](http://dalton-tools.readthedocs.org/en/latest) for the API and more tutorials on the source code.
