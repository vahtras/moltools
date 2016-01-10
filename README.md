[![Documentation Status](https://readthedocs.org/projects/moltools/badge/?version=docs)](https://readthedocs.org/projects/moltools/?badge=docs)
[![Build Status](https://travis-ci.org/fishstamp82/moltools.svg?branch=master)](https://travis-ci.org/fishstamp82/moltools)
[![Coverage Status](https://coveralls.io/repos/fishstamp82/moltools/badge.svg?branch=master&service=github)](https://coveralls.io/github/fishstamp82/moltools?branch=master)
[![DOI](https://zenodo.org/badge/19666/fishstamp82/moltools.svg)](https://zenodo.org/badge/latestdoi/19666/fishstamp82/moltools)


#Welcome to moltools!

	Code purpose: Wrap DALTON LoProp calculation into convinient functions
	callable by the Molecule instance using IPython or python scripts.

## Current features:

	Features include obtaining LoProp properties for solvent molecules/ligands 
	or for proteins and polymers that are covalently bonded via the MFCC procedure.

	By integrating with the particles module, applequist equations are directly 
	solvable for a system of classical molecules using damped 
	charges/dipole-moments 	directly from QM-obtainable properties.

	For localized Beta, this requires the latest development source of 
	DALTON installed.

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
`in [3]: w2 = Water.get_standard()`

##### Translate water 2 by 2.5 AU in the z-axis
`in [4]: w2.translate_by_r( [0, 0, 2.5] )`

##### Add them together into a Cluster
`in [5]: c = Cluster( w1, w2 )`

##### You can always make a quick visualization of a Molecule / Cluster
`in [6]: c.plot()`

##### Attach some properties to the waters (The rotation of properties will be taken care of )

`#See template.py for all available templates`

`in [8]: c.attach_properties( model = 'tip3p', method = 'HF', basis ='ANOPVDZ' )`

##### Output the atomic/ molecular/ cluster propertiy via the .Property keyword, or via the quick-wrapper .p (.d for dipole, .a alpha .etc )

`in [9]: print c.p.a`
`[ 15.02184   0.        0.       11.48016   0.       13.72182]`

##### Calculate each waters properties from ab-initio using DALTON, and put those properties on each atom using LoProp in one step:

`In [10]: c.props_from_qm( tmpdir = '/tmp', dalpath = $PATH_TO_DALTON_SCRIPT )`

##### If the dalton version is the development master branch, localized hyperpolarizabilities are obtainable:

`In [11]: c.props_from_qm( method = 'b3lypqua', tmpdir = '/tmp', dalpath = $PATH_TO_DALTON_SCRIPT )`

******


## Extra features:

	These include uncommenting "#from mayavi import mlab" in 
	src/pdbreader.py and an installation of mayavi2.
	This enables plotting of the beta tensor around molecules and clusters.


######Visit [the documentation](http://moltools.readthedocs.org/en/latest) for the API and more tutorials on the source code. Work in progress and most stuff are outdated.


