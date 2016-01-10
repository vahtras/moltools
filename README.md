
# Running a second order Applequist interaction between two water molecules
This tutorial will describe the steps needed to calculate the classical polarizability and hyperpolarizability of two water molecules interacting, using the localized properties of each molecule obtained seperately from quantum mechanics.


The two interacting water molecules will be of the TIP3P model, separated at 5 bohr.

------------


# Step 1
## Obtaining the DALTON source code

The DALTON source can be obtained from http://www.daltonprogram.org/. It uses an LGPL license.

-----------

To run a DALTON calculation once it is installed, execute the dalton runscript:

```bash
>>> dalton -get input.dal molecule.mol
```

To see a list of options such as MPI core usage and temporary directory setup, execute:

```bash
>>> dalton
```

#### The following input file should be used for a quadratic response calculation

# dalton.inp
**DALTON INPUT
.RUN RESPONSE
.DIRECT
.PARALLELL
**WAVE FUNCTION
.HF
**RESPONSE
.PROPAV
XDIPLEN
.PROPAV
YDIPLEN
.PROPAV
ZDIPLEN
*QUADRATIC
.DIPLEN
**END OF DALTON INPUT
#### The following molecule file specifies one TIP3P water molecule in the xz-plane with its dipole pointing in the positive z-axis direction
## molecule.inp
ATOMBASIS
Comment 1
Comment 2
Atomtypes=2 Charge=0 Nosymm
Charge=8.0 Atoms=1 Basis=ano-1 4 3 1
O       0.00000   0.00000   0.00000
Charge=1.0 Atoms=2 Basis=ano-1 2
H       1.43043   0.00000   1.10716
H      -1.43043   0.00000   1.10716
## Step 2 obtain properties from LoProp

Once the input_molecule.tar.gz file is obtained, obtain the source code to do the Applequist calculations in an easy way:

```bash
>>> git clone --recursive https://github.com/vahtras/loprop.git
```

With the loprop/loprop.py script in your PATH, execute 

```bash
>>> loprop.py -l 1 -a 2 -B 1
```


```python

```
