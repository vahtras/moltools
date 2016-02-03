import unittest
import numpy as np

import warnings
warnings.simplefilter('error')
from nose.plugins.attrib import attr
from moltools import Water
from moltools import utilz

HF_FILE = """


     ************************************************************************
     *************** Dalton - An Electronic Structure Program ***************
     ************************************************************************

    This is output from DALTON 2016.alpha
   ----------------------------------------------------------------------------
    NOTE:
     
    Dalton is an experimental code for the evaluation of molecular
    properties using (MC)SCF, DFT, CI, and CC wave functions.
    The authors accept no responsibility for the performance of
    the code or for the correctness of the results.
     
    The code (in whole or part) is provided under a licence and
    is not to be reproduced for further distribution without
    the written permission of the authors or their representatives.
     
    See the home page "http://daltonprogram.org" for further information.
     
    If results obtained with this code are published,
    the appropriate citations would be both of:
     
       K. Aidas, C. Angeli, K. L. Bak, V. Bakken, R. Bast,
       L. Boman, O. Christiansen, R. Cimiraglia, S. Coriani,
       J. Cukras, P. Dahle, E. K. Dalskov, U. Ekstroem,
       T. Enevoldsen, J. J. Eriksen, P. Ettenhuber, B. Fernandez,
       L. Ferrighi, H. Fliegl, L. Frediani, K. Hald, A. Halkier,
       C. Haettig, H. Heiberg, T. Helgaker, A. C. Hennum,
       H. Hettema, E. Hjertenaes, S. Hoest, I.-M. Hoeyvik,
       M. F. Iozzi, B. Jansik, H. J. Aa. Jensen, D. Jonsson,
       P. Joergensen, M. Kaminski, J. Kauczor, S. Kirpekar,
       T. Kjaergaard, W. Klopper, S. Knecht, R. Kobayashi, H. Koch,
       J. Kongsted, A. Krapp, K. Kristensen, A. Ligabue,
       O. B. Lutnaes, J. I. Melo, K. V. Mikkelsen, R. H. Myhre,
       C. Neiss, C. B. Nielsen, P. Norman, J. Olsen,
       J. M. H. Olsen, A. Osted, M. J. Packer, F. Pawlowski,
       T. B. Pedersen, P. F. Provasi, S. Reine, Z. Rinkevicius,
       T. A. Ruden, K. Ruud, V. Rybkin, P. Salek, C. C. M. Samson,
       A. Sanchez de Meras, T. Saue, S. P. A. Sauer,
       B. Schimmelpfennig, K. Sneskov, A. H. Steindal,
       K. O. Sylvester-Hvid, P. R. Taylor, A. M. Teale,
       E. I. Tellgren, D. P. Tew, A. J. Thorvaldsen, L. Thoegersen,
       O. Vahtras, M. A. Watson, D. J. D. Wilson, M. Ziolkowski
       and H. Agren,
       "The Dalton quantum chemistry program system",
       WIREs Comput. Mol. Sci. 2013. (doi: 10.1002/wcms.1172)
    
    and
    
       Dalton, a Molecular Electronic Structure Program,
       Release Dalton2016.alpha (2015), see http://daltonprogram.org
   ----------------------------------------------------------------------------

    Authors in alphabetical order (major contribution(s) in parenthesis):

  Kestutis Aidas,           Vilnius University,           Lithuania   (QM/MM)
  Celestino Angeli,         University of Ferrara,        Italy       (NEVPT2)
  Keld L. Bak,              UNI-C,                        Denmark     (AOSOPPA, non-adiabatic coupling, magnetic properties)
  Vebjoern Bakken,          University of Oslo,           Norway      (DALTON; geometry optimizer, symmetry detection)
  Radovan Bast,             KTH Stockholm,                Sweden      (DALTON installation and execution frameworks)
  Pablo Baudin,             University of Valencia,       Spain       (Cholesky excitation energies)
  Linus Boman,              NTNU,                         Norway      (Cholesky decomposition and subsystems)
  Ove Christiansen,         Aarhus University,            Denmark     (CC module)
  Renzo Cimiraglia,         University of Ferrara,        Italy       (NEVPT2)
  Sonia Coriani,            University of Trieste,        Italy       (CC module, MCD in RESPONS)
  Janusz Cukras,            University of Trieste,        Italy       (MChD in RESPONS)
  Paal Dahle,               University of Oslo,           Norway      (Parallelization)
  Erik K. Dalskov,          UNI-C,                        Denmark     (SOPPA)
  Thomas Enevoldsen,        Univ. of Southern Denmark,    Denmark     (SOPPA)
  Janus J. Eriksen,         Aarhus University,            Denmark     (Polarizable embedding model, TDA)
  Berta Fernandez,          U. of Santiago de Compostela, Spain       (doublet spin, ESR in RESPONS)
  Lara Ferrighi,            Aarhus University,            Denmark     (PCM Cubic response)
  Heike Fliegl,             University of Oslo,           Norway      (CCSD(R12))
  Luca Frediani,            UiT The Arctic U. of Norway,  Norway      (PCM)
  Bin Gao,                  UiT The Arctic U. of Norway,  Norway      (Gen1Int library)
  Christof Haettig,         Ruhr-University Bochum,       Germany     (CC module)
  Kasper Hald,              Aarhus University,            Denmark     (CC module)
  Asger Halkier,            Aarhus University,            Denmark     (CC module)
  Erik D. Hedegaard,        Univ. of Southern Denmark,    Denmark     (Polarizable embedding model, QM/MM)
  Hanne Heiberg,            University of Oslo,           Norway      (geometry analysis, selected one-electron integrals)
  Trygve Helgaker,          University of Oslo,           Norway      (DALTON; ABACUS, ERI, DFT modules, London, and much more)
  Alf Christian Hennum,     University of Oslo,           Norway      (Parity violation)
  Hinne Hettema,            University of Auckland,       New Zealand (quadratic response in RESPONS; SIRIUS supersymmetry)
  Eirik Hjertenaes,         NTNU,                         Norway      (Cholesky decomposition)
  Maria Francesca Iozzi,    University of Oslo,           Norway      (RPA)
  Brano Jansik              Technical Univ. of Ostrava    Czech Rep.  (DFT cubic response)
  Hans Joergen Aa. Jensen,  Univ. of Southern Denmark,    Denmark     (DALTON; SIRIUS, RESPONS, ABACUS modules, London, and much more)
  Dan Jonsson,              UiT The Arctic U. of Norway,  Norway      (cubic response in RESPONS module)
  Poul Joergensen,          Aarhus University,            Denmark     (RESPONS, ABACUS, and CC modules)
  Maciej Kaminski,          University of Warsaw,         Poland      (CPPh in RESPONS)
  Joanna Kauczor,           Linkoeping University,        Sweden      (Complex polarization propagator (CPP) module)
  Sheela Kirpekar,          Univ. of Southern Denmark,    Denmark     (Mass-velocity & Darwin integrals)
  Wim Klopper,              KIT Karlsruhe,                Germany     (R12 code in CC, SIRIUS, and ABACUS modules)
  Stefan Knecht,            ETH Zurich,                   Switzerland (Parallel CI and MCSCF)
  Rika Kobayashi,           Australian National Univ.,    Australia   (DIIS in CC, London in MCSCF)
  Henrik Koch,              NTNU,                         Norway      (CC module, Cholesky decomposition)
  Jacob Kongsted,           Univ. of Southern Denmark,    Denmark     (Polarizable embedding model, QM/MM)
  Andrea Ligabue,           University of Modena,         Italy       (CTOCD, AOSOPPA)
  Nanna H. List             Univ. of Southern Denmark,    Denmark     (Polarizable embedding model)
  Ola B. Lutnaes,           University of Oslo,           Norway      (DFT Hessian)
  Juan I. Melo,             University of Buenos Aires,   Argentina   (LRESC, Relativistic Effects on NMR Shieldings)
  Kurt V. Mikkelsen,        University of Copenhagen,     Denmark     (MC-SCRF and QM/MM)
  Rolf H. Myhre,            NTNU,                         Norway      (Cholesky, subsystems and ECC2)
  Christian Neiss,          Univ. Erlangen-Nuernberg,     Germany     (CCSD(R12))
  Christian B. Nielsen,     University of Copenhagen,     Denmark     (QM/MM)
  Patrick Norman,           Linkoeping University,        Sweden      (Cubic response and complex response in RESPONS)
  Jeppe Olsen,              Aarhus University,            Denmark     (SIRIUS CI/density modules)
  Jogvan Magnus H. Olsen,   Univ. of Southern Denmark,    Denmark     (Polarizable embedding model, QM/MM)
  Anders Osted,             Copenhagen University,        Denmark     (QM/MM)
  Martin J. Packer,         University of Sheffield,      UK          (SOPPA)
  Filip Pawlowski,          Kazimierz Wielki University,  Poland      (CC3)
  Morten N. Pedersen,       Univ. of Southern Denmark,    Denmark     (Polarizable embedding model)
  Thomas B. Pedersen,       University of Oslo,           Norway      (Cholesky decomposition)
  Patricio F. Provasi,      University of Northeastern,   Argentina   (Analysis of coupling constants in localized orbitals)
  Zilvinas Rinkevicius,     KTH Stockholm,                Sweden      (open-shell DFT, ESR)
  Elias Rudberg,            KTH Stockholm,                Sweden      (DFT grid and basis info)
  Torgeir A. Ruden,         University of Oslo,           Norway      (Numerical derivatives in ABACUS)
  Kenneth Ruud,             UiT The Arctic U. of Norway,  Norway      (DALTON; ABACUS magnetic properties and much more)
  Pawel Salek,              KTH Stockholm,                Sweden      (DALTON; DFT code)
  Claire C. M. Samson       University of Karlsruhe       Germany     (Boys localization, r12 integrals in ERI)
  Alfredo Sanchez de Meras, University of Valencia,       Spain       (CC module, Cholesky decomposition)
  Trond Saue,               Paul Sabatier University,     France      (direct Fock matrix construction)
  Stephan P. A. Sauer,      University of Copenhagen,     Denmark     (SOPPA(CCSD), SOPPA prop., AOSOPPA, vibrational g-factors)
  Bernd Schimmelpfennig,    Forschungszentrum Karlsruhe,  Germany     (AMFI module)
  Kristian Sneskov,         Aarhus University,            Denmark     (Polarizable embedding model, QM/MM)
  Arnfinn H. Steindal,      UiT The Arctic U. of Norway,  Norway      (parallel QM/MM, Polarizable embedding model)
  Casper Steinmann,         Univ. of Southern Denmark,    Denmark     (QFIT, Polarizable embedding model)
  K. O. Sylvester-Hvid,     University of Copenhagen,     Denmark     (MC-SCRF)
  Peter R. Taylor,          VLSCI/Univ. of Melbourne,     Australia   (Symmetry handling ABACUS, integral transformation)
  Andrew M. Teale,          University of Nottingham,     England     (DFT-AC, DFT-D)
  David P. Tew,             University of Bristol,        England     (CCSD(R12))
  Olav Vahtras,             KTH Stockholm,                Sweden      (triplet response, spin-orbit, ESR, TDDFT, open-shell DFT)
  David J. Wilson,          La Trobe University,          Australia   (DFT Hessian and DFT magnetizabilities)
  Hans Agren,               KTH Stockholm,                Sweden      (SIRIUS module, RESPONS, MC-SCRF solvation model)
 --------------------------------------------------------------------------------

     Date and time (Linux)  : Tue Jun 23 22:14:19 2015
     Host name              : archer                                  

 * Work memory size             :    64000000 =  488.28 megabytes.
 + memory for in-core integrals :   100000000

 * Directories for basis set searches:
   1) /home/ignat/test/water
   2) /home/ignat/repos/dalton/build_gnu/basis


Compilation information
-----------------------

 Who compiled             | ignat
 Host                     | archer
 System                   | Linux-4.0.5-1-ARCH
 CMake generator          | Unix Makefiles
 Processor                | x86_64
 64-bit integers          | OFF
 MPI                      | On
 Fortran compiler         | /usr/bin/mpif90
 Fortran compiler version | GNU Fortran (GCC) 5.1.0
 Fortran flags            | -DVAR_GFORTRAN -DGFORTRAN=445 -ffloat-store -fcray
                          | -pointer -m64  -O0 -g -fbacktrace -fcray-pointer -
                          | Wuninitialized
 C compiler               | /usr/bin/mpicc
 C compiler version       | gcc (GCC) 5.1.0
 C flags                  | -std=c99 -DRESTRICT=restrict -DFUNDERSCORE=1 -DHAV
                          | E_NO_LSEEK64 -ffloat-store -Wall -m64 -O0 -g3
 C++ compiler             | /usr/bin/mpicxx
 C++ compiler version     | unknown
 C++ flags                | -g -Wall -fno-rtti -fno-exceptions -m64 -march=nat
                          | ive -O0 -g3
 BLAS                     | /usr/lib/libblas.so
 LAPACK                   | /usr/lib/liblapack.so
 Static linking           | OFF
 Last Git revision        | 9e6893dfe1675186f4e5a0c5ca97afe725e7638b
 Git branch               | master
 Configuration time       | 2015-06-13 18:52:22.411053

 * MPI run using 4 processes.


   Content of the .dal input file
 ----------------------------------

**DALTON INPUT                                    
.RUN RESPONSE                                     
.DIRECT                                           
.PARALLELL                                        
**WAVE FUNCTION                                   
.HF                                               
.INTERFACE                                        
**INTEGRAL                                        
.DIPLEN                                           
.SECMOM                                           
**RESPONSE                                        
.PROPAV                                           
XDIPLEN                                           
.PROPAV                                           
YDIPLEN                                           
.PROPAV                                           
ZDIPLEN                                           
*QUADRATIC                                        
.QLOP                                             
.DIPLEN                                           
**END OF DALTON INPUT                             


   Content of the .mol file
 ----------------------------

ATOMBASIS                                                                      
                                                                               
                                                                               
Atomtypes=2 Charge=0 Nosymm                                                    
Charge=8.0 Atoms=1 Basis=ano-1 4 3 1                                           
O                 0.00000   0.00000   0.00000                                  
Charge=1.0 Atoms=2 Basis=ano-1 2                                               
H                 1.43043   0.00000   1.10716                                  
H                -1.43043   0.00000   1.10716                                  


       *******************************************************************
       *********** Output from DALTON general input processing ***********
       *******************************************************************

 --------------------------------------------------------------------------------
   Overall default print level:    0
   Print level for DALTON.STAT:    1

    Parallel calculation using MPI
    AO-direct calculation (in sections where implemented)
    HERMIT 1- and 2-electron integral sections will be executed
    "Old" integral transformation used (limited to max 255 basis functions)
    Wave function sections will be executed (SIRIUS module)
    Dynamic molecular response properties section will be executed (RESPONSE module)
 --------------------------------------------------------------------------------


   ****************************************************************************
   *************** Output of molecule and basis set information ***************
   ****************************************************************************


    The two title cards from your ".mol" input:
    ------------------------------------------------------------------------
 1:                                                                         
 2:                                                                         
    ------------------------------------------------------------------------

  Atomic type no.    1
  --------------------
  Nuclear charge:   8.00000
  Number of symmetry independent centers:    1
  Number of basis sets to read;    2
  The basis set is "ano-1 4 3 1" from the basis set library.
  Basis set file used for this atomic type with Z =   8 :
     "/home/ignat/repos/dalton/build_gnu/basis/ano-1"

  Atomic type no.    2
  --------------------
  Nuclear charge:   1.00000
  Number of symmetry independent centers:    2
  Number of basis sets to read;    2
  The basis set is "ano-1 2" from the basis set library.
  Basis set file used for this atomic type with Z =   1 :
     "/home/ignat/repos/dalton/build_gnu/basis/ano-1"


                         SYMGRP: Point group information
                         -------------------------------

@    Point group: C1 


                                 Isotopic Masses
                                 ---------------

                           O          15.994915
                           H           1.007825
                           H           1.007825

                       Total mass:    18.010565 amu
                       Natural abundance:  99.730 %

 Center-of-mass coordinates (a.u.):    0.000000    0.000000    0.123908


  Atoms and basis sets
  --------------------

  Number of atom types :    2
  Total number of atoms:    3

  label    atoms   charge   prim   cont     basis
  ----------------------------------------------------------------------
  O           1    8.0000    61    18      [14s9p4d|4s3p1d]                                   
  H           2    1.0000     8     2      [8s|2s]                                            
  ----------------------------------------------------------------------
  total:      3   10.0000    77    22
  ----------------------------------------------------------------------
  Spherical harmonic basis used.

  Threshold for neglecting AO integrals:  1.00D-12


  Cartesian Coordinates (a.u.)
  ----------------------------

  Total number of coordinates:    9
  O       :     1  x   0.0000000000    2  y   0.0000000000    3  z   0.0000000000
  H       :     4  x   1.4304300000    5  y   0.0000000000    6  z   1.1071600000
  H       :     7  x  -1.4304300000    8  y   0.0000000000    9  z   1.1071600000


   Interatomic separations (in Angstrom):
   --------------------------------------

            O           H           H     
            ------      ------      ------
 O     :    0.000000
 H     :    0.957201    0.000000
 H     :    0.957201    1.513902    0.000000


  Max    interatomic separation is    1.5139 Angstrom (    2.8609 Bohr)
  between atoms    3 and    2, "H     " and "H     ".

  Min HX interatomic separation is    0.9572 Angstrom (    1.8088 Bohr)


  Bond distances (Angstrom):
  --------------------------

                  atom 1     atom 2       distance
                  ------     ------       --------
  bond distance:  H          O            0.957201
  bond distance:  H          O            0.957201


  Bond angles (degrees):
  ----------------------

                  atom 1     atom 2     atom 3         angle
                  ------     ------     ------         -----
  bond angle:     H          O          H            104.520




 Principal moments of inertia (u*A**2) and principal axes
 --------------------------------------------------------

   IA       0.614459          1.000000    0.000000    0.000000
   IB       1.154917          0.000000    0.000000    1.000000
   IC       1.769375          0.000000    1.000000    0.000000


 Rotational constants
 --------------------

@    The molecule is planar.

               A                   B                   C

         822478.2742         437589.1937         285625.6621 MHz
           27.434922           14.596404            9.527447 cm-1


@  Nuclear repulsion energy :    9.194951107924 Hartree


                     .---------------------------------------.
                     | Starting in Integral Section (HERMIT) |
                     `---------------------------------------'



 ***************************************************************************************
 ****************** Output from **INTEGRALS input processing (HERMIT) ******************
 ***************************************************************************************



    *************************************************************************
    ****************** Output from HERMIT input processing ******************
    *************************************************************************


 Default print level:        1

 * Nuclear model: Point charge

 Calculation of one-electron Hamiltonian integrals.

 The following one-electron property integrals are calculated as requested:
          - overlap integrals
          - dipole length integrals
          - second moment integrals

 Center of mass  (bohr):      0.000000000000      0.000000000000      0.123907664973
 Operator center (bohr):      0.000000000000      0.000000000000      0.000000000000
 Gauge origin    (bohr):      0.000000000000      0.000000000000      0.000000000000
 Dipole origin   (bohr):      0.000000000000      0.000000000000      0.000000000000


     ************************************************************************
     ************************** Output from HERINT **************************
     ************************************************************************

 >>>  Time used in ONEDRV     is   0.15 seconds
 >>>  Time used in QUADRUP    is   0.27 seconds
 >>>  Time used in KINENE     is   0.28 seconds
 >>>  Time used in SECMOM     is   0.27 seconds
 >>>  Time used in GABGEN     is   0.28 seconds
 >>>> Total CPU  time used in HERMIT:   1.49 seconds
 >>>> Total wall time used in HERMIT:   1.49 seconds


                        .----------------------------------.
                        | End of Integral Section (HERMIT) |
                        `----------------------------------'



                   .--------------------------------------------.
                   | Starting in Wave Function Section (SIRIUS) |
                   `--------------------------------------------'


 *** Output from Huckel module :

     Using EWMO model:          T
     Using EHT  model:          F
     Number of Huckel orbitals each symmetry:    7

 EWMO - Energy Weighted Maximum Overlap - is a Huckel type method,
        which normally is better than Extended Huckel Theory.
 Reference: Linderberg and Ohrn, Propagators in Quantum Chemistry (Wiley, 1973)

 Huckel EWMO eigenvalues for symmetry :  1
          -20.684968      -1.611697      -0.778263      -0.688371      -0.616200
           -0.232270      -0.168131

 **********************************************************************
 *SIRIUS* a direct, restricted step, second order MCSCF program       *
 **********************************************************************

 
     Date and time (Linux)  : Tue Jun 23 22:14:21 2015
     Host name              : archer                                  

 Title lines from ".mol" input file:
                                                                             
                                                                             

 Print level on unit LUPRI =   2 is   0
 Print level on unit LUW4  =   2 is   5

@    Restricted, closed shell Hartree-Fock calculation.

@    Time-dependent Hartree-Fock calculation (random phase approximation).
 Fock matrices are calculated directly and in parallel without use of integrals on disk.

 Initial molecular orbitals are obtained according to
 ".MOSTART EWMO  " input option

     Wave function specification
     ============================
@    Wave function type        >>> HF <<<
@    Number of closed shell electrons          10
@    Number of electrons in active shells       0
@    Total charge of the molecule               0

@    Spin multiplicity and 2 M_S                1         0
@    Total number of symmetries                 1 (point group: C1 )
@    Reference state symmetry                   1 (irrep name : A  )

     Orbital specifications
     ======================
@    Abelian symmetry species          All |    1
@                                          |  A  
                                       --- |  ---
@    Occupied SCF orbitals               5 |    5
@    Secondary orbitals                 17 |   17
@    Total number of orbitals           22 |   22
@    Number of basis functions          22 |   22

     Optimization information
     ========================
@    Number of configurations                 1
@    Number of orbital rotations             85
     ------------------------------------------
@    Total number of variables               86

     Maximum number of Fock   iterations      0
     Maximum number of DIIS   iterations     60
     Maximum number of QC-SCF iterations     60
     Threshold for SCF convergence     1.00D-05


 ***********************************************
 ***** DIIS acceleration of SCF iterations *****
 ***********************************************

 C1-DIIS algorithm; max error vectors =    8

 Iter      Total energy        Error norm    Delta(E)  DIIS dim.
 -----------------------------------------------------------------------------
@  1    -75.8864462763        1.48374D+00   -7.59D+01    1
      Virial theorem: -V/T =      1.997526
@      MULPOP O      -0.74; H       0.37; H       0.37; 
 -----------------------------------------------------------------------------
@  2    -76.0298087217        3.23206D-01   -1.43D-01    2
      Virial theorem: -V/T =      2.002410
@      MULPOP O      -0.68; H       0.34; H       0.34; 
 -----------------------------------------------------------------------------
@  3    -76.0358019237        8.47837D-02   -5.99D-03    3
      Virial theorem: -V/T =      1.999778
@      MULPOP O      -0.68; H       0.34; H       0.34; 
 -----------------------------------------------------------------------------
@  4    -76.0364968396        4.21539D-02   -6.95D-04    4
      Virial theorem: -V/T =      2.001194
@      MULPOP O      -0.66; H       0.33; H       0.33; 
 -----------------------------------------------------------------------------
@  5    -76.0365782749        6.52162D-03   -8.14D-05    5
      Virial theorem: -V/T =      2.000434
@      MULPOP O      -0.67; H       0.33; H       0.33; 
 -----------------------------------------------------------------------------
@  6    -76.0365858837        1.62052D-03   -7.61D-06    6
      Virial theorem: -V/T =      2.000420
@      MULPOP O      -0.67; H       0.34; H       0.34; 
 -----------------------------------------------------------------------------
@  7    -76.0365862987        2.28849D-04   -4.15D-07    7
      Virial theorem: -V/T =      2.000445
@      MULPOP O      -0.67; H       0.34; H       0.34; 
 -----------------------------------------------------------------------------
@  8    -76.0365863067        3.19997D-05   -8.07D-09    8
      Virial theorem: -V/T =      2.000444
@      MULPOP O      -0.67; H       0.34; H       0.34; 
 -----------------------------------------------------------------------------
@  9    -76.0365863069        4.79125D-06   -1.43D-10    8

@ *** DIIS converged in   9 iterations !
@     Converged SCF energy, gradient:    -76.036586306879    4.79D-06
    - total time used in SIRFCK :              0.00 seconds


 *** SCF orbital energy analysis ***

 Number of electrons :   10
 Orbital occupations :    5

 Sym       Hartree-Fock orbital energies

1 A     -20.57082761    -1.35552717    -0.72515978    -0.58965453    -0.51438127
          0.06726777     0.19936705     0.28276809     0.30520709     0.30854233
          0.41052344     0.75137306     0.93835825     1.97004726     1.97082070
          2.03471323     2.04464288     2.07836885     2.10854643     2.27493764
          3.08575297     3.64296060

    E(LUMO) :     0.06726777 au (symmetry 1)
  - E(HOMO) :    -0.51438127 au (symmetry 1)
  ------------------------------------------
    gap     :     0.58164903 au

 >>> Writing SIRIFC interface file

 >>>> CPU and wall time for SCF :       0.793       0.792


                       .-----------------------------------.
                       | >>> Final results from SIRIUS <<< |
                       `-----------------------------------'


@    Spin multiplicity:           1
@    Spatial symmetry:            1 ( irrep  A   in C1  )
@    Total charge of molecule:    0

@    Final HF energy:             -76.036586306879                 
@    Nuclear repulsion:             9.194951107924
@    Electronic energy:           -85.231537414803

@    Final gradient norm:           0.000004791249

 
     Date and time (Linux)  : Tue Jun 23 22:14:22 2015
     Host name              : archer                                  

File label for MO orbitals:  23Jun15   FOCKDIIS

 (Only coefficients >0.0100 are printed.)

 Molecular orbitals for symmetry species 1  (A  )
 ------------------------------------------------

    Orbital         1        2        3        4        5        6        7
   1 O   :1s    -1.0000   0.0232   0.0000   0.0289  -0.0000  -0.1739  -0.0000
   2 O   :1s    -0.0004  -0.7853   0.0000   0.4549  -0.0000  -1.4539  -0.0000
   3 O   :1s     0.0001   0.0782   0.0000   0.0961   0.0000  -1.3926  -0.0000
   4 O   :1s    -0.0003   0.0313   0.0000   0.0130   0.0000  -0.4518  -0.0000
   5 O   :2px    0.0000  -0.0000  -0.7141   0.0000  -0.0000  -0.0000   0.7781
   6 O   :2py    0.0000  -0.0000  -0.0000   0.0000   0.9980  -0.0000  -0.0000
   7 O   :2pz   -0.0007  -0.0956   0.0000  -0.8242   0.0000  -0.4563   0.0000
   8 O   :2px    0.0000  -0.0000   0.0365  -0.0000  -0.0000  -0.0000   0.3204
   9 O   :2py    0.0000  -0.0000  -0.0000   0.0000   0.0487  -0.0000  -0.0000
  10 O   :2pz    0.0007   0.0546   0.0000   0.0115   0.0000  -0.3844  -0.0000
  11 O   :2px    0.0000  -0.0000   0.0160  -0.0000  -0.0000  -0.0000  -0.0634
  12 O   :2py    0.0000  -0.0000  -0.0000   0.0000  -0.0220  -0.0000  -0.0000
  13 O   :2pz   -0.0008   0.0205   0.0000   0.0254   0.0000  -0.0635  -0.0000
  15 O   :3d1-  -0.0000   0.0000   0.0000  -0.0000   0.0335   0.0000  -0.0000
  16 O   :3d0    0.0000  -0.0058  -0.0000  -0.0347   0.0000  -0.0068   0.0000
  17 O   :3d1+   0.0000   0.0000  -0.0529  -0.0000   0.0000  -0.0000   0.0258
  18 O   :3d2+   0.0002  -0.0153   0.0000  -0.0079   0.0000  -0.0128  -0.0000
  19 H   :1s    -0.0003  -0.1834  -0.2759  -0.2306  -0.0000   1.1625  -1.5677
  20 H   :1s     0.0002   0.0112   0.1117   0.0144  -0.0000   0.9376  -1.2823
  21 H   :1s    -0.0003  -0.1834   0.2759  -0.2306   0.0000   1.1625   1.5677
  22 H   :1s     0.0002   0.0112  -0.1117   0.0144  -0.0000   0.9376   1.2823

    Orbital         8        9       10       11       12       13       14
   1 O   :1s    -0.2139   0.1097   0.0000   0.0000   0.1569  -0.0000   0.0000
   2 O   :1s    -1.7441   0.9149   0.0000   0.0000   1.3858  -0.0000   0.0000
   3 O   :1s    -3.0951   1.1760  -0.0000   0.0000   0.8992  -0.0000  -0.0000
   4 O   :1s    -1.3323   0.4375  -0.0000   0.0000   0.3096  -0.0000  -0.0000
   5 O   :2px   -0.0000   0.0000   0.0000   1.3472  -0.0000  -0.0278   0.0000
   6 O   :2py    0.0000  -0.0000  -0.0343   0.0000  -0.0000  -0.0000  -0.0470
   7 O   :2pz   -0.4287   0.2238  -0.0000   0.0000   0.9047  -0.0000  -0.0000
   8 O   :2px   -0.0000  -0.0000   0.0000   1.9871  -0.0000   1.3420   0.0000
   9 O   :2py    0.0000  -0.0000   0.9016   0.0000  -0.0000  -0.0000   0.4206
  10 O   :2pz   -0.5821   1.2365  -0.0000   0.0000   0.4394  -0.0000  -0.0000
  11 O   :2px   -0.0000  -0.0000   0.0000   0.5497  -0.0000   0.3521   0.0000
  12 O   :2py    0.0000  -0.0000   0.4311   0.0000  -0.0000  -0.0000  -0.8809
  13 O   :2pz   -0.1081   0.4990  -0.0000   0.0000   0.0307  -0.0000  -0.0000
  15 O   :3d1-  -0.0000   0.0000  -0.0052   0.0000   0.0000   0.0000   0.2121
  16 O   :3d0   -0.0032  -0.0033  -0.0000  -0.0000   0.0195  -0.0000  -0.0000
  17 O   :3d1+  -0.0000   0.0000  -0.0000  -0.0172  -0.0000   0.1293  -0.0000
  18 O   :3d2+  -0.0125   0.0014  -0.0000   0.0000  -0.1222   0.0000  -0.0000
  19 H   :1s     1.4254  -0.7456   0.0000  -2.3719  -1.2763  -1.2532   0.0000
  20 H   :1s     1.1629  -0.5317   0.0000  -1.4834   0.0461  -1.9624   0.0000
  21 H   :1s     1.4254  -0.7456  -0.0000   2.3719  -1.2763   1.2532  -0.0000
  22 H   :1s     1.1629  -0.5317   0.0000   1.4834   0.0461   1.9624   0.0000

    Orbital        15
   2 O   :1s    -0.0704
   3 O   :1s    -0.2560
   4 O   :1s     0.1429
   7 O   :2pz   -0.0684
  10 O   :2pz    0.3265
  13 O   :2pz   -0.7630
  16 O   :3d0    0.4462
  18 O   :3d2+  -0.2304
  19 H   :1s     0.0476
  20 H   :1s     0.0728
  21 H   :1s     0.0476
  22 H   :1s     0.0728



 >>>> Total CPU  time used in SIRIUS :      0.81 seconds
 >>>> Total wall time used in SIRIUS :      0.81 seconds

 
     Date and time (Linux)  : Tue Jun 23 22:14:22 2015
     Host name              : archer                                  


                     .---------------------------------------.
                     | End of Wave Function Section (SIRIUS) |
                     `---------------------------------------'



                 .------------------------------------------------.
                 | Starting in Dynamic Property Section (RESPONS) |
                 `------------------------------------------------'


 ------------------------------------------------------------------------------
  RESPONSE  -  an MCSCF, MC-srDFT, DFT, and SOPPA response property program
 ------------------------------------------------------------------------------


 <<<<<<<<<< OUTPUT FROM RESPONSE INPUT PROCESSING >>>>>>>>>>




 CHANGES OF DEFAULTS FOR RSPINP:
 -------------------------------


 AO-direct Fock matrix calculations.

 Default : Using Fock type decoupling of the two-electron density matrix :
    Add DV*(FC+FV) instead of DV*FC to E[2] approximate orbital diagonal



 Quadratic Response calculation
 ------------------------------

 First hyperpolarizability calculation : HYPCAL= T

 Spin of operator A , ISPINA=    0
 Spin of operator B , ISPINB=    0
 Spin of operator C , ISPINC=    0

  1 B-frequencies  0.000000D+00
  1 C-frequencies  0.000000D+00

 Print level                                    : IPRHYP =   2
 Maximum number of iterations in lin.rsp. solver: MAXITL =  60
 Threshold for convergence of linear resp. eq.s : THCLR  = 1.000D-03
 Maximum iterations in optimal orbital algorithm: MAXITO =   5
 Direct one-index transformation                : DIROIT = T

    3 A OPERATORS OF SYMMETRY NO:    1 AND LABELS:

          XDIPLEN 
          YDIPLEN 
          ZDIPLEN 

    3 B OPERATORS OF SYMMETRY NO:    1 AND LABELS:

          XDIPLEN 
          YDIPLEN 
          ZDIPLEN 

    3 C OPERATORS OF SYMMETRY NO:    1 AND LABELS:

          XDIPLEN 
          YDIPLEN 
          ZDIPLEN 


   SCF energy         :      -76.036586306878860
 -- inactive part     :      -85.231537414802574
 -- nuclear repulsion :        9.194951107923721


                     ***************************************
                     *** RHF response calculation (TDHF) ***
                     ***************************************



             Calculation of electronic one-electron expectation values
            ----------------------------------------------------------

 (Note that to get e.g. a dipole moment you must multiply the
  electronic number by -1 and add the nuclear contribution.)

 *** Individual non-zero orbital contributions
 *** to the expectation value for property XDIPLEN  :

     Inactive       1  1 in sym 1 :    -0.00000000
     Inactive       2  2 in sym 1 :     0.00000000
     Inactive       3  3 in sym 1 :     0.00000000
     Inactive       4  4 in sym 1 :     0.00000000
     Inactive       5  5 in sym 1 :    -0.00000000

     XDIPLEN  inactive part: 7.44252169D-15
     XDIPLEN  active part  : 0.00000000D+00
     XDIPLEN  total        : 7.44252169D-15

 *** Individual non-zero orbital contributions
 *** to the expectation value for property YDIPLEN  :

     Inactive       1  1 in sym 1 :    -0.00000000
     Inactive       2  2 in sym 1 :     0.00000000
     Inactive       3  3 in sym 1 :    -0.00000000
     Inactive       4  4 in sym 1 :     0.00000000
     Inactive       5  5 in sym 1 :    -0.00000000

     YDIPLEN  inactive part: 7.53885075D-16
     YDIPLEN  active part  : 0.00000000D+00
     YDIPLEN  total        : 7.53885075D-16

 *** Individual non-zero orbital contributions
 *** to the expectation value for property ZDIPLEN  :

     Inactive       1  1 in sym 1 :     0.00049027
     Inactive       2  2 in sym 1 :     0.64776114
     Inactive       3  3 in sym 1 :     0.83917250
     Inactive       4  4 in sym 1 :    -0.20936129
     Inactive       5  5 in sym 1 :     0.08212857

     ZDIPLEN  inactive part:     1.36019118
     ZDIPLEN  active part  :     0.00000000
     ZDIPLEN  total        :     1.36019118


 Linear response calculations for quadratic response
 - singlet property operator of symmetry    1  ( A  )


 Perturbation symmetry.     KSYMOP:       1
 Perturbation spin symmetry.TRPLET:       F
 Orbital variables.         KZWOPT:      85
 Configuration variables.   KZCONF:       0
 Total number of variables. KZVAR :      85


 QRLRVE -- linear response calculation for symmetry  1  ( A  )
 QRLRVE -- operator label : XDIPLEN 
 QRLRVE -- operator spin  :   0
 QRLRVE -- frequencies :  0.000000



 <<<  SOLVING SETS OF LINEAR EQUATIONS FOR LINEAR RESPONSE PROPERTIES >>>

 Operator symmetry = 1  ( A  ); triplet =   F


 *** THE REQUESTED    1 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00D-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   10)
 RSP solution vector no.    1; norm of residual   5.09D-04

 *** RSPCTL MICROITERATIONS CONVERGED

@ QRLRVE: SINGLET SOLUTION FOR SYMMETRY    1  ( A  )    LABEL   XDIPLEN     FREQUENCY   0.000000D+00

@ QRLRVE:  << XDIPLEN  ; XDIPLEN  >> (   0.00000):     6.88797139440    

@ QRLRVE:  << YDIPLEN  ; XDIPLEN  >> (   0.00000):   -2.457085259705E-15

@ QRLRVE:  << ZDIPLEN  ; XDIPLEN  >> (   0.00000):    4.993539086389E-14


 QRLRVE -- linear response calculation for symmetry  1  ( A  )
 QRLRVE -- operator label : YDIPLEN 
 QRLRVE -- operator spin  :   0
 QRLRVE -- frequencies :  0.000000



 <<<  SOLVING SETS OF LINEAR EQUATIONS FOR LINEAR RESPONSE PROPERTIES >>>

 Operator symmetry = 1  ( A  ); triplet =   F


 *** THE REQUESTED    1 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00D-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   12)
 RSP solution vector no.    1; norm of residual   5.04D-05

 *** RSPCTL MICROITERATIONS CONVERGED

@ QRLRVE: SINGLET SOLUTION FOR SYMMETRY    1  ( A  )    LABEL   YDIPLEN     FREQUENCY   0.000000D+00

@ QRLRVE:  << XDIPLEN  ; YDIPLEN  >> (   0.00000):   -2.463163106656E-15

@ QRLRVE:  << YDIPLEN  ; YDIPLEN  >> (   0.00000):     5.18329883914    

@ QRLRVE:  << ZDIPLEN  ; YDIPLEN  >> (   0.00000):    1.063418866696E-15


 QRLRVE -- linear response calculation for symmetry  1  ( A  )
 QRLRVE -- operator label : ZDIPLEN 
 QRLRVE -- operator spin  :   0
 QRLRVE -- frequencies :  0.000000



 <<<  SOLVING SETS OF LINEAR EQUATIONS FOR LINEAR RESPONSE PROPERTIES >>>

 Operator symmetry = 1  ( A  ); triplet =   F


 *** THE REQUESTED    1 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00D-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   12)
 RSP solution vector no.    1; norm of residual   2.92D-04

 *** RSPCTL MICROITERATIONS CONVERGED

@ QRLRVE: SINGLET SOLUTION FOR SYMMETRY    1  ( A  )    LABEL   ZDIPLEN     FREQUENCY   0.000000D+00

@ QRLRVE:  << XDIPLEN  ; ZDIPLEN  >> (   0.00000):    4.355582821759E-14

@ QRLRVE:  << YDIPLEN  ; ZDIPLEN  >> (   0.00000):    1.521124169556E-15

@ QRLRVE:  << ZDIPLEN  ; ZDIPLEN  >> (   0.00000):     5.94876530901    
 
 ======================================================================
 >>>>>>>>    L I N E A R   R E S P O N S E   F U N C T I O N S   <<<<<<<<
 ======================================================================

 The -<<A;B>>(omega_b) functions from vectors generated
 in a *QUADRA calculation of <<A;B,C>>(omega_b,omega_c)

 Note: the accuracy of off-diagonal elements will be linear
 in the convergence threshold THCLR =  1.00D-03


 Perturbation symmetry.     KSYMOP:       1
 Perturbation spin symmetry.TRPLET:       F
 Orbital variables.         KZWOPT:      85
 Configuration variables.   KZCONF:       0
 Total number of variables. KZVAR :      85

@ Singlet linear response function in a.u.

@ A operator, symmetry, frequency:   XDIPLEN    1 -0.000000
@ B operator, symmetry, frequency:   XDIPLEN    1  0.000000

@ Value of linear response -<<A;B>>(omega):       6.887971394397

@ Singlet linear response function in a.u.

@ A operator, symmetry, frequency:   YDIPLEN    1 -0.000000
@ B operator, symmetry, frequency:   XDIPLEN    1  0.000000

@ Value of linear response -<<A;B>>(omega):      -0.000000000000

@ Singlet linear response function in a.u.

@ A operator, symmetry, frequency:   ZDIPLEN    1 -0.000000
@ B operator, symmetry, frequency:   XDIPLEN    1  0.000000

@ Value of linear response -<<A;B>>(omega):       0.000000000000

@ Singlet linear response function in a.u.

@ A operator, symmetry, frequency:   XDIPLEN    1 -0.000000
@ B operator, symmetry, frequency:   YDIPLEN    1  0.000000

@ Value of linear response -<<A;B>>(omega):      -0.000000000000

@ Singlet linear response function in a.u.

@ A operator, symmetry, frequency:   YDIPLEN    1 -0.000000
@ B operator, symmetry, frequency:   YDIPLEN    1  0.000000

@ Value of linear response -<<A;B>>(omega):       5.183298839138

@ Singlet linear response function in a.u.

@ A operator, symmetry, frequency:   ZDIPLEN    1 -0.000000
@ B operator, symmetry, frequency:   YDIPLEN    1  0.000000

@ Value of linear response -<<A;B>>(omega):       0.000000000000

@ Singlet linear response function in a.u.

@ A operator, symmetry, frequency:   XDIPLEN    1 -0.000000
@ B operator, symmetry, frequency:   ZDIPLEN    1  0.000000

@ Value of linear response -<<A;B>>(omega):       0.000000000000

@ Singlet linear response function in a.u.

@ A operator, symmetry, frequency:   YDIPLEN    1 -0.000000
@ B operator, symmetry, frequency:   ZDIPLEN    1  0.000000

@ Value of linear response -<<A;B>>(omega):       0.000000000000

@ Singlet linear response function in a.u.

@ A operator, symmetry, frequency:   ZDIPLEN    1 -0.000000
@ B operator, symmetry, frequency:   ZDIPLEN    1  0.000000

@ Value of linear response -<<A;B>>(omega):       5.948765309008


  Results from quadratic response calculation
 --------------------------------------------



 CRLRV3 -- linear response calc for sym:  1
 CRLRV3 -- operator label1: XDIPLEN 
 CRLRV3 -- operator label2: XDIPLEN 
 CRLRV3 -- freqr1 :   0.000000D+00
 CRLRV3 -- freqr2 :   0.000000D+00



 <<<  SOLVING SETS OF LINEAR EQUATIONS FOR LINEAR RESPONSE PROPERTIES >>>

 Operator symmetry = 1  ( A  ); triplet =   F


 *** THE REQUESTED    1 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00D-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   10)
 RSP solution vector no.    1; norm of residual   9.37D-04

 *** RSPCTL MICROITERATIONS CONVERGED
    XDIPLEN     XDIPLEN      freq1     freq2                Norm
 ---------------------------------------------------------------
                          0.000000  0.000000          9.56741326
@ B-freq = 0.000000  C-freq = 0.000000     beta(X;X,X) =     -0.00000000
@ B-freq = 0.000000  C-freq = 0.000000     beta(Y;X,X) =     -0.00000000
@ B-freq = 0.000000  C-freq = 0.000000     beta(Z;X,X) =    -12.10518652
@ B-freq = 0.000000  C-freq = 0.000000     beta(X;Y,X) = beta(Y,X,X)


 CRLRV3 -- linear response calc for sym:  1
 CRLRV3 -- operator label1: YDIPLEN 
 CRLRV3 -- operator label2: XDIPLEN 
 CRLRV3 -- freqr1 :   0.000000D+00
 CRLRV3 -- freqr2 :   0.000000D+00



 <<<  SOLVING SETS OF LINEAR EQUATIONS FOR LINEAR RESPONSE PROPERTIES >>>

 Operator symmetry = 1  ( A  ); triplet =   F


 *** THE REQUESTED    1 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00D-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:    8)
 RSP solution vector no.    1; norm of residual   3.95D-04

 *** RSPCTL MICROITERATIONS CONVERGED
    YDIPLEN     XDIPLEN      freq1     freq2                Norm
 ---------------------------------------------------------------
                          0.000000  0.000000          7.26114949
@ B-freq = 0.000000  C-freq = 0.000000     beta(X;Y,X) =     -0.00000000
@ B-freq = 0.000000  C-freq = 0.000000     beta(Y;Y,X) =      0.00000000
@ B-freq = 0.000000  C-freq = 0.000000     beta(Z;Y,X) =      0.00000000
@ B-freq = 0.000000  C-freq = 0.000000     beta(X;Z,X) = beta(Z,X,X)
@ B-freq = 0.000000  C-freq = 0.000000     beta(Y;Z,X) = beta(Z,Y,X)


 CRLRV3 -- linear response calc for sym:  1
 CRLRV3 -- operator label1: ZDIPLEN 
 CRLRV3 -- operator label2: XDIPLEN 
 CRLRV3 -- freqr1 :   0.000000D+00
 CRLRV3 -- freqr2 :   0.000000D+00



 <<<  SOLVING SETS OF LINEAR EQUATIONS FOR LINEAR RESPONSE PROPERTIES >>>

 Operator symmetry = 1  ( A  ); triplet =   F


 *** THE REQUESTED    1 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00D-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:    8)
 RSP solution vector no.    1; norm of residual   5.90D-04

 *** RSPCTL MICROITERATIONS CONVERGED
    ZDIPLEN     XDIPLEN      freq1     freq2                Norm
 ---------------------------------------------------------------
                          0.000000  0.000000         10.50779485
@ B-freq = 0.000000  C-freq = 0.000000     beta(X;Z,X) =    -12.09975413
@ B-freq = 0.000000  C-freq = 0.000000     beta(Y;Z,X) =      0.00000000
@ B-freq = 0.000000  C-freq = 0.000000     beta(Z;Z,X) =      0.00000000
@ B-freq = 0.000000  C-freq = 0.000000     beta(X;X,Y) = beta(Y,X,X)
@ B-freq = 0.000000  C-freq = 0.000000     beta(Y;X,Y) = beta(Y,Y,X)
@ B-freq = 0.000000  C-freq = 0.000000     beta(Z;X,Y) = beta(Z,Y,X)
@ B-freq = 0.000000  C-freq = 0.000000     beta(X;Y,Y) = beta(Y,Y,X)


 CRLRV3 -- linear response calc for sym:  1
 CRLRV3 -- operator label1: YDIPLEN 
 CRLRV3 -- operator label2: YDIPLEN 
 CRLRV3 -- freqr1 :   0.000000D+00
 CRLRV3 -- freqr2 :   0.000000D+00



 <<<  SOLVING SETS OF LINEAR EQUATIONS FOR LINEAR RESPONSE PROPERTIES >>>

 Operator symmetry = 1  ( A  ); triplet =   F


 *** THE REQUESTED    1 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00D-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   12)
 RSP solution vector no.    1; norm of residual   2.33D-04

 *** RSPCTL MICROITERATIONS CONVERGED
    YDIPLEN     YDIPLEN      freq1     freq2                Norm
 ---------------------------------------------------------------
                          0.000000  0.000000         13.91330727
@ B-freq = 0.000000  C-freq = 0.000000     beta(X;Y,Y) =      0.00000000
@ B-freq = 0.000000  C-freq = 0.000000     beta(Y;Y,Y) =     -0.00000000
@ B-freq = 0.000000  C-freq = 0.000000     beta(Z;Y,Y) =      2.18839797
@ B-freq = 0.000000  C-freq = 0.000000     beta(X;Z,Y) = beta(Z,Y,X)
@ B-freq = 0.000000  C-freq = 0.000000     beta(Y;Z,Y) = beta(Z,Y,Y)


 CRLRV3 -- linear response calc for sym:  1
 CRLRV3 -- operator label1: ZDIPLEN 
 CRLRV3 -- operator label2: YDIPLEN 
 CRLRV3 -- freqr1 :   0.000000D+00
 CRLRV3 -- freqr2 :   0.000000D+00



 <<<  SOLVING SETS OF LINEAR EQUATIONS FOR LINEAR RESPONSE PROPERTIES >>>

 Operator symmetry = 1  ( A  ); triplet =   F


 *** THE REQUESTED    1 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00D-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   10)
 RSP solution vector no.    1; norm of residual   2.49D-04

 *** RSPCTL MICROITERATIONS CONVERGED
    ZDIPLEN     YDIPLEN      freq1     freq2                Norm
 ---------------------------------------------------------------
                          0.000000  0.000000          9.08114434
@ B-freq = 0.000000  C-freq = 0.000000     beta(X;Z,Y) =      0.00000000
@ B-freq = 0.000000  C-freq = 0.000000     beta(Y;Z,Y) =      2.18809724
@ B-freq = 0.000000  C-freq = 0.000000     beta(Z;Z,Y) =     -0.00000000
@ B-freq = 0.000000  C-freq = 0.000000     beta(X;X,Z) = beta(Z,X,X)
@ B-freq = 0.000000  C-freq = 0.000000     beta(Y;X,Z) = beta(Z,Y,X)
@ B-freq = 0.000000  C-freq = 0.000000     beta(Z;X,Z) = beta(Z,Z,X)
@ B-freq = 0.000000  C-freq = 0.000000     beta(X;Y,Z) = beta(Z,Y,X)
@ B-freq = 0.000000  C-freq = 0.000000     beta(Y;Y,Z) = beta(Z,Y,Y)
@ B-freq = 0.000000  C-freq = 0.000000     beta(Z;Y,Z) = beta(Z,Z,Y)
@ B-freq = 0.000000  C-freq = 0.000000     beta(X;Z,Z) = beta(Z,Z,X)
@ B-freq = 0.000000  C-freq = 0.000000     beta(Y;Z,Z) = beta(Z,Z,Y)


 CRLRV3 -- linear response calc for sym:  1
 CRLRV3 -- operator label1: ZDIPLEN 
 CRLRV3 -- operator label2: ZDIPLEN 
 CRLRV3 -- freqr1 :   0.000000D+00
 CRLRV3 -- freqr2 :   0.000000D+00



 <<<  SOLVING SETS OF LINEAR EQUATIONS FOR LINEAR RESPONSE PROPERTIES >>>

 Operator symmetry = 1  ( A  ); triplet =   F


 *** THE REQUESTED    1 SOLUTION VECTORS CONVERGED

 Convergence of RSP solution vectors, threshold = 1.00D-03
 ---------------------------------------------------------------
 (dimension of paired reduced space:   10)
 RSP solution vector no.    1; norm of residual   9.42D-04

 *** RSPCTL MICROITERATIONS CONVERGED
    ZDIPLEN     ZDIPLEN      freq1     freq2                Norm
 ---------------------------------------------------------------
                          0.000000  0.000000         10.30701610
@ B-freq = 0.000000  C-freq = 0.000000     beta(X;Z,Z) =     -0.00000000
@ B-freq = 0.000000  C-freq = 0.000000     beta(Y;Z,Z) =     -0.00000000
@ B-freq = 0.000000  C-freq = 0.000000     beta(Z;Z,Z) =     -2.64089412

 >>>> Total CPU  time used in RESPONSE:   1.02 seconds
 >>>> Total wall time used in RESPONSE:   1.02 seconds


                   .-------------------------------------------.
                   | End of Dynamic Property Section (RESPONS) |
                   `-------------------------------------------'

 >>>> Total CPU  time used in DALTON:   3.34 seconds
 >>>> Total wall time used in DALTON:   3.35 seconds

 
     Date and time (Linux)  : Tue Jun 23 22:14:23 2015
     Host name              : archer                                  
"""

@attr(speed = 'fast' )
class UtilzTestCase( unittest.TestCase ):

    def setUp(self):
        pass

    def test_reflect_point_by_three_points(self):
        """Plane is just xy plane, should give [1, 1, -1] """
        p = np.array( [1, 1, 1] )
        p1 = np.array( [0,0,0] )
        p2 = np.array( [1,0,0] )
        p3 = np.array( [0,1,0] )
        np.testing.assert_allclose( [1, 1, -1], 
                utilz.reflect_point_by_three_points( p, p1, p2, p3 ))

        p = np.array( [-5, 0, 3] )
        p1 = np.array( [0,0,0] )
        p2 = np.array( [1,0,0] )
        p3 = np.array( [0,1,0] )
        np.testing.assert_allclose( [-5, 0, -3], 
                utilz.reflect_point_by_three_points( p, p1, p2, p3 ))

    def test_splitter(self):
        a = ['A', 'B', 'B', 'C', 'C', 'C' ]
        u = utilz.splitter( a )
        assert len(u) == 3
        assert len(u[0]) == 1
        assert len(u[1]) == 2
        assert len(u[2]) == 3

    def test_unique(self):
        a = ['A', 'B', 'C', 'C', 'C' ]
        u = utilz.unique( a )
        assert len( u ) == 3

        a = [ {'A': 1}, {'A' : 1}, {'A' : 44}, {'A' : 'dogg'}, {'A': 5} ]
        u = utilz.unique( a, lambda x: x['A'] )
        assert len( u ) == 4


    def test_center_and_xz_1(self):
        p1 = np.array( [5, 5, 5] )
        p2 = np.array( [5, 5, 6] )
        p3 = np.array( [6, 5, 5] )
        t_v, r1, r2, r3 = utilz.center_and_xz( p1, p2, p3 )
        np.testing.assert_allclose( t_v, np.full( (3,), -5.0 ) , atol =1e-10 )
        np.testing.assert_allclose( r1, 0.0 , atol =1e-10 )
        np.testing.assert_allclose( r2, 0.0 , atol =1e-10 )
        np.testing.assert_allclose( r3, 0.0 , atol =1e-10 )


    def test_get_euler(self):
        v1, v2 = map(lambda x: np.array(x), [[0,0,1], [1,0,0]] )
        r1, r2, r3 = utilz.get_euler( v1, v2 )
        np.testing.assert_allclose( [r1,r2,r3], [0,0,0] )

    def test_rotate_point_by_two_points(self):
        p = np.array( [ 1,1,1] )
        p1 = np.array( [0, 0, 0] )
        p2 = np.array( [0, 0, 1] )
        theta = np.pi/2
        p_out = utilz.rotate_point_by_two_points(p, p1, p2, theta)
        np.testing.assert_allclose( p_out, np.array( [-1,1,1] ), atol=1e-10 )

    def test_rotate_point_by_two_points_2(self):
        p = np.array( [ -1, -1, 0] )
        p1 = np.array( [ 0, 0,  1] )
        p2 = np.array( [ 0, 0, -1] )
        theta = np.pi * 3/2
        p_out = utilz.rotate_point_by_two_points(p, p1, p2, theta)
        np.testing.assert_allclose( p_out, np.array( [ 1, -1, 0] ), atol=1e-10 )

    def test_rotate_point_by_two_points_3(self):
        p = np.array( [ 0, 3, 1] )
        p1 = np.array( [ 3, 3,  0] )
        p2 = np.array( [ 3, 3,  1] )
        theta = np.pi * 3/2
        p_out = utilz.rotate_point_by_two_points(p, p1, p2, theta)
        np.testing.assert_allclose( p_out, np.array( [ 3, 6, 1] ), atol=1e-10 )

    def test_rotate_point_around_cross(self):
        p1 = np.array( [ 1, 0, 0] )
        p2 = np.array( [ 0, 0,  0] )
        p3 = np.array( [ 0, 1,  0] )
        theta = np.pi/2
        p_out = utilz.rotate_point_around_cross(p1, p2, p3, theta)
        np.testing.assert_allclose( p_out, np.array( [ 0, 1, 0] ), atol=1e-10 )

    def test_get_t_and_rho(self):
        w = Water.get_standard() 

        t, r1, r2, r3 = utilz.get_t_and_rho( w.o.r, (w.h1.r - w.h2.r)/2 + w.h2.r, w.h1.r,
                plane = 'xy' )

    def test_converged(self):
        ret = utilz.converged( HF_FILE )
        assert ret == True

    def test_dipole_iso(self):
        d = np.array( [ -1, 2, -2 ] )
        d_iso = utilz.dipole_iso( d )
        np.testing.assert_allclose( d_iso, 3, atol = 1e-7 )

    def test_beta_par(self):
        w = Water.get_standard()
        w.attach_properties()
        d = w.p.d
        b = utilz.ut2s( w.p.b )

        b_para1 = utilz.b_para( b, d ) 
        b_para2 = 1.0/5.0*( (np.einsum( 'ijj', b ) + np.einsum( 'jij', b ) +  np.einsum( 'jji', b )).sum() )
        np.testing.assert_allclose( b_para1, b_para2 )


    def test_get_rotation(self):
        p1 = np.array( [ 0, 0, 0] )
        p2 = np.array( [ 0, 0, 1] )
        p3 = np.array( [ 0, 1, 0] )
        R = utilz.get_rotation( p1, p2, p3 )
        R = np.einsum( 'ij->ji', R )
        Rz = utilz.Rz( np.pi/2.0 )
        np.testing.assert_allclose( Rz, R, atol =1e-7 )


    def test_center_of_nuclei_charge(self):
        _type = np.float
        p_n = np.array( [ [0, 0, -1], [0, 0, 1] ], dtype = _type )
        q_n = np.array( [ 1, 1 ], dtype = _type )
        coc = utilz.center_of_nuclei_charge( p_n, q_n )
        np.testing.assert_allclose( coc, np.zeros( (3,) ) )

    def test_monopole_moment(self):
        _type = np.float
        q_e = np.array( [ -1, -0.5 -0.5 ], dtype = _type )
        m = utilz.electric_monopole_moment( q_e )
        np.testing.assert_allclose( -2, m )

    def test_dipole_moment_inv(self):
        """Same dipole moment should be obtained after translating 
        coordinate frame """
        _type = np.float
        p_n = np.random.random( (5, 3,) )
        q_n = np.random.random( (5, ) )
        p_e = np.random.random( (5, 3,) )
        q_e = np.random.random( (5, ) )

        m1 = utilz.electric_dipole_moment( p_n, q_n, p_e, q_e )
        t_vector = np.random.random( (3,) ) 

        p_n += t_vector
        p_e += t_vector

        m2 = utilz.electric_dipole_moment( p_n, q_n, p_e, q_e )
        np.testing.assert_allclose( m1, m2 )

    def test_alpha_aniso2(self):
        test = np.random.random( 6 )
        a1 = utilz.alpha_aniso( test )
        a2= utilz.alpha_aniso2( test )
        np.testing.assert_allclose( a1, a2 , atol =1e-7 )



if __name__ == '__main__':
    unittest.main()
