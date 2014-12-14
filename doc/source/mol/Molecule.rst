.. _Molecule:

Molecule
========================

Calc distance between to atoms::
   >>> H1 = Molecule( **{ 'element' : H } )
   >>> H2 = Molecule( **{ 'element' : H, 'z' : 1.0 } )
   >>> print H1.dist_to_atom( H2 )
   1.0
