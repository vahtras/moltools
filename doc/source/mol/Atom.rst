.. _Atom:

Atom
========================

Calc distance between to atoms::
   >>> H1 = Atom( **{ 'element' : H } )
   >>> H2 = Atom( **{ 'element' : H, 'z' : 1.0 } )
   >>> print H1.dist_to_atom( H2 )
   1.0
