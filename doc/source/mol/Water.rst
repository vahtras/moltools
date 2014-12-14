.. _Water:

Water
========================

Calc distance between to atoms::
   >>> H1 = Water( **{ 'element' : H } )
   >>> H2 = Water( **{ 'element' : H, 'z' : 1.0 } )
   >>> print H1.dist_to_atom( H2 )
   1.0
