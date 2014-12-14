.. _Rotator:

Rotator
========================

Calc distance between to atoms::
   >>> H1 = Rotator( **{ 'element' : H } )
   >>> H2 = Rotator( **{ 'element' : H, 'z' : 1.0 } )
   >>> print H1.dist_to_atom( H2 )
   1.0
