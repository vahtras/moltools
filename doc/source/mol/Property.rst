.. _Property:

Property
========================

Calc distance between to atoms::
   >>> H1 = Property( **{ 'element' : H } )
   >>> H2 = Property( **{ 'element' : H, 'z' : 1.0 } )
   >>> print H1.dist_to_atom( H2 )
   1.0
