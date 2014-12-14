.. _Cluster:

Cluster
========================

Calc distance between to atoms::
   >>> H1 = Cluster( **{ 'element' : H } )
   >>> H2 = Cluster( **{ 'element' : H, 'z' : 1.0 } )
   >>> print H1.dist_to_atom( H2 )
   1.0
