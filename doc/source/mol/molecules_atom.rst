.. _Atom:

Atom
=====================


.. code-block:: python

   Class Atom( object ):
      __init__(self, *args, **kwargs):
         pass


**Class Atom** is an atom representation.

Typical usage, create an oxygen atom in origo using atomic units:

.. code-block:: python

   >>> o = Atom( {'x' : 0.0, "y" :0.0, 'z' : 0.0, "element" : "O", AA = False}

