.. _Molecule:

Molecule
=====================


.. code-block:: python

   Class Molecule( list ):
      __init__(self, *args, **kwargs):
         pass

**Class Molecule** is a molecule representation, consisting of **class Atom** members in it's list.

Serves as a general molecule for specific molecules to enherit its generic methods.

.. code-block:: python

   >>> o = Atom( {'x' : 0.0, "y" :0.0, 'z' : 0.0, "element" : "O", AA = False}


