molecules.py/Property
=====================

Purpose:



.. code-block:: python

   Class Property( list ):
      __init__(self):
         pass


**Class property** which is attached to an atom.
Has support for multi-pole moments up to quadrupoles,
upper-triangular alpha
and upper-triangular beta


Most typical usage:

Property.add_prop_from_template( class Atom a, class Template t):

Will generate a property from Template t, which is found in template.py/Template, and make it the Property of atom A, documented in molecules.py/Atom.
