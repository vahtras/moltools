molecules.py
=====================

Purpose:

   Contains various classes that go together in order to create input and read output from DALTON calculations.

**Defined classes**

.. code-block:: python 

   Class Property( list ):
      __init__(self):
         pass

Used to model _properties which are attached to a class Atom object.


.. code-block:: python 

   Class Atom( object ):
      __init__(self, *args, **kwargs)

.. code-block:: python 

   Class Molecule( list ):
      __init__(self, *args, **kwargs)

.. code-block:: python 

   Class Water( Molecule ):
      __init__(self)


