molecules.py
==============

Class Property
--------------------

.. code-block:: python 

   Class Property( dict ):
      __init__(self):
         pass

Used to model :ref:`Property`.


Class Atom
--------------------

.. code-block:: python 

   Class Atom( object ):
      __init__(self, *args, **kwargs)

Used to model :ref:`Atom`.

Class Molecule
--------------------

.. code-block:: python 

   Class Molecule( list ):
      __init__(self, *args, **kwargs)

Used to model :ref:`Molecule`.

Class Water
--------------------

.. code-block:: python 

   Class Water( Molecule ):
      __init__(self)

Used to model :ref:`Water`.


Class Rotator
--------------------

.. code-block:: python 

   Class Rotator( object ):
      __init__(self)
Used to model :ref:`Rotator`.


Class Cluster
--------------------

.. code-block:: python 

   Class Cluster( list ):
      __init__(self)
Used to model :ref:`Cluster`.

