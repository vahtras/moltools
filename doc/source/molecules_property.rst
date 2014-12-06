.. _Property:

Property
=====================


.. code-block:: python

   Class Property( list ):
      __init__(self):
         pass

**Class property** is attached to an **class Atom**.

Has support for the following properties:
   - Charge
   - Dipole moment
   - Quadrupole moment
   - Upper triangular Alpha
   - Upper triangular Beta

Most typical usage:

   1) Read in 10 water molecules from a .mol file (DALTON molecule format),
written in atomic units.

   2) Import a fix water property template given:
       - Model
       - Method
       - Basis
       - LoProp on/off
       - Frequency

   3) Rotate the properties for each water molecule to the current coordinate frame.

in code:

.. code-block:: python

   >>> from molecules import Property, Atom, Water
   >>> from template import Template
   >>> molfile = "10_tip3p.mol"
   >>> waters = Water.read_waters( mol = molfile,
   ...          in_AA = False, out_AA = False, N_waters = 10 )
   >>> for wat in waters:
   ...     t1, t2, t3 = wat.get_euler()
   ...     kw_dict = Template().get( *("TIP3P", "HF", "ANOPVDZ",
   ...         dist = True, "0,0")
   ...     for at in wat:
   ...         Property.add_prop_from_template( at, kw_dict )
   ...         at.Property.transform_ut_properties( t1, t2, t3 )
