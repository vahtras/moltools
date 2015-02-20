import numpy as np
import molecules

from matplotlib import pyplot as plt

a0 = 0.52917721092

class Cellnp( np.ndarray ):

    def __new__(cls, 
            my_min = [0.0, 0.0, 0.0],
            my_max = [10.0, 10.0, 10.0],
            my_cutoff = 1.5,
            AA = False,
        ):

        xdim = int( np.ceil ( (my_max[0] - my_min[0])/my_cutoff ))
        ydim = int( np.ceil ( (my_max[1] - my_min[1])/my_cutoff ))
        zdim = int( np.ceil ( (my_max[2] - my_min[2])/my_cutoff ))

        if xdim == 0:
            xdim = 1
        if ydim == 0:
            ydim = 1
        if zdim == 0:
            zdim = 1
        shape = (xdim, ydim, zdim)

        obj = np.zeros(shape, dtype = object).view(cls)
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self._I = None

class Cell( np.ndarray ):

    def __new__(cls, 
            my_min = [0.0, 0.0, 0.0],
            my_max = [10.0, 10.0, 10.0],
            my_cutoff = 1.5,
        ):

        xdim = int( np.ceil ( (my_max[0] - my_min[0])/my_cutoff ))
        ydim = int( np.ceil ( (my_max[1] - my_min[1])/my_cutoff ))
        zdim = int( np.ceil ( (my_max[2] - my_min[2])/my_cutoff ))

        if xdim == 0:
            xdim = 1
        if ydim == 0:
            ydim = 1
        if zdim == 0:
            zdim = 1
        shape = (xdim, ydim, zdim)
        obj = np.zeros(shape, dtype = object )
        return obj

    def __init__(self, 
            my_min = [0.0, 0.0, 0.0],
            my_max = [10.0, 10.0, 10.0],
            my_cutoff = 1.5,
            AA = False):
        """docstring for __init__"""

        self.AA = AA

        self.my_xmin = my_min[0]
        self.my_ymin = my_min[1]
        self.my_zmin = my_min[2]
        
        self.my_xmax = my_max[0]
        self.my_ymax = my_max[1]
        self.my_zmax = my_max[2]

        self.my_cutoff = my_cutoff

        self.xdim = int( np.ceil ( (self.my_xmax - self.my_xmin)/my_cutoff ))
        self.ydim = int( np.ceil ( (self.my_ymax - self.my_ymin)/my_cutoff ))
        self.zdim = int( np.ceil ( (self.my_zmax - self.my_zmin)/my_cutoff ))

        self[:,:,:] = []

    def __array_finalize__(self, obj):
        if obj is None:
            return

    @staticmethod
    def from_xyz( fil, co = 2.0, in_AA = False, out_AA = False ):

        ats = []

        if not in_AA:
            co /= a0

        for f_ in open(fil ).readlines()[2:]:
            el, x, y, z = f_.split()
            x, y, z = map(float, [x,y,z] )
            if in_AA and not out_AA:
                x, y, z = map( lambda x: x/a0, [x,y,z] )
            ats.append( molecules.Atom( element = el, x=  x, y = y, z = z ))

        cell = Cell( my_min = [ min( ats, key = lambda x: x.x ).x,
                         min( ats, key = lambda x: x.y ).y,
                         min( ats, key = lambda x: x.z ).z],
              my_max = [ max( ats, key = lambda x: x.x ).x,
                         max( ats, key = lambda x: x.y ).y,
                         max( ats, key = lambda x: x.z ).z],
              my_cutoff = co,
              AA = out_AA,
              )
        for at in ats:
            cell.add(at)

        for at in cell:
            if len(at.Molecule) == 0:
                at.Molecule.append( at )
            cell.build_molecules( current = at, closeby = cell.get_closest(at) )
        return cell

    def build_molecules(self, current = None, visited = [],
            closeby = [], max_dist = 1.46, AA = False):
        visited.append( current )
        if not self.AA:
            max_dist /= a0

        if closeby == []:
            return

        for at in closeby:
            if at in current.Molecule:
                continue

            if max_dist < current.dist_to_atom( at ):
                continue

            current.Molecule.append( at )
            at.Molecule = current.Molecule
            close = [a for a in self.get_closest( at ) if a not in visited ]
            self.build_molecules( current = at, closeby = close, visited = visited )
 
    def add_atom( self, atom ):
        assert type( atom ) == molecules.Atom
        self.add( atom )
    def add_molecule( self, mol ):
        assert isinstance( mol, molecules.Molecule )
        for at in mol:
            self.add( at )

    def add_cluster( self, clus ):
        assert isinstance( clus, molecules.Cluster )
        for item in [at for mol in clus for at in mol ]:
            self.add( item )

    def add(self, item ):
        x_ind, y_ind, z_ind = self.get_index( item )

        if item not in self[ x_ind, y_ind, z_ind ]:
            self[ x_ind, y_ind, z_ind ].append( item )

    def __iter__(self):
        for i in range(len(self)):
            for j in range(len(self[i])):
                for k in range(len(self[i][j])):
                    for at in self[i][j][k]:
                        yield at

    def get_closest( self, item ):
        """
Return the closest items, 
to iterate not over whole cell box but closest

    >>> c = Cell( my_cutoff = 1.5 )
    >>> a1 = Atom( element = 'H', x = 1.4 ) #in the x index 0
    >>> a2 = Atom( element = 'O', x = 1.6 ) #in the x index 1
    >>> c.add( a1 )
    >>> c.add( a2 )
    >>> c.get_closest( a1 ) #Will return list where only the Oxygen exists
    [<molecules.Atom at 0x0xah5ah3h5] 
        """
        x_ind, y_ind, z_ind = self.get_index( item )
        tmp_list = []
        for i in range( x_ind -1, x_ind + 2 ):
            for j in range( y_ind -1, y_ind + 2 ):
                for k in range( z_ind -1, z_ind + 2 ):
                    try:
                        for at in self[i][j][k]:
                            if at == item or at in tmp_list:
                                continue
                            tmp_list.append( at )
                    except IndexError:
                        pass
        return tmp_list

    def update(self):

        tmp = []
        for x in range(len(self)):
            for y in range(len(self[x])):
                for z in range(len(self[x][y])):
                    for item in self[x][y][z]:
                        tmp.append( item )
                        self[x][y][z].remove( item )
        for item in tmp:
            self.add( item )


    def plot(self):
        """
Plot all Atoms in the cell

.. code:: python

    >>> cell = Cell()
    >>> cell.add( Atom(element = 'H' ))
    >>> cell.plot()
"""
#Plot water molecule in green and  nice xyz axis
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d' )
        ax.plot( [0, 1, 0, 0, 0, 0], [0,0 ,0,1,0,0], [0,0,0,0,0,1] )
        ax.text( 1.1, 0, 0, "X", color = 'red' )
        ax.text( 0, 1.1, 0, "Y", color = 'red' )
        ax.text( 0, 0, 1.1, "Z", color = 'red' )

        for at in self:
            ax.plot( [at.x], [at.y], [at.z], at.Molecule.style[at.element], linewidth= at.Molecule.linewidth[at.element] )
        ax.set_zlim3d( -5,5)
        plt.xlim(-5,5)
        plt.ylim(-5,5)
        plt.show()

    def get_index( self, item ):
        """
Return the x, y, and z index for cell for this item,

    >>> c = Cell( my_cutoff = 1.5 )
    >>> a1 = Atom( element = 'H', x = 1.4 ) #in the x index 0
    >>> print c.get_index( a1 )
    (0, 0, 0,)
"""
        assert self.my_xmin <= item[0] <= self.my_xmax
        assert self.my_ymin <= item[1] <= self.my_ymax
        assert self.my_zmin <= item[2] <= self.my_zmax

        tmp_xmin = item[0] - self.my_xmin
        tmp_ymin = item[1] - self.my_ymin
        tmp_zmin = item[2] - self.my_zmin

        x_ind = int( np.floor( tmp_xmin /  self.my_cutoff))
        y_ind = int( np.floor( tmp_ymin /  self.my_cutoff))
        z_ind = int( np.floor( tmp_zmin /  self.my_cutoff))
        
        return (x_ind, y_ind, z_ind)
