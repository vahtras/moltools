import numpy as np

class Cell( list ):

    def __init__(self, 
            my_min = [0.0, 0.0, 0.0],
            my_max = [10.0, 10.0, 10.0],
            my_cutoff = 1.5):
        """docstring for __init__"""
        self.my_xmin = my_min[0]
        self.my_ymin = my_min[1]
        self.my_zmin = my_min[2]
        
        self.my_xmax = my_max[0]
        self.my_ymax = my_max[1]
        self.my_zmax = my_max[2]

        self.my_cutoff = my_cutoff

        xdim = int( np.ceil ( (self.my_xmax - self.my_xmin)/my_cutoff ))
        ydim = int( np.ceil ( (self.my_ymax - self.my_ymin)/my_cutoff ))
        zdim = int( np.ceil ( (self.my_zmax - self.my_zmin)/my_cutoff ))


        for i in range( xdim ):
            self.append( [] )
            for j in range( ydim ):
                self[i].append( [] )
                for k in range( zdim ):
                    self[i][j].append( [] )

    def add(self, item ):
        assert self.my_xmin <= item[0] <= self.my_xmax
        assert self.my_ymin <= item[1] <= self.my_ymax
        assert self.my_zmin <= item[2] <= self.my_zmax

        tmp_xmin = item[0] - self.my_xmin
        tmp_ymin = item[1] - self.my_ymin
        tmp_zmin = item[2] - self.my_zmin

        x_ind = int( np.floor( tmp_xmin /  self.my_cutoff))
        y_ind = int( np.floor( tmp_ymin /  self.my_cutoff))
        z_ind = int( np.floor( tmp_zmin /  self.my_cutoff))

        self[ x_ind ][ y_ind ][ z_ind ].append( item )
