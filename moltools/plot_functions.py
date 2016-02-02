#!/usr/bin/env python
# -*- coding: utf-8 -*-

__all__ = []

from molecules import Molecule

def plot_maya( mol, key = lambda x: (x[0].r, x[1].r, x[2].r),
        copy = True ):
    """Wrapper function to plot quivers of beta 

    For Residue
    
    at1 is the first atom to center around
    at2 is the atom which will lie in the z axis together with at1
    at3 will be projected to lie in the zx-axis
    """

    assert isinstance( mol, Molecule )
    from mayavi import mlab
    
    if copy:
        copy = mol.copy()
    else:
        copy = mol
    copy.populate_bonds()
    p1, p2, p3 = key( copy )
    v, t1, t2, t3 = utilz.center_and_xz( p1, p2, p3 )

# plotting section
    copy.translate_by_r( v )
    copy.inv_rotate( t1, t2, t3 )

    x, y, z = utilz.E_at_sphere(r_points=5)
    bv = utilz.b_at_sphere( utilz.ut2s(copy.p.b) , x, y, z )

    mlab.figure(figure=None, bgcolor=(1,1,1), fgcolor=None, engine=None, size=(400, 350))
    mlab.quiver3d( x, y, z, bv[...,0], bv[...,1], bv[...,2],
            colormap = 'BrBG' )
#Plot bonds
    for each in copy.bond_dict:
        for key in copy.bond_dict[ each ]:
            mlab.plot3d( [key.x, each.x],
                     [key.y, each.y],
                     [key.z, each.z], color = (0,0,0,) )
    for i in copy:
        mlab.points3d([i.x], [i.y], [i.z], 
                color = color_dict[ i.element ],
                resolution = 50,
                scale_factor = scale_factor_dict[ i.element] )
