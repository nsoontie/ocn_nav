# Script to rotate the riops forecast files that did not work with nco tools.

import os
import subprocess

import netCDF4 as nc
import numpy as np

source='/data/hdd/riops/riopsf/'
grid='/data/hdd/grids/riops/grid_angle.nc'
savedir='/data/hdd/riops/riopsf/rotated/'
files = ['2016071712_006_3D_nps.nc',
         '2016071718_003_3D_nps.nc',
         '2016071818_006_3D_nps.nc',
         '2016072012_003_3D_nps.nc',
         #'2016072400_009_3D_nps.nc' # Skipping this file because it is fully corrupt and not used in the navigator.
]

problem_files = ['2016071818_006_3D_nps.nc',
                 '2016072012_003_3D_nps.nc',
                 '2016072400_009_3D_nps.nc']
sample_file='2016072400_006_3D_nps.nc'   

g = nc.Dataset(grid)
sin_alpha=g.variables['sin_alpha'][:]
cos_alpha=g.variables['cos_alpha'][:]

# Shape
f = nc.Dataset(os.path.join(source, sample_file))
usample=f.variables['vozocrtx'][:]
f.close()

for fname in files:
    print('Processing {}'.format(fname))
    new='{}_cardinal_velocity.nc'.format(fname[:-3])
    savefile=os.path.join(savedir,new)
    print(savefile)
    subprocess.call(['cp', os.path.join(source, fname), savefile])
    f = nc.Dataset(savefile,'r+')
    if fname in problem_files:
        ev = np.ma.masked_array(np.full_like(usample, 1.e20, dtype=np.float64), mask=np.ones_like(usample))
        nv = np.ma.masked_array(np.full_like(usample, 1.e20, dtype=np.float64), mask=np.ones_like(usample))
    else:
        u = f.variables['vozocrtx'][:]
        v = f.variables['vomecrty'][:]
        mask=u.mask
        ev = np.ma.masked_array(u*cos_alpha - v*sin_alpha, mask=mask)
        nv = np.ma.masked_array(u*sin_alpha + v*cos_alpha, mask=mask)    
    east_vel=f.createVariable('east_vel', np.float64,
                              ('time', 'depth', 'yc', 'xc'),
                               fill_value=1.e20)
    north_vel=f.createVariable('north_vel', np.float64,
                               ('time', 'depth', 'yc','xc'),
                               fill_value=1.e20)
    east_vel[:]=ev
    north_vel[:]=nv
    # Meta data
    east_vel.cell_methods="time: mean (interval: 3 hours)"
    east_vel.coordinates="longitude latitude time"
    east_vel.grid_mapping="polar_stereographic"
    east_vel.missing_value=1.e20
    east_vel.nomvar="UUW"
    east_vel.stdoffset=0.
    east_vel.stdscale=1.
    east_vel.units="m s-1"
    east_vel.valid_max=20.
    east_vel.valid_min=-20.
    east_vel.long_name="ocean current in eastward direction"
    east_vel.short_name="east_vel"
    east_vel.standard_name="ocean current eastward"

    north_vel.cell_methods="time: mean (interval: 3 hours)"
    north_vel.coordinates="longitude latitude time"
    north_vel.grid_mapping="polar_stereographic"
    north_vel.missing_value=1.e20
    north_vel.nomvar="VVW"
    north_vel.stdoffset=0.
    north_vel.stdscale=1.
    north_vel.units="m s-1"
    north_vel.valid_max=20.
    north_vel.valid_min=-20.
    north_vel.long_name="ocean current in northward direction"
    north_vel.short_name="north_vel"
    north_vel.standard_name="ocean current northward"
                                                        
    f.close()
    
    # Delete other variables and compress
    delvar='votemper,vosaline,vozocrtx,vomecrty,polar_stereographic'
    subprocess.call(['ncks', '-O', '-x', '-v', delvar, savefile, savefile])
    subprocess.call(['ncks', '-4', '-L4', '-O', savefile, savefile])
    delatts=['source','institution','contact','history']
    for att in delatts:
        subprocess.call(['ncatted', '-h', '-a', '{},global,d,,'.format(att),
                         savefile])
    comment='DERIVED PRODUCT - created with script riops_bad.py'
    subprocess.call(['ncatted', '-h', '-a',
                     'comment,global,o,c,{}'.format(comment),
                     savefile])
    
g.close()
