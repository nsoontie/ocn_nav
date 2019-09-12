# Script to 'unstagger' velocity fields and write to new netcdf files
# Also copies T/W/bio files to Agrid dir and links back to original fdir
# Will not operate on a file if it has already been unstaggered/linked
# Author: Nancy Soontiens Feb 2018


import netCDF4 as nc
import numpy as np
import os
import glob
import sys
import subprocess

MODEL_PATH = '/data/hdd/salishsea/model_orig'
SAVEDIR = '/data/hdd/salishsea/model_Agrid'

def main(ufiles, vfiles, tfiles, wfiles, biofiles, savedir=SAVEDIR):
    """ Unstagger and save U velocity in ufiles nd V velocity in vfiles"""
    for fname in ufiles:
        print(fname)
        build_unstaggered_U_file(fname, savedir=savedir)
    for fname in vfiles:
        print(fname)
        build_unstaggered_V_file(fname, savedir=savedir)
    move_other_files(tfiles,wfiles,biofiles)

def build_unstaggered_U_file(fname, savedir=SAVEDIR):
    # Copy file
    basename = os.path.basename(fname)
    dirname = os.path.dirname(fname)
    datedir = os.path.basename(os.path.dirname(fname))
    newdir = os.path.join(savedir, datedir)
    if not os.path.exists(newdir):
        os.makedirs(newdir)
    tfile = os.path.join(dirname, '{}grid_T.nc'.format(basename[:-9]))
    newname = os.path.join(newdir, '{}_Agrid.nc'.format(basename[:-3]))
    if not os.path.exists(newname):
        subprocess.call(['cp', fname, newname])
        # Open file and load U
        fU = nc.Dataset(newname, 'r+')
        uold = fU.variables['vozocrtx']
        # Unstagger
        unew = np.empty_like(uold[:])
        unew[...,1:] = np.add(uold[...,:-1], uold[...,1:]) / 2
        uold[:] = unew[:]
        uold.comment = 'unstaggered to A grid (T cells)'
        # Update other variables (nav_lon, nav_lat, depthu, area)
        fT = nc.Dataset(tfile)
        lonT = fT.variables['nav_lon']
        latT = fT.variables['nav_lat']
        areaT = fT.variables['area']
        depthT = fT.variables['deptht']
        bounds_lonT = fT.variables['bounds_lon']
        bounds_latT = fT.variables['bounds_lat']
        lon = fU.variables['nav_lon']
        lon[:] = lonT[:]
        lat = fU.variables['nav_lat']
        lat[:] = latT[:]
        area = fU.variables['area']
        area[:] = areaT[:]
        depth = fU.variables['depthu']
        depth[:] = depthT[:]
        depth.long_name = 'Vertical T levels'
        depth.bounds='deptht_bounds'
        bounds_lon = fU.variables['bounds_lon']
        bounds_lon[:] = bounds_lonT[:]
        bounds_lat = fU.variables['bounds_lat']
        bounds_lat[:] = bounds_latT[:]
        fT.close()
        # Update global attributes
        fU.comment1 = 'U unstaggered to A-grid'
        fU.comment2 = 'Created with script {}'.format(os.path.join(
                                                      os.getcwd(),sys.argv[0]))
        # Save (Close file)
        fU.close()
        # Use nco tools to rename
        subprocess.call(['ncrename', '-d', 'depthu,deptht',
                         '-v', 'depthu,depth',
                         '-v', 'depthu_bounds,deptht_bounds',
                          newname])
        subprocess.call(['ncatted', '-a', 'name,depth,d,,',
                         '-a',
                         'coordinates,vozocrtx,m,c,time_centered depth nav_lat nav_lon',
                         newname])
        

def build_unstaggered_V_file(fname, savedir=SAVEDIR):
    # Copy file
    basename = os.path.basename(fname)
    dirname = os.path.dirname(fname)
    datedir = os.path.basename(os.path.dirname(fname))
    newdir = os.path.join(savedir, datedir)
    if not os.path.exists(newdir):
        os.makedirs(newdir)            
    tfile = os.path.join(dirname, '{}grid_T.nc'.format(basename[:-9]))
    newname = os.path.join(newdir, '{}_Agrid.nc'.format(basename[:-3]))
    if not os.path.exists(newname):
        subprocess.call(['cp', fname, newname])
        # Open file and load V
        fV = nc.Dataset(newname, 'r+')
        vold = fV.variables['vomecrty']
        # Unstagger
        vnew = np.empty_like(vold[:])
        vnew[...,1:, :] = np.add(vold[...,:-1, :], vold[...,1:, :]) / 2
        vold[:] = vnew[:]
        vold.comment = 'unstaggered to A grid (T cells)'
        # Update other variables (nav_lon, nav_lat, depthv, area, ...)
        fT = nc.Dataset(tfile)
        lonT = fT.variables['nav_lon']
        latT = fT.variables['nav_lat']
        areaT = fT.variables['area']
        depthT = fT.variables['deptht']
        bounds_lonT = fT.variables['bounds_lon']
        bounds_latT = fT.variables['bounds_lat']
        lon = fV.variables['nav_lon']
        lon[:] = lonT[:]
        lat = fV.variables['nav_lat']
        lat[:] = latT[:]
        area = fV.variables['area']
        area[:] = areaT[:]
        depth = fV.variables['depthv']
        depth[:] = depthT[:]
        depth.long_name = 'Vertical T levels'
        depth.bounds='deptht_bounds'       
        bounds_lon = fV.variables['bounds_lon']
        bounds_lon[:] = bounds_lonT[:]
        bounds_lat = fV.variables['bounds_lat']
        bounds_lat[:] = bounds_latT[:]
        fT.close()
        # Update global attributes
        fV.comment1 = 'V unstaggered to A-grid'
        fV.comment2 = 'Created with script {}'.format(os.path.join(
                                                      os.getcwd(),sys.argv[0]))
        # Save (Close file)
        fV.close()
        # Use nco tools to rename
        subprocess.call(['ncrename', '-d', 'depthv,deptht',
                         '-v', 'depthv,depth',
                         '-v', 'depthv_bounds,deptht_bounds',
                          newname])
        subprocess.call(['ncatted', '-a', 'name,depth,d,,',
                         '-a',
                         'coordinates,vomecrty,m,c,time_centered depth nav_lat nav_lon',
                         newname])

        
def move_other_files(tfiles, wfiles, biofiles, savedir=SAVEDIR):
    """Move files to Agrid dir for Ocean Navigator but link back to orig
    Ocean Navigator thredds won't work with links"""
    for fnames in [tfiles, wfiles, biofiles]:
        for fname in fnames:
            if not os.path.islink(fname):
                basename = os.path.basename(fname)
                dirname = os.path.dirname(fname)
                datedir = os.path.basename(os.path.dirname(fname))
                newname = os.path.join(savedir,datedir,basename)
                subprocess.call(['mv', fname, newname])
                subprocess.call(['ln', '-s', newname, dirname])
                                         
    
        
if __name__ == '__main__':
    ufiles = glob.glob(os.path.join(MODEL_PATH, '*/*_grid_U.nc'))
    vfiles = glob.glob(os.path.join(MODEL_PATH, '*/*_grid_V.nc'))
    tfiles = glob.glob(os.path.join(MODEL_PATH, '*/*_grid_T.nc'))
    wfiles = glob.glob(os.path.join(MODEL_PATH, '*/*_grid_W.nc'))  
    biofiles = glob.glob(os.path.join(MODEL_PATH, '*/*_ptrc_T.nc'))
    main(ufiles, vfiles, tfiles, wfiles, biofiles)
