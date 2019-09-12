# Script to interpolate from Agrid to Cgrid
# Also links T/W/bio files to Cgrid_interp dir
# Usage: python make_Cgrid_interp.pu YYYYMMDD 
# Author: Nancy Soontiens May 2018


import netCDF4 as nc
import numpy as np
import os
import glob
import sys
import subprocess
import scipy.interpolate as interp

ORIG_PATH = '/data/hdd/salishsea/model_orig'
AGRID_PATH='/data/hdd/salishsea/model_Agrid/'
SAVEDIR = '/data/hdd/salishsea/model_Cgrid_interp'

MESH='/data/hdd/salishsea/grid/mesh_mask201702.nc'

def interp_to_Cgrid(date, lonU, latU, lonV, latV, lonT, latT,
                    umask, vmask):
    # Load
    AfU = nc.Dataset(glob.glob(os.path.join(AGRID_PATH, date,
                                            '*{}*U*.nc'.format(date)))[0])
    AfV = nc.Dataset(glob.glob(os.path.join(AGRID_PATH, date,
                                            '*{}*V*.nc'.format(date)))[0])
    Uagrid = AfU.variables['vozocrtx'][:]
    Vagrid = AfV.variables['vomecrty'][:]

    # Construct save files names
    dirsave = os.path.join(SAVEDIR,date)
    subprocess.call(['mkdir', '-p', dirsave])
    #U
    uorig = glob.glob(os.path.join(ORIG_PATH, date,
                                   '*{}*U*.nc'.format(date)))[0]
    print(uorig)
    basename = os.path.basename(uorig)
    saveu = os.path.join(SAVEDIR,date,basename)
    subprocess.call(['cp', uorig, saveu])
    # V
    vorig = glob.glob(os.path.join(ORIG_PATH, date,
                                   '*{}*V*.nc'.format(date)))[0]
    basename = os.path.basename(vorig)
    savev = os.path.join(SAVEDIR,date,basename)
    subprocess.call(['cp', vorig, savev])
    # Linking other files
    # T
    torig = glob.glob(os.path.join(ORIG_PATH, date,
                                   '*{}*grid_T*.nc'.format(date)))[0]
    basename = os.path.basename(torig)
    savet = os.path.join(SAVEDIR,date,basename)
    subprocess.call(['ln', '-fs', torig, savet])
    # W
    worig = glob.glob(os.path.join(ORIG_PATH, date,
                                   '*{}*grid_W*.nc'.format(date)))[0]
    basename = os.path.basename(worig)
    savew = os.path.join(SAVEDIR,date,basename)
    subprocess.call(['ln', '-fs', worig, savew])
    # bio
    borig = glob.glob(os.path.join(ORIG_PATH, date,
                                   '*{}*ptrc_T*.nc'.format(date)))[0]
    basename = os.path.basename(borig)
    saveb = os.path.join(SAVEDIR,date,basename)
    subprocess.call(['ln', '-fs', borig, saveb])              

    print('Beginning interpolation {}'.format(date))
    # interpolate
    U_cgrid = np.zeros_like(Uagrid)
    V_cgrid = np.zeros_like(Vagrid)
    points = (lonT.flatten(), latT.flatten())
    for t in range(U_cgrid.shape[0]):
        print('{}/{}'.format(t, U_cgrid.shape[0]))
        for k in range(U_cgrid.shape[1]):
            U_cgrid[t,k,:,:] = interp.griddata(points,
                                               Uagrid[t,k,:,:].flatten(),
                                               (lonU, latU))
            V_cgrid[t,k,:,:] = interp.griddata(points,
                                               Vagrid[t,k,:,:].flatten(),
                                               (lonV, latV))
    print('Finished interpolation {}'.format(date))
    AfU.close()
    AfV.close()
    # Saving
    CfU = nc.Dataset(saveu, 'r+')
    U = CfU.variables['vozocrtx']
    U[:] = np.ma.masked_array(U_cgrid[:],
                              mask=np.logical_not(umask + np.zeros_like(U_cgrid)))
    U.comment = 'Interpolated from A grid to C grid'
    CfU.comment = 'Created with script {}'.format(os.path.join(os.getcwd(),
                                                               sys.argv[0]))
    CfU.close()

    CfV = nc.Dataset(savev,'r+')
    V = CfV.variables['vomecrty']
    V[:] = np.ma.masked_array(V_cgrid[:],
                              mask=np.logical_not(vmask + np.zeros_like(V_cgrid)))
    V.comment = 'Interpolated from A grid to C grid'
    CfV.comment = 'Created with script {}'.format(os.path.join(os.getcwd(),
                                                               sys.argv[0]))
    CfV.close()
        
if __name__ == '__main__':
    date = sys.argv[1]
    mesh = nc.Dataset(MESH)
    lonu = mesh.variables['glamu'][0,:]
    lonv = mesh.variables['glamv'][0,:]
    lont = mesh.variables['glamt'][0,:]
    latu = mesh.variables['gphiu'][0,:]
    latv = mesh.variables['gphiv'][0,:]
    latt = mesh.variables['gphit'][0,:]
    umask = mesh.variables['umask'][:]
    vmask = mesh.variables['vmask'][:]
    interp_to_Cgrid(date, lonu, latu, lonv, latv, lont, latt, umask, vmask)
