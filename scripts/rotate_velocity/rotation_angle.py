# Script to write angle of rotation between model x-grid lines and East-West.
# Usage: python rotation_angle.py gridfile varlon varlat xdimname ydimname \
#                                 method orca_grid outfile
# gridfile - name of netcdf file that contains latitude and longitude
# varlon - name of longitude variable in gridfile
# varlat - name of latitude variable in gridfile
# xdimname - name of x dimension in gridfile
# ydimname - name of y dimension in gridfile
# method - method for angle calculation - opa or ps
#   opa - based on NEMO code for angle calculation (geo2ocean.F90). Can only be
#   used when gridfile is a mesh file with lat/lon on u,v,t,f grids
#   ps - transform coordinates to polar stereographic.
#   Calculate angle between vector connecting grid corner and north and vector
#   connecting previous x point with grid corner. This method is not very
#   accurate and should only be used in absence of mesh files
#   opa is based on NEMO code for angle calculation (geo2ocean.F90). Can only
#   be used when gridfile is a mesh file with lat/lon on u,v,t,f grids
# orca_grid - indicates whether or not grid is an ORCA grid for
#  for special treatment of north fold. T - orca grid, F not orca grid
# Nancy Soontiens March 2018

import os
import netCDF4 as nc
import numpy as np
import sys
import subprocess
import opa_angle


def main(gridfile, varlon, varlat, xdimname, ydimname, method, orca, outfile,
         fill_value=-99.):
    """
    :arg str gridfile: name of netcdf file for reading grid information
    :arg str varlon: name of longitude variable in gridfile
    :arg str varlat: name of latitude variable in gridfile
    :arg str xdimname: name of x dimnesion in gridfile
    :arg str ydimname: name of y dimension in gridfile   
    :arg str method: opa or ps
    :arg str orca: T or F
    :arg str outfile: name netcdf file for saving results

    """
    f = nc.Dataset(gridfile)
    lons = f.variables[varlon]
    lats = f.variables[varlat]
    londata = lons[:].data if np.ma.is_masked(lons[:]) else lons[:]
    latdata = lata[:].data if np.ma.is_masked(lats[:]) else lats[:]
    if method=='opa':
        sintheta, _, _, _,\
        costheta, _, _, _ = opa_angle.angles(gridfile, orca, fill_value)
        theta = np.rad2deg(np.arctan2(sintheta, costheta))
        inds = np.where(np.logical_and(sintheta==fill_value, costheta==fill_value))
        theta[inds]=fill_value
    elif method=='ps':
        sintheta, costheta = grid_angle_ps(np.squeeze(londata),
                                           np.squeeze(latdata))
        theta = np.rad2deg(np.arctan2(sintheta, costheta))
    else:
        print('Invalid option: method={}'.format(method))
        exit
    save_netcdf(gridfile, varlon, varlat, xdimname, ydimname, outfile,
                lons, lats, theta, costheta, sintheta, method, orca,
                fill_value=fill_value)
    f.close()
    # Netcdf deflaction
    subprocess.call(['ncks', '-4', '-L4', '-O', outfile, outfile])


def save_netcdf(gridfile, varlon, varlat, xdimname, ydimname, outfile,
                lons, lats, theta, costheta, sintheta, method, orca,
                fill_value=-999.):
    """Save rotation angle and other metadate into netcdf file.
    :arg str gridfile: name of netcdf file for reading grid information
    :arg str varlon: name of longitude variable in gridfile
    :arg str varlat: name of latitude variable in gridfile
    :arg str xdimname: name of x dimnesion in gridfile
    :arg str ydimname: name of y dimension in gridfile
    :arg str outfile: name netcdf file for saving results 
   
    :arg lons: array of longitudes
    :type lons: numpy array (2D)

    :arg lats: array of latitudes
    :type lats: numpy array (2D)

    :arg theta: array of grid angle (degrees)
    :type theta: numpy array (2D)

    :arg costheta: array of cosine of grid angle
    :type costheta: numpy array (2D)

    :arg sintheta: array of sine of grid angle
    :type sintheta: nump array (2D)

    :arg str method: opa or ps
    :arg str orca: T or F  
    """
    f = nc.Dataset(outfile, 'w')
    # Create dimensions
    xdim = f.createDimension(xdimname, lons[:].shape[-1])
    ydim = f.createDimension(ydimname, lons[:].shape[-2])
    # Create variables
    latnc = f.createVariable(varlat, np.float32, (ydimname,xdimname,))
    lonnc = f.createVariable(varlon, np.float32, (ydimname,xdimname,))
    anglenc = f.createVariable('alpha', np.float32, (ydimname,xdimname,),fill_value=fill_value)
    cosalpha = f.createVariable('cos_alpha', np.float32, (ydimname,xdimname,),fill_value=fill_value)
    sinalpha = f.createVariable('sin_alpha', np.float32, (ydimname,xdimname,),fill_value=fill_value)
    # Write variables and attributes
    latdata = lats[:].data if np.ma.is_masked(lats[:]) else lats[:]
    latnc[:]= np.squeeze(latdata)
    if lats.ncattrs():
        latnc.setncatts({k: lats.getncattr(k) for k in lats.ncattrs()
                         if k != '_FillValue' and k != 'valid_min'
                         and k != 'valid_max'})
    else:
        latnc.standard_name = varlat
        latnc.long_name = 'latitude'
        latnc.units = 'degrees_north'
    londata = lons[:].data if np.ma.is_masked(lons[:]) else lons[:]
    lonnc[:]= np.squeeze(londata)
    if lons.ncattrs():
        lonnc.setncatts({k: lons.getncattr(k) for k in lons.ncattrs()
                         if k != '_FillValue' and k != 'valid_min'
                         and k != 'valid_max'})
    else:
        lonnc.standard_name = varlon
        lonnc.long_name = 'longitude'
        lonnc.units = 'degrees_east'
                                        
    anglenc[:] = theta[:]
    anglenc.standard_name = 'alpha'
    anglenc.long_name = 'Angle from East to model x-gridlines'
    anglenc.missing_value=fill_value
    if method=='opa':
        comment = 'Angle derived from opa angle using coordinates in {}'.format(gridfile)
    else:
        comment =('Angle derived from polar stereographic transformation using' 
                  'coordinates '
                  '{} {} in {}'.format(varlon, varlat, gridfile))
    anglenc.comment = comment
    anglenc.units = 'degrees CCW from East'

    cosalpha[:] = costheta[:]
    cosalpha.standard_name = 'cos_alpha'
    cosalpha.long_name = 'Cosine of grid angle alpha'
    cosalpha.missing_value=fill_value
    
    sinalpha[:] = sintheta[:]
    sinalpha.standard_name = 'sin_alpha'
    sinalpha.long_name = 'Sine of grid angle alpha'
    sinalpha.missing_value=fill_value
                
    # Global attributes
    f.comment = ('Created with call:' \
                 'python {} {} {} {} {} {} {} {} {}'.format(sys.argv[0],
                                                         gridfile,
                                                         varlon,
                                                         varlat,
                                                         xdimname,
                                                         ydimname,
                                                         method, orca,
                                                         outfile))
    f.references = os.path.realpath(__file__)
    f.close()

    
def grid_angle(lons, lats):
    """
    Find the angle between longitude grid lines and East.
    theta = pi/2 - (atan2(sin(dlong).cos(lat2),
           cos(lat1).sin(lat2) − sin(lat1).cos(lat2).cos(dlong)))
    theta can be used to rotate velocity to cardinal directions
    Reference: https://www.movable-type.co.uk/scripts/latlong.html
    ueast = u_x*cos(theta) - v_y*sin(theta)
    vnorth = u_x*sin(theta) + v_y*cos(theta)

    :arg lons: longitude in degrees shape (y, x)
    :type lons: 2D numpy array

    :arg lats: latitude in degrees shape (y, x)
    :type lats: 2D numpy array

    :returns: theta, angle measured from East to X gridlines
    """
    dlons = np.radians(lons[:,1:] - lons[:, 0:-1])
    lat1 = np.radians(lats[:, 0:-1])
    lat2 = np.radians(lats[:, 1:])
    theta=np.empty_like(lons)
    x = np.sin(dlons)*np.cos(lat2)
    y = np.cos(lat1) * np.sin(lat2) -\
        (np.sin(lat1)* np.cos(lat2) * np.cos(dlons))
    # Extend theta by copying first column
    theta[:, 1:] = np.arctan2(x,y)
    theta[:, 0] = theta[:,1]
    # Theta is the angle with North so subtract from pi/2 for angle with East
    return np.rad2deg(np.pi/2 -theta)


def grid_angle_jlines(lons, lats):
    """                                                                       
    Find the angle between latitude grid lines and North.                     
    theta = pi/2 - (atan2(sin(dlong).cos(lat2),                                
           cos(lat1).sin(lat2) − sin(lat1).cos(lat2).cos(dlong)))              
    theta can be used to rotate velocity to cardinal directions                
    Reference: https://www.movable-type.co.uk/scripts/latlong.html             
    ueast = u_x*cos(theta) - v_y*sin(theta)                                    
    vnorth = u_x*sin(theta) + v_y*cos(theta)                                   
                                                                               
    :arg lons: longitude in degrees shape (y, x)                          
    :type lons: 2D numpy array                                                 
                                                                               
    :arg lats: latitude in degrees shape (y, x)                                
    :type lats: 2D numpy array                                                 
                                                                               
    :returns: theta, angle measured from East to X gridline                    
    """
    dlons = np.radians(lons[1:,:] - lons[0:-1, :])
    lat1 = np.radians(lats[0:-1, :])
    lat2 = np.radians(lats[1:, :])
    theta=np.empty_like(lons)
    x = np.sin(dlons)*np.cos(lat2)
    y = np.cos(lat1) * np.sin(lat2) -\
        (np.sin(lat1)* np.cos(lat2) * np.cos(dlons))
    # Extend theta by copying first column
    theta[1:, :] = np.arctan2(x,y)
    theta[0, :] = theta[1, :]
    # Theta is a bearing (CW), but want CCW from North
    return np.rad2deg(-theta)


def grid_angle_ps(lons, lats):
    """Caclulate grid angle via polar stereographic projection.
       Angle between vector from grid corner to north and 
       vector from previous x point to grid corner"""
    Ny, Nx = lons.shape
    gsint = np.empty((Ny, Nx))
    gcost = np.empty((Ny, Nx))
    # Useful constants
    rad = np.pi/180.
    rpi = np.pi
                
    for jj in range(1, Ny-1):
        for ji in range(1, Nx):
            # Angle with north pole
            zlam = lons[jj,ji]
            zphi = lats[jj,ji]
            zxnpt = 0. - 2. * np.cos( rad*zlam ) * np.tan( rpi/4. - rad*zphi/2. )
            zynpt = 0. - 2. * np.sin( rad*zlam ) * np.tan( rpi/4. - rad*zphi/2. )
            znnpt = zxnpt*zxnpt + zynpt*zynpt
            # i direction segment
            zlam = lons[jj,ji  ]
            zphi = lats[jj,ji  ]
            zlan = lons[jj,ji-1]
            zphh = lats[jj,ji-1]
            zxuut =  2. * np.cos( rad*zlam ) * np.tan( rpi/4. - rad*zphi/2. ) - \
                     2. * np.cos( rad*zlan ) * np.tan( rpi/4. - rad*zphh/2. )
            zyuut =  2. * np.sin( rad*zlam ) * np.tan( rpi/4. - rad*zphi/2. ) - \
                     2. * np.sin( rad*zlan ) * np.tan( rpi/4. - rad*zphh/2. )
            znuut = np.sqrt( znnpt * ( zxuut*zxuut + zyuut*zyuut )  )
            znuut = max( znuut, 1.e-14 )

            gsint[jj,ji] = ( zxnpt*zxuut + zynpt*zyuut ) / znuut
            gcost[jj,ji] =-( zxnpt*zyuut - zynpt*zxuut ) / znuut
    # Boundaryies
    for var in [gsint, gcost]:
        var[:,0] = var[:,1]
        # Last row copied from second last row
        for var in [ gsint, gcost]:
            var[:,-1] = var[:,-2]
        # First column copied from second column
        for var in [gsint, gcost]:
            var[0,:] = var[1,:]
        # Last column copied from second last column
        for var in [gsint, gcost]:
            var[-1,:] = var[-2,:]
    return gsint, gcost
                                                                               
if __name__ == '__main__':
    gridfile, varlon, varlat, xdimname, ydimname, method, orca, outfile = \
        sys.argv[1:]
    main(gridfile, varlon, varlat, xdimname, ydimname, method, orca, outfile)
