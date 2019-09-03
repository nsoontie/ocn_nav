# Module to reproduce angle and rot_rep code in NEMO
# http://forge.ipsl.jussieu.fr/little_nemo/browser/trunk/NEMOGCM/NEMO/OPA_SRC/SBC/geo2ocean.F90

import numpy as np
import netCDF4 as nc
import math

def angles(mesh_file, orca, fill_value):
    """
    !!----------------------------------------------------------------------
          !!                  ***  ROUTINE angle  ***
          !!
          !! ** Purpose :   Compute angles between model grid lines and the North direction
          !!
          !! ** Method  :
          !!
          !! ** Action  :   Compute (gsint, gcost, gsinu, gcosu, gsinv, gcosv, gsinf, gcosf) arrays:
          !!      sinus and cosinus of the angle between the north-south axe and the
          !!      j-direction at t, u, v and f-points
          !!
          !! History :
          !!   7.0  !  96-07  (O. Marti )  Original code
          !!   8.0  !  98-06  (G. Madec )
          !!   8.5  !  98-06  (G. Madec )  Free form, F90 + opt.
          !!   9.2  !  07-04  (S. Masson)  Add T, F points and bugfix in cos lateral boundary
          !!----------------------------------------------------------------------
    :arg str mesh_file: path to netcdf mesh file
    
    :returns:   
    """    
    # Load mesh
    f = nc.Dataset(mesh_file)
    glamt = f.variables['glamt'][:]
    if glamt.ndim > 2:
        glamt = f.variables['glamt'][0,...]
        glamu = f.variables['glamu'][0,...]
        glamv = f.variables['glamv'][0,...]
        glamf = f.variables['glamf'][0,...]
        gphit = f.variables['gphit'][0,...]
        gphiu = f.variables['gphiu'][0,...]
        gphiv = f.variables['gphiv'][0,...]
        gphif = f.variables['gphif'][0,...]
    else:
        glamt = f.variables['glamt'][:]
        glamu = f.variables['glamu'][:]
        glamv = f.variables['glamv'][:]
        glamf = f.variables['glamf'][:]
        gphit = f.variables['gphit'][:]
        gphiu = f.variables['gphiu'][:]
        gphiv = f.variables['gphiv'][:]
        gphif = f.variables['gphif'][:]
                                                                        
    f.close()
    Ny,Nx = glamt.shape

    # Allocate empty arrays
    gsint = np.empty((Ny,Nx))
    gsinu = np.empty((Ny,Nx))
    gsinv = np.empty((Ny,Nx))
    gsinf = np.empty((Ny,Nx))
    
    gcost = np.empty((Ny,Nx))
    gcosu = np.empty((Ny,Nx))
    gcosv = np.empty((Ny,Nx))
    gcosf = np.empty((Ny,Nx))
            
    # Useful constants
    rad = np.pi/180.
    rpi = np.pi

    #      ! ============================= !
    #      ! Compute the cosinus and sinus !
    #      ! ============================= !
    for jj in range(1, Ny-1):
        for ji in range(1, Nx): 
            # north pole direction & modulous (at t-point)
            zlam = glamt[jj,ji]
            zphi = gphit[jj,ji]
            zxnpt = 0. - 2. * np.cos( rad*zlam ) * np.tan( rpi/4. - rad*zphi/2. )
            zynpt = 0. - 2. * np.sin( rad*zlam ) * np.tan( rpi/4. - rad*zphi/2. )
            znnpt = zxnpt*zxnpt + zynpt*zynpt
            # north pole direction & modulous (at u-point)
            zlam = glamu[jj,ji]
            zphi = gphiu[jj,ji]
            zxnpu = 0. - 2. * np.cos( rad*zlam ) * np.tan( rpi/4. - rad*zphi/2. )
            zynpu = 0. - 2. * np.sin( rad*zlam ) * np.tan( rpi/4. - rad*zphi/2. )
            znnpu = zxnpu*zxnpu + zynpu*zynpu
   
            # north pole direction & modulous (at v-point)
            zlam = glamv[jj,ji]
            zphi = gphiv[jj,ji]
            zxnpv = 0. - 2. * np.cos( rad*zlam ) * np.tan( rpi/4. - rad*zphi/2. )
            zynpv = 0. - 2. * np.sin( rad*zlam ) * np.tan( rpi/4. - rad*zphi/2. )
            znnpv = zxnpv*zxnpv + zynpv*zynpv

            # north pole direction & modulous (at f-point)
            zlam = glamf[jj,ji]
            zphi = gphif[jj,ji]
            zxnpf = 0. - 2. * np.cos( rad*zlam ) * np.tan( rpi/4. - rad*zphi/2. )
            zynpf = 0. - 2. * np.sin( rad*zlam ) * np.tan( rpi/4. - rad*zphi/2. )
            znnpf = zxnpf*zxnpf + zynpf*zynpf
   
            # j-direction: v-point segment direction (around t-point)
            zlam = glamv[jj,ji  ]
            zphi = gphiv[jj,ji  ]
            zlan = glamv[jj-1,ji]
            zphh = gphiv[jj-1,ji]
            zxvvt =  2. * np.cos( rad*zlam ) * np.tan( rpi/4. - rad*zphi/2. ) - \
                     2. * np.cos( rad*zlan ) * np.tan( rpi/4. - rad*zphh/2. )
            zyvvt =  2. * np.sin( rad*zlam ) * np.tan( rpi/4. - rad*zphi/2. ) - \
                     2. * np.sin( rad*zlan ) * np.tan( rpi/4. - rad*zphh/2. )
            znvvt = np.sqrt( znnpt * ( zxvvt*zxvvt + zyvvt*zyvvt )  )
            znvvt = max( znvvt, 1.e-14)
            
            # j-direction: f-point segment direction (around u-point)
            zlam = glamf[jj,ji  ]
            zphi = gphif[jj,ji  ]
            zlan = glamf[jj-1,ji]
            zphh = gphif[jj-1,ji]
            zxffu =  2. * np.cos( rad*zlam ) * np.tan( rpi/4. - rad*zphi/2. ) - \
                     2. * np.cos( rad*zlan ) * np.tan( rpi/4. - rad*zphh/2. )
            zyffu =  2. * np.sin( rad*zlam ) * np.tan( rpi/4. - rad*zphi/2. ) - \
                     2. * np.sin( rad*zlan ) * np.tan( rpi/4. - rad*zphh/2. )
            znffu = np.sqrt( znnpu * ( zxffu*zxffu + zyffu*zyffu )  )
            znffu = max( znffu, 1.e-14 )
    
            # i-direction: f-point segment direction (around v-point)
            zlam = glamf[jj  ,ji]
            zphi = gphif[jj  ,ji]
            zlan = glamf[jj,ji-1]
            zphh = gphif[jj,ji-1]
            zxffv =  2. * np.cos( rad*zlam ) * np.tan( rpi/4. - rad*zphi/2. ) - \
                     2. * np.cos( rad*zlan ) * np.tan( rpi/4. - rad*zphh/2. )
            zyffv =  2. * np.sin( rad*zlam ) * np.tan( rpi/4. - rad*zphi/2. ) - \
                     2. * np.sin( rad*zlan ) * np.tan( rpi/4. - rad*zphh/2. )
            znffv = np.sqrt( znnpv * ( zxffv*zxffv + zyffv*zyffv )  )
            znffv = max( znffv, 1.e-14 )

            # j-direction: u-point segment direction (around f-point)
            zlam = glamu[jj+1,ji]
            zphi = gphiu[jj+1,ji]
            zlan = glamu[jj,ji  ]
            zphh = gphiu[jj,ji  ]
            zxuuf =  2. * np.cos( rad*zlam ) * np.tan( rpi/4. - rad*zphi/2. ) - \
                     2. * np.cos( rad*zlan ) * np.tan( rpi/4. - rad*zphh/2. )
            zyuuf =  2. * np.sin( rad*zlam ) * np.tan( rpi/4. - rad*zphi/2. ) - \
                     2. * np.sin( rad*zlan ) * np.tan( rpi/4. - rad*zphh/2. )
            znuuf = np.sqrt( znnpf * ( zxuuf*zxuuf + zyuuf*zyuuf )  )
            znuuf = max( znuuf, 1.e-14 )
   
            # cosinus and sinus using scalar and vectorial products
            gsint[jj,ji] = ( zxnpt*zyvvt - zynpt*zxvvt ) / znvvt
            gcost[jj,ji] = ( zxnpt*zxvvt + zynpt*zyvvt ) / znvvt

            gsinu[jj,ji] = ( zxnpu*zyffu - zynpu*zxffu ) / znffu
            gcosu[jj,ji] = ( zxnpu*zxffu + zynpu*zyffu ) / znffu

            gsinf[jj,ji] = ( zxnpf*zyuuf - zynpf*zxuuf ) / znuuf
            gcosf[jj,ji] = ( zxnpf*zxuuf + zynpf*zyuuf ) / znuuf
   
            # (caution, rotation of 90 degres)
            gsinv[jj,ji] = ( zxnpv*zxffv + zynpv*zyffv ) / znffv
            gcosv[jj,ji] =-( zxnpv*zyffv - zynpv*zxffv ) / znffv
                                    


    #      ! =============== !
    #      ! Geographic mesh !
    #      ! =============== !
    for jj in range(1, Ny-1):
        for ji in range(1, Nx): 
            if np.mod( np.abs( glamv[jj,ji] - glamv[jj-1,ji] ), 360. ) < 1.e-8 :
                 gsint[jj,ji] = 0.
                 gcost[jj,ji] = 1.
            if np.mod( np.abs( glamf[jj,ji] - glamf[jj-1,ji] ), 360. ) < 1.e-8 :
                gsinu[jj,ji] = 0.
                gcosu[jj,ji] = 1.
            if        np.abs( gphif[jj,ji] - gphif[jj,ji-1] )          < 1.e-8 :
                gsinv[jj,ji] = 0.
                gcosv[jj,ji] = 1.
            if np.mod( np.abs( glamu[jj,ji] - glamu[jj+1,ji] ), 360. ) < 1.e-8 :
                gsinf[jj,ji] = 0.
                gcosf[jj,ji] = 1.
    
    #      ! =========================== !
    #      ! Lateral boundary conditions !
    #      ! =========================== !
   
    #      ! lateral boundary cond.: T-, U-, V-, F-pts, sgn
    #      CALL lbc_lnk( gcost, 'T', -1. )   ;   CALL lbc_lnk( gsint, 'T', -1. )
    #      CALL lbc_lnk( gcosu, 'U', -1. )   ;   CALL lbc_lnk( gsinu, 'U', -1. )
    #      CALL lbc_lnk( gcosv, 'V', -1. )   ;   CALL lbc_lnk( gsinv, 'V', -1. )
    #      CALL lbc_lnk( gcosf, 'F', -1. )   ;   CALL lbc_lnk( gsinf, 'F', -1. )
    if orca == 'T':     
        # Treatment of north fold
        print("Treating north fold")
        sign=-1
        gsint = orca_bc(gsint, sign, 'T', fill_value)
        gcost = orca_bc(gcost, sign, 'T', fill_value)
        gsinu = orca_bc(gsinu, sign, 'U', fill_value)
        gcosu = orca_bc(gcosu, sign, 'U', fill_value)
        gsinv = orca_bc(gsinv, sign, 'V', fill_value)
        gcosv = orca_bc(gcosv, sign, 'V', fill_value)
        gsinf = orca_bc(gsinf, sign, 'F', fill_value)
        gcosf = orca_bc(gcosf, sign, 'F', fill_value)
    else:
        # Set boundaries to closest inner
        # First row copied from second row
        for var in [gsint, gsinu, gsinv, gsinf, \
                    gcost, gcosu, gcosv, gcosf ]:
            var[:,0] = var[:,1]
        # Last row copied from second last row
        for var in [gsint, gsinu, gsinv, gsinf, \
                    gcost, gcosu, gcosv, gcosf ]:
            var[:,-1] = var[:,-2]                    
        # First column copied from second column
        for var in [gsint, gsinu, gsinv, gsinf, \
                    gcost, gcosu, gcosv, gcosf ]:
            var[0,:] = var[1,:]                
        # Last column copied from second last column
        for var in [gsint, gsinu, gsinv, gsinf, \
                    gcost, gcosu, gcosv, gcosf ]:
            var[-1,:] = var[-2,:]                

    return gsint, gsinu, gsinv, gsinf, \
           gcost, gcosu, gcosv, gcosf


def orca_bc(var, sign, case,fill_value):
    """Replicate orca north-fold bc for T-pivot
       Also, east/west cyclic"""
    varT = np.copy(var.T)
    Ny, Nx = var.shape
    # East/West cyclic first
    varT[0, :] = varT[Nx-2, :]
    varT[Nx-1, :] = varT[1, :]
    # South closed
    varT[:, 0] = fill_value
    # North fold
    jpiglo, jpjglo = Nx, Ny
    ipr2dj = 0
    ijpj = Ny
    ijpjm1 = Ny-1
    jpi = Nx
    jpj = Ny
    if case == 'T':
        for jl in  range(ipr2dj+1):
           for ji in range(2, jpiglo+1):
               ijt = jpiglo-ji+2
               varT[ji-1, ijpj+jl-1] = sign * varT[ijt-1, ijpj-2-jl-1]
        varT[0, ijpj-1] = sign * varT[2, ijpj-3]
        for ji in range(int(jpiglo/2)+1, jpiglo+1):
            ijt = jpiglo-ji+2
            varT[ji-1, ijpj-2] = sign * varT[ijt-1, ijpj-2]
    elif case =='U':
        for jl in range(ipr2dj+1):
            for ji in range (1, jpiglo):
                iju = jpiglo-ji+1
                varT[ji-1, ijpj+jl-1] = sign * varT[iju-1, ijpj-2-jl -1]
        varT[0       , ijpj-1] = sign * varT[    1   , ijpj-3]
        varT[jpiglo-1, ijpj-1] = sign * varT[jpiglo-2, ijpj-3]
        varT[0       , ijpj-2] = sign * varT[jpiglo-1, ijpj-2]
        for ji in range( int(jpiglo/2), jpiglo):
            iju = jpiglo-ji+1
            varT[ji-1, ijpjm1-1] = sign * varT[iju-1, ijpjm1-1]
    elif case == 'V':
        for jl in range( -1, ipr2dj +1):
            for ji in range( 2, jpiglo +1):
                ijt = jpiglo-ji+2
                varT[ji-1, ijpj+jl-1] = sign* varT[ijt -1, ijpj-3-jl-1]
        varT[0, ijpj-1]   = sign * varT[ 2 ,ijpj-3]
    elif case == 'F':
        for jl in range( -1, ipr2dj+1):
            for ji in range(1, jpiglo):
                iju = jpiglo-ji+1
                varT[ji-1, ijpj+jl-1] = sign * varT[iju-1, ijpj-3-jl-1]
        varT[   0    , ijpj-1] = sign * varT[    1   , ijpj-4]
        varT[jpiglo-1, ijpj-1] = sign * varT[jpiglo-2, ijpj-4]
        varT[jpiglo-1, ijpj-2] = sign * varT[jpiglo-2, ijpj-3]
        varT[   0    , ijpj-2] = sign * varT[    1   , ijpj-2]
    else:
        print("Case {} is invalid".format(case))
        exit
    return varT.T

