"""
This file generates a new ROMS output file from scratch.
It is netCDF4 CF-compliant.

TODO: Find a way to loop over a list of strings to delete 
             all this crappy infinite command lines. >(
"""

from   netCDF4 import Dataset
from   matplotlib import path 
from   progress.bar import IncrementalBar
import numpy as np
import time
seaseaiceFillVal = 1.e+37

class seaice(object):
    def vars(seaiceOriDir,seaiceNewDir):
        """
        Generates a new WRF output file from scratch.
        
        Parameters
        ----------

        >>> seaseaiceBox = [lon_min, lon_max, lat_min, lat_max]

        >>> seaiceLevel = np.arange(min(seaiceLevel), max(seaiceLevel)+1)

        >>> seaiceTStep = np.arange(min(seaiceTStep), max(seaiceTStep)+1)
        """
        # kwargs
        seaiceBox = kwargs.get('seaiceBox')
        seaiceTStep = kwargs.get('seaiceTStep')        
        age = kwargs.get('age')  
        aice = kwargs.get('aice')  
        hice = kwargs.get('hice')  
        vice = kwargs.get('vice')  
        uice = kwargs.get('uice')  
        snowthick = kwargs.get('snowthick')  
        tisrf = kwargs.get('tisrf')  
        iomflx = kwargs.get('iomflx')  
        ti = kwargs.get('ti')  

        # Original output file.
        iceRawFile             = Dataset(iceOriDir, mode='r')
        iceNewFile             = Dataset(iceNewDir, 'w', format='NETCDF4')
        iceNewFile.title       = "Budgell Sea-ice output file"
        iceNewFile.description = "Created with Ekman Toolbox in " + time.ctime(time.time())
        iceNewFile.link        = "https://github.com/uesleisutil/Ekman"

        if seaiceBox is not None:
            def bbox2ij(lon,lat,seaiceBox=[-160., -155., 18., 23.]):
                """Return indices for i,j that will completely cover the specified bounding box.

                i0,i1,j0,j1 = bbox2ij(lon,lat,seaiceBox)
                
                lon,lat = 2D arrays that are the target of the subset
                seaiceBox = list containing the bounding box: [lon_min, lon_max, lat_min, lat_max]

                Example
                -------  
                >>> i0,i1,j0,j1 = bbox2ij(lon_rho,[-71, -63., 39., 46])
                >>> h_subset = nc.variables['h'][j0:j1,i0:i1]       
                """
                seaiceBox=np.array(seaiceBox)
                mypath=np.array([seaiceBox[[0,1,1,0]],seaiceBox[[2,2,3,3]]]).T
                p = path.Path(mypath)
                points = np.vstack((lon.flatten(),lat.flatten())).T
                n,m = np.shape(lon)
                inside = p.contains_points(points).reshape((n,m))
                ii,jj = np.meshgrid(range(m),range(n))
                return min(ii[inside]),max(ii[inside]),min(jj[inside]),max(jj[inside])

            lon_rho = iceRawFile.variables['lon_rho'][:,:]
            lat_rho = iceRawFile.variables['lat_rho'][:,:]
            i0,i1,j0,j1 = bbox2ij(lon_rho,lat_rho,seaiceBox)
            lon_rho = iceRawFile.variables['lon_rho'][j0:j1, i0:i1]
            lat_rho = iceRawFile.variables['lat_rho'][j0:j1, i0:i1]  
            iceNewFile.createDimension('eta_rho', len(lon_rho[:,0]))    
            iceNewFile.createDimension('xi_rho', len(lon_rho[0,:]))    
            print("Bounding box selected. New domain limits are: Longitude "+str(seaiceBox[0])+"/"+str(seaiceBox[1])+" and Latitude "+str(seaiceBox[2])+"/"+str(seaiceBox[3])+".")       
        else: 
            print("No bounding box selected: Using XLAT and XLONG variables from input file.")           
            lon_rho = iceNewFile.variables['lon_rho'][:,:]
            lat_rho = iceNewFile.variables['lat_rho'][:,:] 
            eta_rho = iceNewFile.dimensions['eta_rho']
            xi_rho  = iceNewFile.dimensions['xi_rho']
            iceNewFile.createDimension('eta_rho', len(eta_rho))    
            iceNewFile.createDimension('xi_rho', len(xi_rho))  
        
        iceNewLon = iceNewFile.createVariable('lon_rho', 'd', ('eta_rho', 'xi_rho'), fill_value=seaiceFillVal)
        iceNewLon.long_name = 'Longitude on RHO-points'
        iceNewLon.units = 'degree_east'
        iceNewLon.standard_name = 'longitude'
        iceNewLon[:,:] = lon_rho

        iceNewLat = iceNewFile.createVariable('lat_rho', 'd', ('eta_rho', 'xi_rho'), fill_value=seaiceFillVal)
        iceNewLat.long_name = 'Latitude on RHO-points'
        iceNewLat.units = 'degree_north'
        iceNewLat.standard_name = 'latitude'
        iceNewLat[:, :] = lat_rho    
     
        # New Sea-ice output file.
        iceNewFile.createDimension('ocean_time', 0)
        ice_time= iceRawFile.variables['ocean_time']
        iceNewOTdim = iceNewFile.createVariable('ocean_time', dtype('double').char, ('ocean_time'))
        iceNewOTdim.long_name = ice_time.units
        iceNewOTdim.units = ice_time.units

        s_rho = iceRawFile.dimensions['s_rho']
        s_w = iceRawFile.dimensions['s_w'] 
            
        if seaiceTStep == True:
            ntimes = seaiceTStep
            print("Time-step selected: Working from time-step "+str(np.argmin(ntimes))+" to "+str(np.argmax(ntimes))+".")
        else:
            ntimes = iceRawFile.variables['ocean_time'][:]
            print("No time-step selected. Working with entire time-step.")

        # If Budgell Sea-Ice has been chosen.                                               
        if age == True:
            print('Working on Budgell Sea-Ice Age.')
            bar = IncrementalBar(max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if seaiceBox == True:
                    if i == np.argmin(ntimes):
                        iceRawVar = iceRawFile.variables['ageice'][ntimes[0]+i,j0:j1, i0:i1]
                        iceNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                        iceNewVar = iceNewFile.createVariable('ageice', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=seaiceFillVal)
                        iceNewVar.long_name = 'Sea-ice age'
                        iceNewVar.units = 'S'
                        iceNewVar[i,:,:] = iceRawVar  
                    else:
                        iceRawFile = iceRawFile.variables['ageice'][ntimes[0]+i,j0:j1,i0:i1]
                        iceNewVar[i,:,:] = iceRawVar                         
                else: 
                    if i == np.argmin(ntimes):
                        iceRawVar = iceRawFile.variables['ageice'][ntimes[0]+i,:,:]
                        iceNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                        iceNewVar = iceNewFile.createVariable('ageice', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=seaiceFillVal)
                        iceNewVar.long_name = 'Sea-ice age'
                        iceNewVar.units = 'S'
                        iceNewVar[i,:,:] = iceRawVar  
                    else:
                        iceRawVar = iceRawFile.variables['ageice'][ntimes[0]+i,:,:]
                        iceNewVar[i,:,:] = iceRawVar                     
                bar.next()
            bar.finish() 

        # If Budgell Sea-Ice Fraction of Cell Covered by Ice has been chosen.                                               
        if aice == True:
            print('Working on Budgell Sea-Ice Fraction of Cell Covered by Ice.')
            bar = IncrementalBar(max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if seaiceBox == True:
                    if i == np.argmin(ntimes):
                        iceRawVar = iceRawFile.variables['aice'][ntimes[0]+ii,j0:j1, i0:i1]
                        iceNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                        iceNewVar = iceNewFile.createVariable('aice', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=seaiceFillVal)
                        iceNewVar.long_name = 'Fraction of Cell Covered by Ice'
                        iceNewVar[i,:,:] = iceRawVar  
                    else:
                        iceRawVar = iceRawFile.variables['aice'][ntimes[0]+i,j0:j1,i0:i1]
                        iceNewVar[i,:,:] = iceRawVar                         
                else: 
                    if i == np.argmin(ntimes):
                        iceRawVar = iceRawFile.variables['aice'][ntimes[0]+i,:,:]
                        iceNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                        iceNewVar = iceNewFile.createVariable('aice', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=seaiceFillVal)
                        iceNewVar.long_name = 'Fraction of Cell Covered by Ice'
                        iceNewVar[i,:,:] = iceRawVar  
                    else:
                        iceRawVar = iceRawFile.variables['aice'][ntimes[0]+i,:,:]
                        iceNewVar[i,:,:] = iceRawVar                     
                bar.next()
            bar.finish() 

        # If Budgell Sea-Ice Average Ice Thickness in Cell has been chosen.                                               
        if hice == True:
            print('Working on Budgell Sea-Ice Average Ice Thickness in Cell.')
            bar = IncrementalBar(max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if seaiceBox == True:
                    if i == np.argmin(ntimes):
                        iceRawVar = iceRawFile.variables['hice'][ntimes[0]+i,j0:j1, i0:i1]
                        iceNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                        iceNewVar = iceNewFile.createVariable('hice', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=seaiceFillVal)
                        iceNewVar.long_name = 'Average Ice Thickness in Cell'
                        iceNewVar.units     = 'meters'
                        iceNewVar[i,:,:] = iceRawVar  
                    else:
                        iceRawVar = iceRawFile.variables['hice'][ntimes[0]+i,j0:j1,i0:i1]
                        iceNewVar[i,:,:] = iceRawVar                         
                else: 
                    if i == np.argmin(ntimes):
                        iceRawVar = iceRawFile.variables['hice'][ntimes[0]+i,:,:]
                        iceNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                        iceNewVar = iceNewFile.createVariable('hice', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=seaiceFillVal)
                        iceNewVar.long_name = 'Average Ice Thickness in Cell'
                        iceNewVar.units     = 'meters'
                        iceNewVar[i,:,:] = iceRawVar  
                    else:
                        iceRawVar = iceRawFile.variables['hice'][ntimes[0]+i,:,:]
                        iceNewVar[i,:,:] = iceRawVar                     
                bar.next()
            bar.finish() 

        # If Budgell Sea-Ice V-Velocity has been chosen.                                               
        if vice == True:
            print('Working on Budgell Sea-Ice V-Velocity.')
            bar = IncrementalBar(max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if seaiceBox == True:
                    if i == np.argmin(ntimes):
                        iceRawVar = iceRawFile.variables['vice'][ntimes[0]+i,j0:j1, i0:i1]
                        iceNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                        iceNewVar = iceNewFile.createVariable('vice', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=seaiceFillVal)
                        iceNewVar.long_name = 'V-component of Ice Velocity'
                        iceNewVar.units     = 'm s-1'
                        iceNewVar[i,:,:] = iceRawVar  
                    else:
                        iceRawVar = iceRawFile.variables['vice'][ntimes[0]+i,j0:j1,i0:i1]
                        iceNewVar[i,:,:] = iceRawVar                         
                else: 
                    if i == np.argmin(ntimes):
                        iceRawVar = iceRawFile.variables['vice'][ntimes[0]+i,:,:]
                        iceNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                        iceNewVar = iceNewFile.createVariable('vice', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=seaiceFillVal)
                        iceNewVar.long_name = 'V-component of Ice Velocity'
                        iceNewVar.units     = 'm s-1'
                        iceNewVar[i,:,:] = iceRawVar  
                    else:
                        iceRawVar = iceRawFile.variables['vice'][ntimes[0]+i,:,:]
                        iceNewVar[i,:,:] = iceRawVar                     
                bar.next()
            bar.finish() 

        # If Budgell Sea-Ice U-Velocity has been chosen.                                               
        if uice == True:
            print('Working on Budgell Sea-Ice U-Velocity.')
            bar = IncrementalBar(max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if seaiceBox == True:
                    if i == np.argmin(ntimes):
                        iceRawVar = iceRawFile.variables['uice'][ntimes[0]+i,j0:j1, i0:i1]
                        iceNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                        iceNewVar = iceNewFile.createVariable('uice', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=seaiceFillVal)
                        iceNewVar.long_name = 'U-component of Ice Velocity'
                        iceNewVar.units     = 'm s-1'
                        iceNewVar[i,:,:] = iceRawVar  
                    else:
                        iceRawVar = iceRawFile.variables['uice'][ntimes[0]+i,j0:j1,i0:i1]
                        iceNewVar[i,:,:] = iceRawVar                         
                else: 
                    if i == np.argmin(ntimes):
                        iceRawVar = iceRawFile.variables['uice'][ntimes[0]+i,:,:]
                        iceNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                        iceNewVar = iceNewFile.createVariable('uice', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=seaiceFillVal)
                        iceNewVar.long_name = 'U-component of Ice Velocity'
                        iceNewVar.units     = 'm s-1'
                        iceNewVar[i,:,:] = iceRawVar  
                    else:
                        iceRawVar = iceRawFile.variables['uice'][ntimes[0]+i,:,:]
                        iceNewVar[i,:,:] = iceRawVar                     
                bar.next()
            bar.finish() 

        # If Budgell Sea-cover Thickness has been chosen.                                               
        if snowthick == True:
            print('Working on Budgell Sea-cover Thickness.')
            bar = IncrementalBar(max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if seaiceBox == True:
                    if i == np.argmin(ntimes):
                        iceRawVar = iceRawFile.variables['snow_thick'][ntimes[0]+i,j0:j1, i0:i1]
                        iceNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                        iceNewVar = iceNewFile.createVariable('snow_thick', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=seaiceFillVal)
                        iceNewVar.long_name = 'Sea-cover Thickness'
                        iceNewVar.units     = 'meter'
                        iceNewVar[i,:,:] = iceRawVar  
                    else:
                        iceRawVar = iceRawFile.variables['snow_thick'][ntimes[0]+i,j0:j1,i0:i1]
                        iceNewVar[i,:,:] = iceRawVar                         
                else: 
                    if i == np.argmin(ntimes):
                        iceRawVar = iceRawFile.variables['snow_thick'][ntimes[0]+i,:,:]
                        iceNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                        iceNewVar = iceNewFile.createVariable('snow_thick', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=seaiceFillVal)
                        iceNewVar.long_name = 'Sea-cover Thickness'
                        iceNewVar.units     = 'meter'
                        iceNewVar[i,:,:] = iceRawVar  
                    else:
                        iceRawVar = iceRawFile.variables['snow_thick'][ntimes[0]+i,:,:]
                        iceNewVar[i,:,:] = iceRawVar                     
                bar.next()
            bar.finish() 

        # If Budgell Sea-ice Surface Temperature has been chosen.                                               
        if tisrf == True:
            print('Working on Budgell Sea-ice Surface Temperature.')
            bar = IncrementalBar(max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if seaiceBox == True:
                    if i == np.argmin(ntimes):
                        iceRawVar = iceRawFile.variables['tisrf'][ntimes[0]+i,j0:j1, i0:i1]
                        iceNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                        iceNewVar = iceNewFile.createVariable('tisrf', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=seaiceFillVal)
                        iceNewVar.long_name = 'Sea-Ice Surface Temperature'
                        iceNewVar.units     = 'Degree Celsius'
                        iceNewVar[i,:,:] = iceRawVar  
                    else:
                        iceRawVar = iceRawFile.variables['tisrf'][ntimes[0]+i,j0:j1,i0:i1]
                        iceNewVar[i,:,:] = iceRawVar                         
                else: 
                    if i == np.argmin(ntimes):
                        iceRawVar = iceRawFile.variables['tisrf'][ntimes[0]+i,:,:]
                        iceNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                        iceNewVar = iceNewFile.createVariable('tisrf', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=seaiceFillVal)
                        iceNewVar.long_name = 'Sea-Ice Surface Temperature'
                        iceNewVar.units     = 'Degree Celsius'
                        iceNewVar[i,:,:] = iceRawVar  
                    else:
                        iceRawVar = iceRawFile.variables['tisrf'][ntimes[0]+i,:,:]
                        iceNewVar[i,:,:] = iceRawVar                     
                bar.next()
            bar.finish() 

        # If Budgell Ice-Ocean Mass Flux has been chosen.                                               
        if iomflx == True:
            print('Working on Budgell Ice-Ocean Mass Flux.')
            bar = IncrementalBar(max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if seaiceBox == True:
                    if i == np.argmin(ntimes):
                        iceRawVar = iceRawFile.variables['iomflx'][ntimes[0]+i,j0:j1, i0:i1]
                        iceNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                        iceNewVar = iceNewFile.createVariable('iomflx', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=seaiceFillVal)
                        iceNewVar.long_name = 'Ice-Ocean Mass Flux'
                        iceNewVar.units     = 'm s-1'
                        iceNewVar[i,:,:] = iceRawVar  
                    else:
                        iceRawVar = iceRawFile.variables['iomflx'][ntimes[0]+i,j0:j1,i0:i1]
                        iceNewVar[i,:,:] = iceRawVar                         
                else: 
                    if i == np.argmin(ntimes):
                        iceRawVar = iceRawFile.variables['iomflx'][ntimes[0]+i,:,:]
                        iceNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                        iceNewVar = iceNewFile.createVariable('iomflx', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=seaiceFillVal)
                        iceNewVar.long_name = 'Ice-Ocean Mass Flux'
                        iceNewVar.units     = 'm s-1'
                        iceNewVar[i,:,:] = iceRawVar  
                    else:
                        iceRawFile = iceRawFile.variables['iomflx'][ntimes[0]+i,:,:]
                        iceNewVar[i,:,:] = iceRawVar                     
                bar.next()
            bar.finish() 

        # If Budgell Interior Ice Temperature has been chosen.                                               
        if ti == True:
            print('Working on Budgell Interior Ice Temperature.')
            bar = IncrementalBar(max=len(ntimes))
            for i in range(np.argmin(ntimes),len(ntimes),1):
                if seaiceBox == True:
                    if i == np.argmin(ntimes):
                        iceRawVar = iceRawFile.variables['ti'][ntimes[0]+i,j0:j1, i0:i1]
                        iceNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                        iceNewVar = iceNewFile.createVariable('ti', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=seaiceFillVal)
                        iceNewVar.long_name = 'Interior Ice Temperature'
                        iceNewVar.units     = 'Degree Celcius'
                        iceNewVar[i,:,:] = iceRawVar  
                    else:
                        iceRawVar = iceRawFile.variables['ti'][ntimes[0]+i,j0:j1,i0:i1]
                        iceNewVar[i,:,:] = iceRawVar                         
                else: 
                    if i == np.argmin(ntimes):
                        iceRawVar = iceRawFile.variables['ti'][ntimes[0]+i,:,:]
                        iceNewVar = np.zeros([len(ntimes),len(lat_rho), len(lon_rho)])
                        iceNewVar = iceNewFile.createVariable('ti', 'f', ('ocean_time', 'eta_rho', 'xi_rho'), fill_value=seaiceFillVal)
                        iceNewVar.long_name = 'Interior Ice Temperature'
                        iceNewVar.units     = 'Degree Celcius'
                        iceNewVar[i,:,:] = iceRawVar  
                    else:
                        iceRawVar = iceRawFile.variables['ti'][ntimes[0]+i,:,:]
                        iceNewVar[i,:,:] = iceRawVar                     
                bar.next()
            bar.finish() 