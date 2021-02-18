"""
File name: run_all_tests.py
Author: Ueslei Adriano Sutil
Created: 07 Apr 2020
Last modified:  03 Jan 2021

"""
# Import libraries.
import os
import numpy as np
from ekman.conf import createFolders
import numpy as np
from ekman.ocean import ocean
from ekman.atmos import atmos
from ekman.seaice import seaice
import pytest

projectName = 'Test_01'
projectAuthor = 'Ueslei Adriano Sutil'
selectRomsVars = True
selectWrfVars = False
selectSeaIceVars = False

def test_run():
    # Check if project already exists. If not, create a new folder.
    checkFolder = os.path.isdir(projectName) 
    getFolder   = os.getcwd()
    print('Initializing project: '+projectName)

    if checkFolder == False:
        createFolders(projectName)
        print('Project folders created. Now store raw output files in: '+getFolder+'/'+projectName)
        raise SystemExit(0)
    else:
        print('Project folder already created. Continuing.')
        pass

    if selectRomsVars == True:
        romsOriginalFilename = "roms.nc"
        romsNewFilename = "roms_new.nc"
        romsBox = [-53, -40, -32, -23] # [lon_min, lon_max, lat_min, lat_max]
        romsLevel = np.arange(1,1+1)
        romsTStep = np.arange(159, 160+1)

        romsOriDir = getFolder+'/'+projectName+'/'+romsOriginalFilename
        romsNewDir = getFolder+'/'+projectName+'/'+romsNewFilename
        ocean.vars(romsOriDir,romsNewDir,romsBox=romsBox,romsLevel=romsLevel,romsTStep=romsTStep,temp=True)

    if selectWrfVars == True:  
        wrfOriginalFilename = "wrf.nc"
        wrfNewFilename = "wrf_new.nc"
        wrfBox = [-53, -40, -32, -23] # [lon_min, lon_max, lat_min, lat_max]
        wrfLevel = np.arange(29, 29+1)
        wrfTStep = np.arange(99, 100+1)

        wrfOriDir = getFolder+'/'+projectName+'/'+wrfOriginalFilename
        wrfNewDir = getFolder+'/'+projectName+'/'+wrfNewFilename
        atmos.vars(wrfOriDir,wrfNewDir,wrfBox=wrfBox,wrfLevel=wrfLevel,wrfTStep=wrfTStep,latent=True,sensible=True)

    if selectSeaIceVars == True:
        seaiceOriginalFilename = "roms.nc"
        seaiceNewFilename = "roms_seaice_new.nc"
        seaiceBox = [-53, -40, -32, -23] # [lon_min, lon_max, lat_min, lat_max]
        seaiceLevel = np.arange(29, 29+1)
        seaiceTStep = np.arange(159, 160+1)
        
        seaiceOriDir = getFolder+'/projects/'+projectName+'/'+seaiceOriginalFilename
        seaiceNewDir = getFolder+'/projects/'+projectName+'/'+seaiceNewFilename
        seaice.vars(seaiceOriDir,seaiceNewDir,seaiceBox,seaiceLevel,seaiceTStep)

    # Program finished.
    print('Program finished.')

test_run()