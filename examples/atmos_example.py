# Import libraries.
import os
import numpy as np
from ekman.conf import createFolders
from ekman.atmos import atmos

projectName = 'Test_01'
projectAuthor = 'Ueslei Adriano Sutil'
selectWrfVars = True

def run():
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
  
    if selectWrfVars == True:  
        print('------------------------')        
        print('WRF section initialized')
        print('------------------------')  
        wrfOriginalFilename = "wrf.nc"
        wrfNewFilename = "wrf_new.nc"
        wrfBox = [-53, -40, -32, -23] # [lon_min, lon_max, lat_min, lat_max]
        wrfLevel = np.arange(29, 29+1)
        wrfTStep = np.arange(159, 160+1)

        wrfOriDir = getFolder+'/'+projectName+'/'+wrfOriginalFilename
        wrfNewDir = getFolder+'/'+projectName+'/'+wrfNewFilename
        atmos.vars(wrfOriDir,wrfNewDir,wrfBox=wrfBox,wrfLevel=wrfLevel,latent=True)

run()