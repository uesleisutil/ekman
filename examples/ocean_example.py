# Import libraries.
import os
import numpy as np
from ekman.conf import createFolders
from ekman.ocean import ocean

projectName = 'Test_01'
projectAuthor = 'Your name'
selectRomsVars = True

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

    if selectRomsVars == True:

        print('------------------------')        
        print('ROMS section initialized')
        print('------------------------')  
        romsOriginalFilename = "roms.nc"
        romsNewFilename = "roms_new.nc"
        romsBox = [-53, -40, -32, -23] # [lon_min, lon_max, lat_min, lat_max]
        romsLevel = np.arange(29, 29+1)
        romsTStep = np.arange(159, 160+1)

        romsOriDir = getFolder+'/'+projectName+'/'+romsOriginalFilename
        romsNewDir = getFolder+'/'+projectName+'/'+romsNewFilename
        ocean.vars(romsOriDir,romsNewDir,romsBox=romsBox,romsLevel=romsLevel,zeta=False,salt=True,temp=True)

    # Program finished.
    print('Program finished.')

run()