# Import libraries.
import os
from ekman.conf import createFolders
import numpy as np
from ekman.seaice import seaice

projectName = 'Test_01'
projectAuthor = 'Your name'
selectSeaIceVars = True

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

    if selectSeaIceVars == True:
        print('---------------------------')  
        print('Sea-Ice section initialized')
        print('---------------------------')  
        seaiceOriDir = getFolder+'/projects/'+projectName+'/'+iceOriginalFilename
        seaiceNewDir = getFolder+'/projects/'+projectName+'/'+iceNewFilename
        seaice(seaiceOriDir,seaiceNewDir)
    # Program finished.
    print('Program finished.')

run()