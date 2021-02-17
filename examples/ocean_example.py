import numpy as np

romsOriginalFilename = 'ocean_his.nc'
romsNewFilename = 'ocean_his_new.nc'

selectRomsBox = True
romsBox = [-53, -40, -32, -23] # [lon_min, lon_max, lat_min, lat_max]

selectRomsLevel = True
romsLevel = np.arange(29, 29+1)

selectRomsTimeStep   = True
romsTimeStep = np.arange(159, 160+1)

selectRomsVars = True
romsSST = True
romsTemp = False
romsSalt = False
romsTKE = False
romsRho = False
romsZeta = False
romsW = False
romsOmega = False
romsLatent = False
romsSensible = False
romsLWRad = False
romsSWRad = False
romsEminusP = False
romsEvaporation = False
romsUwind = False
romsVwind = False
romsU = False
romsV = False
romsUbar = False
romsVbar = False