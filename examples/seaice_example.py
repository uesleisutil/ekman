import numpy as np

iceOriginalFilename = 'ocean_his.nc'
iceNewFilename = 'sea_ice_his.nc'

selectIceBox = True
iceBox = [-53, -40, -32, -23]

selectIceTimeStep   = True
iceTimeStep = np.arange(200, 255+1)

selectIceVars = False
iceAge = False
iceA = False
iceH = False
iceV = False
iceU = False
iceSnowThick = False
iceSurfaceTemp = False
iceOceanMassFlux = False
iceInteriorTemp = False
