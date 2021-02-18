# Ekman
[![Build Status](https://travis-ci.org/uesleisutil/ekman.svg?branch=main)](https://travis-ci.org/uesleisutil/ekman)
[![MIT Licence](https://badges.frapsoft.com/os/mit/mit.png?v=103)](https://opensource.org/licenses/mit-license.php)
[![FOSSA Status](https://app.fossa.com/api/projects/git%2Bgithub.com%2Fuesleisutil%2FEkman.svg?type=shield)](https://app.fossa.com/projects/git%2Bgithub.com%2Fuesleisutil%2FEkman?ref=badge_shield)
[![codecov](https://codecov.io/gh/uesleisutil/ekman/branch/main/graph/badge.svg?token=J1xWME34Fq)](https://codecov.io/gh/uesleisutil/ekman)

Ekman is an open-source Python package to postprocess, data analyze and visualize Weather Research Forecast (WRF), Regional Ocean Modeling System (ROMS) and Budgell's Sea Ice outputs files.

## Installation

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install.

```bash
pip install ekman
```

> **Note**: `Numpy`, `netCDF4`, `wrf-python` and `matplotlib` are required to run this program.


## Usage

```python
import ekman

from ekman.conf import createFolders
from ekman.ocean import ocean
from ekman.atmos import atmos
from ekman.ice import ice
```

## Examples
Check the [Examples] folder.


### Ocean module

The oceanic module currently supports the [ROMS] model. The arguments accepted by the library are:

| **kwargs  | Details |
| ------    | ------ |
| romsBox   | Crop the grid into a smaller one. Similar to the `sellonlatbox` function from [CDO]. Use as: `[lon_min, lon_max, lat_min, lat_max]` |
| romsLevel | Select vertical levels. Use as: `np.arange(min(romsLevel), max(romsLevel)+1)` |
| romsTStep | Select vertical step: Use as: ` np.arange(min(romsTStep), max(romsTStep)+1)` |
| temp      | Potential Temperature (DegC) |
| salt      | Salinity (PSU) |
| tke       | Turbulent Kinectic Energ (m2 s-2)|
| rho       | Density Anomaly (kg m-3) |
| w         | Vertical Momentum Component (m s-1) |
| omega     | S-coordinate Vertical Momentum Component (m s-1) |
| zeta      | Density Anomaly (kg m-3) |
| latent    | Latent Heat Flux (W m-2) |
| sensible  | Sensible Heat Flux (W m-2) |
| lwrad     | Net Longwave Radiation Flux (W m-2) |
| swrad     | Net Shortwave Radiation Flux (W m-2) |
| eminusp   | Bulk Flux Surface Net Freshwater Flux (m s-1) |
| evapo     | Evaporation Rate (kg m-2 s-1') |
| uwind     | Surface U-wind Component on Mass points (m s-1) |
| vwind     | Surface V-wind Component on Mass points (m s-1) |
| u         | U-wind Component on U points (m s-1) |
| v         | V-wind Component on V points (m s-1) |
| ubar      | Vertically Integrated U-momentum Component on V points (m s-1) |
| vbar      | Vertically Integrated V-momentum Component on V points (m s-1) |

### Atmosphere module

The atmospheric module currently supports the [WRF] model. The arguments accepted by the library are:

| **kwargs  | Details |
| ------    | ------ |
| wrfBox   | Crop the grid into a smaller one. Similar to the `sellonlatbox` function from [CDO] . Use as: `[lon_min, lon_max, lat_min, lat_max]`  |
| wrfLevel | Select vertical levels. Use as: `np.arange(min(romsLevel), max(romsLevel)+1)` |
| wrfTStep | Select vertical step: Use as: ` np.arange(min(romsTStep), max(romsTStep)+1)` |
| temp     | Temperature (DegC) |
| rh       | Potential Temperature (DegC) |
| td       | Relative Humidity (%) |
| twb      | Wet Bulb Temperature (DegC) |
| tv       | Virtual Temperature (DegC) |
| pres     | Full Model Pressure (hPa)|
| avo      | Absolute Vorticity (10-5 s-1) |
| pvo      | Potential Vorticity (10-5 s-1) |
| dbz      | Reflectivity (dBZ) |
| geopt    | Geopotential Height for the Mass Grid (m2 s-2) |
| omega    | Omega (Pa s-1) |
| unstagu  | U-component of Wind on Mass Points (m s-1) |
| unstagv  | V-component of Wind on Mass Points (m s-1) |
| unstagw  | W-component of Wind on Mass Points (m s-1) |
| uvmet    | U and V Components of Wind Rotated to Earth Coordinates (m s-1) |
| uvmet10m | 10 m U and V Components of Wind Rotated to Earth Coordinates (m s-1) |
| latent   | Latent Heat Flux (W m-2) |
| sensible | Sensible Heat Flux (W m-2) |
| slp      | Sea Level Pressure (hPa) |
| rh2      | 2m Relative Humidity (%) |
| terrain  | Model Terrain Height (m) |
| td2      | 2m Dew Point Temperature (DegC) |
| landmask | Land mask (1=land, 0=water) |

### Sea-ice module

The sea-ice module currently supports the [Paul Budgell's] sea-ice model. The arguments accepted by the library are:

| **kwargs  | Details |
| ------    | ------ |
| seaiceBox   | Crop the grid into a smaller one. Similar to the `sellonlatbox` function from [CDO] . Use as: `[lon_min, lon_max, lat_min, lat_max]`  |
| seaiceTStep | Select vertical step: Use as: ` np.arange(min(romsTStep), max(romsTStep)+1)` |
| age     | Sea-ice age (s) |
| aice       | Sea-Ice Fraction of Cell Covered by Ice |
| hice       | Average Ice Thickness in Cell (m) |
| vice      | V-component of Ice Velocity (m s-1) |
| uice | U-component of Ice Velocity (m s-1) |
| snowthick | Sea-cover Thickness (m)|
| tisrf      | Sea-Ice Surface Temperature (DegC) |
| iomflx      |  Ice-Ocean Mass Flux (m s-1) |
| ti      |  Interior Ice Temperature (DegC) |

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

### Next upcoming updates:
 - Find a way to loop over a list of strings to delete all the crappy infinite command lines in `ocean.py`, `atmos.py` and `ice.py`.
 - Add a wave module for SWAN and WaveWatch III models.
 - Add staggered U, V and W wind components in `atmos.py`.

### Know your bugs 


## License
The code in this project is released under the [MIT](https://choosealicense.com/licenses/mit/) license.

[![FOSSA Status](https://app.fossa.com/api/projects/git%2Bgithub.com%2Fuesleisutil%2FEkman.svg?type=large)](https://app.fossa.com/projects/git%2Bgithub.com%2Fuesleisutil%2FEkman?ref=badge_large)

[//]: # (http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)

[Examples]: <https://github.com/uesleisutil/ekman/tree/main/examples>
[ROMS]: <https://www.myroms.org/>
[CDO]: <https://code.mpimet.mpg.de/projects/cdo/>
[Paul Budgell's]: <https://link.springer.com/article/10.1007/s10236-005-0008-3>
[WRF]: <https://www.mmm.ucar.edu/weather-research-and-forecasting-model> 
