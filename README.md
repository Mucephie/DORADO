# DRACO

Data Reduction at the Allan I. Carswell Observatory (DRACO) is a python based expansion of astropy(and affiliated packages) that aims to be a simple and common framework for data reduction tailored for life at the Allan I. Carswell Observatory at York university, Toronto, Ontario, Canada

### ☣⚠ Note: DRACO is in alpha phase ⚠☣

DRACO is in alpha phase and as such is highly experimental and still a work in progress. While still in active development DRACO may still be used and experimented with; but a stable channel will not be created till DRACO enters a beta phase in the future.

Use of DRACO in the alpha phase will not result in any risk to your computer or hardware, but may affect stability of software that are written with DRACO or the quality of data that is produced using DRACO.

During the alpha phase, documentation and code commenting may be subpar as many things will being changing regularly with new code being written and rewritten and in-line testing appearing temporarily.

## draco

The `dracoOP2.py` file is DRACO's core and holds the main functionality needed to make a reduction script. Currently to install DRACO, you will need to download the `dracoOP2.py` file into your project folder/path. In the future DRACO will be released as a PIP or conda installable package.

To make a reduction script first write:  

```python
import dracoOP2 as draco
```

into the header/imports of your reduction script. Next you will most likely need to import a few other tools that DRACO is founded upon.

```python
 from astropy.io import fits  
 import numpy as np  
 import matplotlib.pyplot as plt  
 from astropy.table import Table
```

DRACO uses `astropy.io.fits` and `astropy.table` to handle a lot of the heavy lifting when it comes to opening, reading, writing, and exporting of common astronomical data types such as `.fits` files and `.csv` files. Apart from [astropy](https://www.astropy.org/index.html)'s offerings; DRACO also uses [numpy](http://www.numpy.org/) to store and operate on large series of data sets at the same time and to write the bulk of DRACO's logic.

The import of [matplotlib](https://matplotlib.org/) may be optional for some users if they wish to use DRACO's built in plots (which utilize matplotlib). For users who wish to create their own specialized plots in their reduction scripts, importing of matplotlib is required.


## Dependencies

Draco aims to be as lightweight as possible and utilize as little dependencies as it can. In the future DRACO will be updated so that some dependencies, such as astroscrappy, are optional and only needed for the functions in which they are utilized.

Currently DRACO relies on:  
1. [numpy](http://www.numpy.org/)
1. [matplotlib](https://matplotlib.org/)
1. [astropy](https://www.astropy.org/index.html)
1. [CCDPROC](https://ccdproc.readthedocs.io/en/latest/index.html#)
1. [scipy](https://www.scipy.org/)
1. [scikit-image](https://scikit-image.org/)

Future planned updates to DRACO may include dependencies such as:  

1. [Ginga](https://ejeschke.github.io/ginga/)
1. [imexam](https://imexam.readthedocs.io/en/latest/index.html)
1. [Photutils](https://photutils.readthedocs.io/en/stable/index.html)

## Variable star research and Timeseries analysis

DRACO's current focus is the reduction of a raw observing session into an organized timeseries data object which can be further configured and analyzed by DRACO to produce lightcurves and power spectral density (PSD) plots for Fourier based frequency analysis. DRACO has the ability to process target pixel files (TPF's) of variable stars from the [transiting exoplanet survey satellite](https://tess.mit.edu/) (TESS) and construct and 'Fourier fit' lightcurves from the TPF produced timeseries. 

This capability is currently only meant for internal use as it is very naorrow and experimental by nature; and not yet ready for scientific use.

## Contact

If you have any questions or would like to contribute to DRACO, please contact [@mucephie](https://github.com/Mucephie) at <mucephie@my.yorku.ca>. If you are interested in the Allan I. Carswell observatory at York university, you can find more information at our [website](http://observatory.info.yorku.ca/).  

---

<p align="justify">
🔥🌈🎇✨🔭💻💾💽📷📡📺📓📚🔎📀🚀🛰🛸🌌🪐🌎🏳‍🌈🌒☄💫🕳💬☢🔥 
</p>
