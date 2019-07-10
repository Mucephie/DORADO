# DRACO

Data Reduction at the Allan I. Carswell Observatory (DRACO) is a python based expansion of Astropy(and affiliated packages) that aims to be a simple and common framework for data reduction tailored for life at the Allan I. Carswell Observatory at York university, Toronto, Ontario, Canada

## draco

---

The `draco.py` file is DRACO's core and holds the main functionality needed to make a reduction script. To make a reduction script first write:  

> `import draco`

in the header of your reduction script. Next you will most likely need to import a few other tools that DRACO is founded upon.

> `from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table`

DRACO uses `astropy.io.fits` and `astropy.table` to handle alot of the heavy lifting when it comes to opening, reading, writing, and exporting of common astronomical datatypes such as `.fits` files and `.csv` files. Apart from [astropy](https://www.astropy.org/index.html)'s offerings; DRACO also uses [numpy](http://www.numpy.org/) to store and operate on large series of datasets at the same time and to write the bulk of DRACO's logic.

The import of [matplotlib](https://matplotlib.org/) may be optional for some users if they wish to use DRACO's built in plots (which utilize matplotlib). For users who wish to create their own specialized plots in their reduction scripts, importing of matplotlib is required.

Internally, DRACO uses [astroscappy](https://github.com/astropy/astroscrappy), an astropy affiliated package, for the removal of cosmic rays when opening a `.fits` image. Future versions of DRACO aim to make this feature optional and to reduce the time it takes to search for and remove cosmic rays from large series of data.

## Dependencies

---

Draco aims to be as lightweight as possible and utilize as little dependencies as it can. In the future DRACO will be updated so that some dependencies, such as astroscrappy, are optional and only needed for the functions in which they are utilized.

Currently DRACO relies on:  
1. [numpy](http://www.numpy.org/)
1. [matplotlib](https://matplotlib.org/)
1. [astropy](https://www.astropy.org/index.html)
1. [astroscappy](https://github.com/astropy/astroscrappy)

Future planned updates to DRACO may include dependencies such as:  
1. [CCDPROC](https://ccdproc.readthedocs.io/en/latest/index.html#)
1. [Ginga](https://ejeschke.github.io/ginga/)
1. [imexam](https://imexam.readthedocs.io/en/latest/index.html)
1. [Photutils](https://photutils.readthedocs.io/en/stable/index.html)
