# DRACO

Data Reduction at the Allan I. Carswell Observatory (DRACO) is a python based expansion of astropy(and affiliated packages) that aims to be a simple and common framework for data reduction tailored for life at the Allan I. Carswell Observatory at York university, Toronto, Ontario, Canada

### ☣⚠ Note: DRACO is in alpha phase ⚠☣

DRACO is in alpha phase and as such is highly experimental and still a work in progress. While still in active development DRACO may still be used and experimented with; but a stable channel will not be created till DRACO enters a beta phase in the future.

Use of DRACO in the alpha phase will not result in any risk to your computer or hardware, but may affect stability of software that are written with DRACO or the quality of data that is produced using DRACO.

During the alpha phase, documentation and code commenting may be subpar as many things will being changing regularly with new code being written and rewritten and in-line testing appearing temporarily.

## draco

The `draco.py` file is DRACO's core and holds the main functionality needed to make a reduction script. Currently to install DRACO, you will need to download the `draco.py` file into your project folder/path. In the future DRACO will be released as a PIP or conda installable package.

To make a reduction script first write:  

```python
import draco
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

Internally, DRACO uses [astroscappy](https://github.com/astropy/astroscrappy), an astropy affiliated package, for the removal of cosmic rays when opening a `.fits` image. Future versions of DRACO aim to make this feature optional and to reduce the time it takes to search for and remove cosmic rays from large series of data.

## Dependencies

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

## im_reduc

The file `im_reduc.py` is an example of how simple it is to make a tailored reduction script for any task. `im_reduc.py` utilizes a very basic image reduction procedure to remove cosmic rays, bias, sky bias, and flat fields from series taken in multiple filter bands to produce series of fully reduced astronomical images. These images can then be opened in [DS9](http://ds9.si.edu/site/Home.html), [Ginga](https://ejeschke.github.io/ginga/), or plotted through DRACO/[matplotlib](https://matplotlib.org/) (see example below).  

![Draco image plot example; Comet 46P Wirtanen.](/assets/images/im_reduc_example.png)

> _DRACO/im_reduc example output using DRACO/matplotlib display and test data taken at the Allan I. Carswell observatory on the 60cm telescope in December 2018 of comet 46P 'Wirtanen'_  

`im_reduc.py` also has the capability to export the reduced frames(and master images) as individual `.fits` files to the directory of choice. These exported `.fits` files can then be opened by other astronomical data systems such as the depreciated [IRAF](https://iraf-community.github.io/) or digital graphics software such as [GIMP](https://www.gimp.org/).

## Contact

If you have any questions or would like to contribute to DRACO, please contact [@mucephie](https://github.com/Mucephie) at <mucephie@my.yorku.ca>. If you are interested in the Allan I. Carswell observatory at York university, you can find more information at our [website](http://observatory.info.yorku.ca/).  

---

<p align="justify">
🔥🌈🎇✨🔭💻💾💽📷📡📺📓📚🔎📀🚀🛰🛸🌌🪐🌎🏳‍🌈🌒☄💫🕳💬☢🔥 
</p>
