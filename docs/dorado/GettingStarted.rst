.. _GettingStarted

*******************
Getting Started with Dorado
*******************

Install
-------


Install dorado by running the command: ``pip install dorado``.

    .. note:: The `dorado <https://pypi.org/project/dorado/>`_ project on `PyPI <https://pypi.org/>`_ . 

Pip should install all relevent dependencies for dorado into the python environment 
pip was run from. The dorado team recommends using the latest version of python unless otherwise 
mentioned for the best compatibility.

Dependencies
-------------

Dorado aims to be as lightweight as possible and utilize as little dependencies as it can. 

Currently DORADO relies on:  

1.  `numpy <http://www.numpy.org/>`_

2.  `matplotlib <https://matplotlib.org/>`_

3.  `astropy <https://www.astropy.org/>`_

4.  `CCDPROC <https://ccdproc.readthedocs.io/>`_

5.  `photutils <https://photutils.readthedocs.io/>`_

6.  `astroalign <https://astroalign.readthedocs.io/>`_

7.  `astroquery <https://astroquery.readthedocs.io/>`_



Babies first Dorado script
==========================

Adding Dorado into your flow is fairly simple and can be broken into 3 main steps; 
import, initialize, and call. 

.. code:: python

        ## import
        # import dorado utilities
        from dorado.clippy import Clippy
        from dorado.zellars import Zellars

        ## initialize
        # make an instance of clippy
        clip = Clippy()
        # create a target object
        target = Zellars('target_name')
        # create a control target object
        # create a series object
        ceres = clip.mkceres('2021-01-01+02', target = target)

        ## call
        # calibrate the series red filter data
        ceres.calibrate('R')
        # align the series red filter data
        ceres.align('R', clip, cache = True)
        # perform differential photometry on the target in the red filter 
        # data using the control as a reference
        ceres.dorphotc('R', target, control, shape = 21)

        ## finish by saving
        # save the resulting data
        clip.savewrk(ceres)
        # record the results
        target.record(clip, ceres)


Next: :doc:`Clippy</dorado/clippy>`
