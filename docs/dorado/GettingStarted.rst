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

    .. note:: CCDPROC requires the ``astroscrappy`` package for install. Currently, astroscrappy requires 
                Visual C++ redistributable for build and function. This dependency is currently having issues 
                with ``python 3.9``, compatability with python 3.9 is incoming. The install of astroscrappy may also 
                require the ``wheel`` python package for some python environments. If you are having trouble 
                installing DORADO, please contact us for assistance.

Babies first Dorado script
==========================

Adding Dorado into your flow is fairly simple and can be broken into 3 main steps; 
import, initialize, and call. 

.. code:: python

        ## import
        # import dorado utilities
        from dorado.filer import Filer
        from dorado.target import Target

        ## initialize
        # make an instance of clippy
        clippy = Filer()
        # create a target object 'target of interest' toi
        toi = Target('target_name')
        # create a control target object
        control = Target('control_target_name')
        # create a series object
        ceres = clippy.mkceres('2021-01-01+02', target = toi)

        ## call
        # calibrate the series red filter data
        ceres.calibrate('R')
        # align the series red filter data
        ceres.align('R', clippy, cache = True)
        # perform differential photometry on the target in the red filter 
        # data using the control as a reference
        ceres.dorphotc('R', toi, control)

        ## finish by saving
        # save the resulting data
        clippy.savewrk(ceres)
        # record the results
        toi.record(clippy, ceres)


Next: :doc:`Filer</dorado/filer>`
