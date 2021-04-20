.. _GettingStarted

*******************
Getting Started with Dorado
*******************

Install
-------


Install dorado by running the command: ``pip install dorado``.

    .. note:: The `dorado <https://pypi.org/project/dorado/>`_ project on `PyPI <https://pypi.org/>`_ . 

Pip should install all relevent dependencies for dorado into the python environment 
pip was run from.

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


Getting to know Clippy
======================

You may remember that condescending paperclip from the 90s that helped a computer iliterate
world cope with CRT displays, cleaning mouse balls, and the blue screen of death; well now 
the ghost of clippy is back to help you manage and organized and effective local astronomical
file system. 

.. image:: /dorado/_static/clippy.jfif