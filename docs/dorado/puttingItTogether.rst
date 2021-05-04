
Using Dorado In The Wild
========================


Lets imagine that we just did an observing run and saved our results under the date folder ``.dorado/data/2021-01-02+03/``, we took 
flatfields for ``R``, ``V``, and ``B`` filters but forgot to obtain bias frames. During the night we observed a target ``McCool's star``
and a nearby star ``Boring star`` for a few hours with many 1-minute exposures in the red filter, as well as a few visible and blue
frames at the end of the night. We would now like to turn that
raw data into a lightcurve of McCool's star using Boring star as our control while performing differential photometry. We've set up
our data to transfer from our telescopes control computer in the observatory control room to our data computer in the main building
storing it in the aformentioned folder when the night is done. We've also set things up so that when the data is done transfering 
a script is called to process our fresh data and neatly store the results for us to enjoy after we go home for some much needed rest.

Let's write that script now so we can get a feel for using Dorado in the wild. The first thing we want to do is import our Dorado
utilities (including Clippy) and then initialize Clippy and construct our Zellars(target) objects.

.. code:: python

        # import Clippy from Dorado
        from dorado.clippy import Clippy
        from dorado.zellars import Zellars

        ## initialize
        # make an instance of clippy
        clip = Clippy()
        # create a target object
        target = Zellars('McCool's star')
        # create a control target object
        target = Zellars('Boring star')

Our next step is to tell Clippy that you want to find and read in last nightss data:

.. code:: python

        # retrieve the previous nights datestring
        night = clip.get_night()
        # create a series object from last nights data
        ceres = clip.mkceres(night, target = target)

Clippy has now read through the directory and scanned for each image type. Calibration frames are stacked
into base calibration frames, which are passed along with the science frames to into a data series class
known as Ceres. The Ceres object is returned to us by Clippy and now we wish to calibrate and align the 
data. This is done via the ``calibrate()`` and ``align()`` Ceres commands. Since there is data for more
than one filter stored in our Ceres object, we need to specify that we are interested in performing these
actions on the red filter.

.. code:: python

        ## call
        # calibrate the series red filter data
        ceres.calibrate('R')
        # align the series red filter data
        ceres.align('R', clip, cache = True)

In the alignment step Dorado called 'Astrometry.Net' to plate solve the aligned field so we have WCS
coordinates for this filters stack of data. This will allow Dorado to quickly find McCool's star and
Boring star in the images and perform photometry on them. This is done using the ``dorphotc()` Ceres
command, again we are specifying that we are interested in the red filter and we are feeding our target
and control objects, along with a PSF fit shape to the command.

.. code:: python

        # perform differential photometry on the target in the red filter 
        # data using the control as a reference
        ceres.dorphotc('R', target, control, shape = 21)

We now have our timeseries data for our lightcurve and wish to save it for review later. Lets save both
our target data and the calibrated/aligned images.

.. code:: python

        ## finish by saving
        # save the resulting data
        clip.savewrk(ceres)
        # record the results
        target.record(clip, ceres)

Yay, it was that easy and this script can be called every morning. Although modifications should be made 
to allow for you to set the target and control without hardcoding like we've done here as its probable that 
you would probably like to work with more than one target over time.

Dorado was built to allow for this sort of task to be done painlessly.
