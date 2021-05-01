Getting to know Clippy
======================

You may remember that condescending paperclip from the 90s that helped a computer iliterate
world cope with CRT displays, cleaning mouse balls, and the blue screen of death; well now 
the ghost of clippy is back to help you manage and organized an effective local astronomical
file system. 

.. image:: /dorado/_static/clippy.png

In Dorado, Clippy is a class object and does not currently manifest omnipresent pop-ups offering help
like the Clippy of old. Instead, you as the user will import and then initialize Clippy for use:

.. code:: python

        # import Clippy from Dorado
        from dorado.clippy import Clippy

        # initialize an instance
        clip = Clippy()

Now that Cippy has essentially been summoned into your script and thus into your python interpreter 
and machines memory at runtime, you may be asking yourself "Why? What does it even do? And why do I
even need this janky Clippy immitation anyways?". Well, despite the name and affinity for assistance, 
Dorado's Clippy and the Clippy of old share nothing more in common. Dorado's Clippy is akin to an automatic
file and directory assistant, Clippy can create, locate, navigate, read, and manage astronomical data
directories.

The Clippy class you now have in your script will serve as the backbone of our example Dorado pipeline.
While Dorado is pretty flexible in terms of functionality and use, Dorado does apply some strict procedures
to keep things organized and universal; this is most apparent within Clippy and manifest itself as the
``.dorado`` directory which will be created either upon install of Dorado or upon the first time Clippy is 
initialized. This directory when created has the following format:

::

    .dorado
    ├── data          
    │   ├── wrk
    │   └── bias
    │   └── flats
    │   └── raw
    │   └── graphical
    │   └── projects
    │   └── targets        
    ├── logs
    ├── cache         
    │   └── astrometryNet

Lets look at the ``/data`` folder, raw astronomical data is stored in subfolders within the ``/raw`` folder.
Subfolders in the raw data directory should take a naming scheme based on the date of observation
``/YYYY-MM-DD+DD`` where the ``DD+DD`` portion denotes the date before and after local midnight during the 
observing night. The inside of the nights folder can take the form of one of two basic file structures; single
directory, where all files are stored within the folder with no subdirectories, and multi level directory,
where there is a folder for lights, bias, and flatfield frames. If multiple filters were used to observe, multilevel 
directory format allows for each filter to have its own subdirectory in the lights and flatfield directories.

The ``/bias`` and ``/flats`` directory store the base bias and flatfield images produced by Dorado for each relevent 
date. In the event that an observing session lacked these calibration frames, Clippy will consult these two directories
to retrieve the closest in date calibration frames to use for calibration. 

    .. note:: In the near future this functionality will be expanded upon to allow constraints on the max date separation 
        allowed between the retrieved calibration frames and the date of observation; as well as the ability to compare these 
        calibration frames over time to see how your imaging setup changes and thus allow you to refine your max date separation
        constraints.

The ```/wrk`` directory contains the working folders for your processed data, this is the default location for any
calibrated, aligned, or otherwise worked with images. Data saved will follow the ``/YYYY-MM-DD+DD`` subdirectory
naming scheme.

At the current epoch, the ``/graphical`` and ``/projects`` subdirectories are unused but in the near future ``/graphical`` will  
store any finalized stacked or colour merged images, graphs, or timeseries plots, and ``/projects`` will be a place to store data 
pertaining to ongoing projects that span more than one observing session and target. In the future the community will be consulted 
as to how best organize and interact with the project directory.

The ``/targets`` directory is a directory used to store target data used by Dorado for future use in th context of the 
:doc:`Zellars</dorado/zellars>` target class object. 

The ``/logs`` folder is th future home of log files produced by Dorado during operation, this may either be during observations
or during subsequent data operations. These logs can then be used to debug pipelines and Dorado or used for future reference
as to what actions were performed on a dataset. 

The ``/cache`` folder is a place for Clippy to temporarily store files to the hard disk while performing certain tasks such as uploading
frames to 'Astrometry.net' and downloading the results suring runtime. The cache may also be used to store temporary log files
or any other file objects that are temporary in nature and need a place to exist where it is easily acessible while also easy to remove
via cache cleaning commands.

Using Clippy
============


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


Previous: :doc:`Getting Started</dorado/GettingStarted>` || Next: :doc:`Ceres</dorado/ceres>`



