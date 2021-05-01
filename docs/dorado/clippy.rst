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
While Dorado is pretty flexile in terms of functionality and use, Dorado does apply some strict procedures
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
    ├── targets          
    ├── logs
    ├── cache         
    │   └── astrometryNet

Lets look at the ``/data`` folder, raw astronomical data is stored in subfolders within the ``/raw`` folder.
Subfolders in the raw data directory should take a naming scheme based on the date of observation
``YYYY-MM-DD+DD`` where the ``DD+DD`` portion denotes the date before and after local midnight during the 
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



Lets imagine that we just did an observing run 



Using Clippy
============



Previous: :doc:`Getting Started</dorado/GettingStarted>` || Next: :doc:`Ceres</dorado/ceres>`



