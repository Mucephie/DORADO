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



Using Clippy
============



Previous: :doc:`Getting Started</dorado/GettingStarted>` || Next: :doc:`Ceres</dorado/ceres>`



