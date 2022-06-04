Ceres and stack classes: Observational data handling
====================================================

Ceres
------

To interface with an observational series, Dorado offers the Ceres class. The Ceres class contains a list 
of Stack objects and the coresponding metadata including :class:`~dorado.Target` objects associated with
the campaign.As stated before, Astronomical data is collected and stored in many different ways and formats;
Dorado aims to be an open and community friendly pipeline tool, to do this, dorado needs a modular and 
personalizable way to integrate data into the dorado system. 

The :class:`~dorado.Ceres` class handles data associated with a nights observation campaign, including 
individual stacks of data held in :class:`~dorado.Stack` data objects, and data on the observations themselves. 
Ceres objects are by convention stored in the core dorado instance and referenced by key strings for ease 
of use.


.. code:: python

        ## import
        import dorado
        ## initialize
        # create a target object 'target of interest' toi
        toi = 'target_name'
        dorado.mktrgt(toi)
        # create a series object
        dorado.mkceres('2021-01-01+02', name = 'ceres_key', target = toi)

        ## call
        # calibrate the series red filter data
        dorado.dorphot.calibrate('ceres_key','filter_str')




Stack
------

The :class:`~dorado.Stack`

Previous: :doc:`Reader</dorado/reader>` || Next: :doc:`Targets</dorado/target>`