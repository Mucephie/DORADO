Ceres and stack classes: Observational data handling
====================================================

To interface with an observational series, Dorado offers the Ceres class. The Ceres class contains a list 
of Stack objects and the coresponding metadata including :class:`~dorado.Target` objects associated with
the campaign.


Astronomical data is collected and stored in many different ways and formats and for dorado to truly be an open and community 
friendly pipeline tool, dorado needs a modular and personalizable way to integrate data into the dorado system. One of the 
universal data objects in dorado is the :class:`~dorado.Ceres` class which handles data associated with a nights observation 
campaign, including individual stacks of data held in :class:`~dorado.Stack` data objects. The reader class facilitates the 
reading of the data into classes like :class:`~dorado.Ceres` and :class:`~dorado.Stack` where it can be operated on in a consistant
manor. 


Previous: :doc:`Reader</dorado/reader>` || Next: :doc:`Targets</dorado/target>`