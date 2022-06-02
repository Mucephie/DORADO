Intro to the reader class system
================================= 

The first iteration of Dorado was known as Draco (Data reduction at the Allan Carswell Observatory) way back in 2018/2019.
It was specific to the way the Allan I. Carswell Observatory (AICO) handled undergraduate variable star research which was
previously handled using IRAF. Dorado grew as AICO grew and sought to accommodate the plethora of changes including data 
operations involving the newly installed Planewave CDK-1000 1m telescope. 

Dorado at current epoch uses a 'modular' data reader to create, locate, navigate, read, and manage astronomical data within 
the dorado directory. Users can create their own data reader specific to their data pipeline needs or use one of the prebuilt
data readers provided in the dorado package such as aico_reader(). <-- **link to the aico_reader docs**

Further information on readers is coming as the reader system matures.

Previous: :doc:`Dorado Core</dorado/core>` || Next: :doc:`Ceres</dorado/ceres>`

