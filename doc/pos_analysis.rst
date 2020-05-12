Toolchain for positioner analysis
=================================
The general sequence of events is:

#) Positioners are moved on the instrument.

   #) This could be at the Mayall or in the lab.

#) Data is selected from the online database (posmovedb), and stored to local CSV files.

#) Derive measured fiber positions:

   #) Desimeter calculates measured fiber positions from FVC images.
   #) Or skip this step and use platemaker's values, grabbed from posmovedb. [1]_

#) 

   
   
   #) Values get stored in the online database (posmovedb) for:
   
      #) ``POS_T`` and ``POS_P``, the internally-tracked (theta, phi) angles of each positioner
      
      #) ``OBS_X`` and ``OBS_Y``, the global (x, y) measured position calculated by platemaker [1]_
      
      #) ``PTL_X``, ``PTL_Y``, and ``PTL_Z``, the platemaker-measured position, transformed to petal coordinates, using each petal's (x, y, z, alpha, beta, gamma) rigid body position as known to the instrument at the time of measurement [1]_
      
      #) ``EXPOSURE_ID`` and ``EXPOSURE_ITER``, index values to identify associated FVC FITS file with that row in posmovedb [1]_
      
      
      
.. [1] As of 2020-05-12, this is just now being brought online.
