Toolchain for positioner analysis
=================================
The general sequence of events is:

#) Positioners are moved on the instrument.

   #) This could be at the Mayall or in the lab.
   #) Values get stored in the online database (posmovedb) for:
   
      #) ``POS_T`` and ``POS_P``, the internally-tracked (theta, phi) angles of each positioner
      #) ``OBS_X`` and ``OBS_Y``, the global (x, y) measured position calculated by platemaker [1]_
      #) ``PTL_X``, ``PTL_Y``, and ``PTL_Z``, the platemaker-measured position, transformed to petal coordinates,
      using each petal's (x, y, z, alpha, beta, gamma) rigid body position as known to the instrument at the
      time of measurement.
      
.. [1] As of 2020-05-12, this is just being brought online.
