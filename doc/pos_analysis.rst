Toolchain for positioner analysis
=================================
The general sequence of events is:

#) Positioners are moved on the instrument.

   #) This could be at the Mayall or in the lab.

#) Calculate measured fiber positions from FVC images:

   #) Platemaker, at time of measurement.
   #) Or desimeter, at any time later.

#) Data is selected from the online database (posmovedb), and stored to local CSV files.

#) Analyze / plot.

Moving positioners + measuring
------------------------------
This happens outside of desimeter. The robots are always moved using ``PetalApp.py``, which is an ICS wrapper around ``petal.py``.

``<desi_svn>/code/online/Positioner/PetalApp/<some_tag_or_branch>/python/PetalApp/PetalApp.py``

``<desi_svn>/code/focalplane/plate_control/<some_tag_or_branch>/petal/petal.py``

Depending on the operational circumstance, there will be different high-level scripts which send the requests for which targets to go to, and when to trigger the FVC camera. These are in general the night time OCS script, or else for "xytests" and "arc calibrations", etc, the various PECS scripts.

``<desi_svn>/code/online/OCS/<some_tag_or_branch>/python/OCS/``

``<desi_svn>/code/focalplane/plate_control/<some_tag_or_branch>/pecs/``

For every move performed, numerous values get stored in the online database (posmovedb). These include:

      #) ``POS_T`` and ``POS_P``, the internally-tracked (theta, phi) angles of each positioner
      
      #) ``OBS_X`` and ``OBS_Y``, the global (x, y) measured position calculated by platemaker
      
      #) ``PTL_X``, ``PTL_Y``, and ``PTL_Z``, the platemaker-measured position, transformed to petal coordinates, using each petal's (x, y, z, alpha, beta, gamma) rigid body position as known to the instrument at the time of measurement
      
      #) ``EXPOSURE_ID`` and ``EXPOSURE_ITER``, index values to identify associated FVC FITS file with that row in posmovedb
      
      #) ``MOVE_CMD``, ``MOVE_VAL1``, ``MOVE_VAL2``, ``LOG_NOTE``, fields that can be parsed to determine what higher level operation was being performed

*(As of 2020-05-12, items (2), (3), (4) are in the process of being brought online.)*

Anytime a calibration value is changed (such as ``LENGTH_R1``, ``LENGTH_R2``, etc), or indeed any other change of positioner state data, that fact will be recorded in the online database.

Re-analysis of fiber positions with desimeter
---------------------------------------------
FVC fits files can be re-analyzed at any time after the fact using desimeter. Also see: `<posmov.rst>`_

Retreiving data from the online database
----------------------------------------
The posmovedb for the instrument at the Mayall is hosted at KPNO, and regularly mirrored to NERSC. There is a similar posmovedb hosted at LBNL on the beyonce server, just for test petals there.

As of 2020-05-12, desimeter provides the following tools:

* ``get_posmov_fvc_data`` ... Matches timestamps between posmovedb and fvc FITS files, saves csv combining posmovedb rows with desimeter’s analysis of the measured positions.

* ``get_posmov_calib`` ... Just grabs values like ``LENGTH_*``, ``OFFSET_*`` etc from posmovedb, saves csv.

Analysis of positioning performance
-----------------------------------

As of 2020-05-12, desimeter provides the following tools:

* ``fit_posparams`` ... Fits positioner calibration parameters, by comparing measured (x,y) to internally-tracked (theta,phi).

* ``plot_posparams`` ... Plots results from fit_posparams per positioner. Also plots cumulative positioner errors over time, as calculated when performing those best-fits.

* ``plot_positioner`` ...

* ``plot_posmove`` ...
	

