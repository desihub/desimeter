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

Depending on the operational circumstance, there will be different high-level scripts which send the requests for which targets to go to, and when to trigger the FVC camera. These are in general:

Night time observations --> OCS script
    ``<desi_svn>/code/online/OCS/<some_tag_or_branch>/python/OCS/``

Scripted move sequences --> run_sequence.py
    ``<desi_svn>/code/focalplane/plate_control/<some_tag_or_branch>/pecs/``

For every move performed, numerous values get stored in the online database (posmovedb). These include:

      #) ``POS_T`` and ``POS_P``, the internally-tracked (theta, phi) angles of each positioner
      
      #) ``OBS_X`` and ``OBS_Y``, the global (x, y) measured position calculated by platemaker
      
      #) ``PTL_X``, ``PTL_Y``, and ``PTL_Z``, the platemaker-measured position, transformed to petal coordinates, using each petal's (x, y, z, alpha, beta, gamma) rigid body position as known to the instrument at the time of measurement
      
      #) ``EXPOSURE_ID`` and ``EXPOSURE_ITER``, index values to identify associated FVC FITS file with that row in posmovedb
      
      #) ``MOVE_CMD``, ``MOVE_VAL1``, ``MOVE_VAL2``, ``LOG_NOTE``, fields that can be parsed to determine what higher level operation was being performed

Anytime a calibration value is changed (such as ``LENGTH_R1``, ``LENGTH_R2``, etc), or indeed any other change of positioner state data, that fact will be recorded in the online database.

Re-analysis of fiber positions with desimeter
---------------------------------------------
FVC fits files can be re-analyzed at any time after the fact using desimeter. Also see: `<posmov.rst>`_

Retreiving data from the online database
----------------------------------------
The posmovedb for the instrument at the Mayall is hosted at KPNO, and regularly mirrored to NERSC. There is a similar posmovedb hosted at LBNL on the beyonce server, just for test petals there.

As of 2020-06-16, desimeter provides the following tools:

* ``get_posmov_fvc_data`` ... Matches timestamps between posmovedb and fvc FITS files, saves csv combining posmovedb rows with desimeter’s analysis of the measured positions.

* ``get_posmoves`` ... Retrieves measured move values values ``PTL_X_*``, ``PTL_Y_*`` etc posmovedb, saves csv.

* ``get_posparams`` ... Retrieves calibration values like ``LENGTH_*``, ``OFFSET_*`` etc from posmovedb, saves csv.

Analysis of positioning performance
-----------------------------------

As of 2020-06-16, desimeter provides the following tools:

* ``fit_posparams`` ... Fits positioner calibration parameters, by comparing measured (x,y) to internally-tracked (theta,phi).

* ``plot_posparams`` ... Plots results from fit_posparams per positioner. Also plots cumulative positioner errors over time, as calculated when performing those best-fits.

Preparation of calibration values for use on the instrument
-----------------------------------------------------------

As of 2020-06-16, desimeter provides the following tool:

* ``prepare_posparams_for_instr`` ... Takes a csv file produced by fit_posparams and interactively guides user through validating and selecting safe values for use on instrument. Generates an output csv file which can be ingested by ``pecs/set_calibrations.py``

The procedure for measuring and updating calibrations is given in DESI-5732.

As on overview of the role desimter plays in this process, those basic steps are:

1. ``get_posmoves`` ... get tracked (t,p) and measured (x,y) from online DB

2. ``fit_posparams`` ... best-fit calib params which map (t,p) to (x,y)

3. ``merge_posparams`` ... gather fit result files into one table

4. ``prepare_posparams_for_instr`` ... validate parameters and generate modified table

Finally after desimeter has prepared the new calibration parameters, we use ``pecs/set_calibrations.py`` (managed in DESI's svn repo, *not* desihub), which pushes data to the online DB.
	

