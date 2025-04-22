====================
desimeter change Log
====================

0.7.2 (unreleased)
------------------

* Switch integration testing from CirclCI to GitHub actions [PR `#197`_].

.. _`#197`: https://github.com/desihub/desimeter/pull/197

0.7.1 (2024-08-15)
------------------

* Add support for test Petal0 in lab (direct commits, not PR).

0.7.0 (2023-01-12)
------------------

* Add turbulence correction mode using stationary positioners [PR `#186, #187`_].
* Refactor desi_fvc_proc and add option to run directly w/o subprocess
  [PR `#189, #192`_].
* Clip outlying stationary positioners from turbulence fits at 10 microns
  [PR `#190, #193`_].
* Change print to log [PR `#195`_]

.. _`#186`: https://github.com/desihub/desimeter/pull/186
.. _`#187`: https://github.com/desihub/desimeter/pull/187
.. _`#189`: https://github.com/desihub/desimeter/pull/189
.. _`#190`: https://github.com/desihub/desimeter/pull/190
.. _`#192`: https://github.com/desihub/desimeter/pull/192
.. _`#193`: https://github.com/desihub/desimeter/pull/193
.. _`#195`: https://github.com/desihub/desimeter/pull/195

0.6.7 (2021-05-28)
------------------

* Use vectorized cross product for speed [PR `#185`_].

.. _`#185`: https://github.com/desihub/desimeter/pull/185

0.6.6 (2021-05-13)
------------------

* Add empirical desimeter to platemaker transforms and use for
  fiberassign xy2radec and radec2xy [PR `#184`_]

.. _`#184`: https://github.com/desihub/desimeter/pull/184

0.6.5 (2021-05-11)
------------------

* add ``get_posmoves --latest`` option [PRs `#179, #181, #182`_]
* add ``analyze_pos_performance`` options to filter based upon positioner
  movement scale factors [PR `#180`_]
* add ``desimeter.fiberassign`` with x,y <-> ra,dec transform predictions
  for fiberassign (i.e. not using as-observed GFA data) [PR `#183`_]

.. _`#179`: https://github.com/desihub/desimeter/pull/179
.. _`#180`: https://github.com/desihub/desimeter/pull/180
.. _`#181`: https://github.com/desihub/desimeter/pull/181
.. _`#182`: https://github.com/desihub/desimeter/pull/182
.. _`#183`: https://github.com/desihub/desimeter/pull/183

0.6.4 (2021-03-31)
------------------

* posperf fix for normal tiles [PR `#177`_]
* Analyze pos options for large datasets [PR `#178`_]

.. _`#177`: https://github.com/desihub/desimeter/pull/177
.. _`#178`: https://github.com/desihub/desimeter/pull/178

0.6.3 (2021-02-15)
------------------

* Update posflags [PR `#173`_]
* Petal 0 in LBL 6040 [PR `#174`_]
* Put row deletion in conditional [PR `#175`_]
* Bright focal plane image analysis [PR `#176`_]

.. _`#173`: https://github.com/desihub/desimeter/pull/173
.. _`#174`: https://github.com/desihub/desimeter/pull/174
.. _`#175`: https://github.com/desihub/desimeter/pull/175
.. _`#176`: https://github.com/desihub/desimeter/pull/176

0.6.2 (2020-12-18)
------------------

* Update to field rotation model [PR `#172`_]

.. _`#172`: https://github.com/desihub/desimeter/pull/172

0.6.1 (2020-12-11)
------------------

* Fixes in prepare_posparams used for FP calibration of 20201210 [PR `#167`_]

.. _`#167`: https://github.com/desihub/desimeter/pull/167

0.6.0 (2020-12-10)
------------------

* Update to prepare_posparams [PR `#166`_]
* Several changes for KPNO focal plane calibration [PR `#165`_]
* Increase max bin size to accommodate corner case [PR `#164`_]
* Make new fvc2fp parameters optional. [PR `#163`_]
* Disable bad performers [PR `#160`_]
* Spots coordinates statistical uncertainties [PR `#159`_]
* Unmask POSTSCRIPT and CALIB_NOTE fields when necessary [PR `#158`_]
* Get enabled history [PR `#157`_]
* Better handle ADC settings in unusual configurations. [PR `#156`_]
* Empirical field rotation rate [PR `#155`_]
* Remove parsing 'freeze' in log_note string from inclusion as "FROZEN" [PR `#154`_]
* Fix some ADC angle computations [PR `#153`_]
* Update to spotmatch usage [PR `#152`_]
* Let desimeter routines make directories for output if needed. [PR `#151`_]
* Update X_FP, Y_FP for GFA pixel corners to match target_ra/dec. [PR `#150`_]
* Accepted targets [PR `#149`_]
* Update fp-metrology file based on fiducial systematics script. [PR `#148`_]
* Add DEVICE_ID,BUS_ID... to metrology table using online positioner index table. [PR `#146`_]
* Set dtype to object for cols that can carry None or int [PR `#144`_]
* Clear tmp before analyze bright [PR `#143`_]
* plot_posparams -- skip files where data is missing [PR `#142`_]
* Add option --force-recent-rehome to fit_posparams [PR `#141`_]
* Let fetched posmovedb rows have same timestamp as calib [PR `#140`_]
* Pos analysis updates [PR `#139`_]
* Posmoves tp updates [PR `#138`_]
* Identify requests [PR `#137`_]
* Dither analysis [PR `#136`_]
* Proper motion [PR `#135`_]
* Take into account proper motions [PR `#134`_]
* Bug fix when reading hexapod rotation in fits header [PR `#127`_]
* Turbulence opt [PR `#126`_]
* Turbulence modeling code. [PR `#125`_]
* Use desimeter logger [PR `#124`_]
* Analyze pos perf updates [PR `#123`_]
* Doc xy analysis [PR `#122`_]
* Use device centers [PR `#121`_]
* Spotmatch [PR `#120`_]

.. _`#120`: https://github.com/desihub/desimeter/pull/120
.. _`#121`: https://github.com/desihub/desimeter/pull/121
.. _`#122`: https://github.com/desihub/desimeter/pull/122
.. _`#123`: https://github.com/desihub/desimeter/pull/123
.. _`#124`: https://github.com/desihub/desimeter/pull/124
.. _`#125`: https://github.com/desihub/desimeter/pull/125
.. _`#126`: https://github.com/desihub/desimeter/pull/126
.. _`#127`: https://github.com/desihub/desimeter/pull/127
.. _`#134`: https://github.com/desihub/desimeter/pull/134
.. _`#135`: https://github.com/desihub/desimeter/pull/135
.. _`#136`: https://github.com/desihub/desimeter/pull/136
.. _`#137`: https://github.com/desihub/desimeter/pull/137
.. _`#138`: https://github.com/desihub/desimeter/pull/138
.. _`#139`: https://github.com/desihub/desimeter/pull/139
.. _`#140`: https://github.com/desihub/desimeter/pull/140
.. _`#141`: https://github.com/desihub/desimeter/pull/141
.. _`#142`: https://github.com/desihub/desimeter/pull/142
.. _`#143`: https://github.com/desihub/desimeter/pull/143
.. _`#144`: https://github.com/desihub/desimeter/pull/144
.. _`#146`: https://github.com/desihub/desimeter/pull/146
.. _`#148`: https://github.com/desihub/desimeter/pull/148
.. _`#149`: https://github.com/desihub/desimeter/pull/149
.. _`#150`: https://github.com/desihub/desimeter/pull/150
.. _`#151`: https://github.com/desihub/desimeter/pull/151
.. _`#152`: https://github.com/desihub/desimeter/pull/152
.. _`#153`: https://github.com/desihub/desimeter/pull/153
.. _`#154`: https://github.com/desihub/desimeter/pull/154
.. _`#155`: https://github.com/desihub/desimeter/pull/155
.. _`#156`: https://github.com/desihub/desimeter/pull/156
.. _`#157`: https://github.com/desihub/desimeter/pull/157
.. _`#158`: https://github.com/desihub/desimeter/pull/158
.. _`#159`: https://github.com/desihub/desimeter/pull/159
.. _`#160`: https://github.com/desihub/desimeter/pull/160
.. _`#163`: https://github.com/desihub/desimeter/pull/163
.. _`#164`: https://github.com/desihub/desimeter/pull/164
.. _`#165`: https://github.com/desihub/desimeter/pull/165
.. _`#166`: https://github.com/desihub/desimeter/pull/166

0.5.1 (2020-09-01)
------------------

* Add xytest analysis code [PR `#119`_]
* Fix xy2tp phi wrapping error [PR `#118`_]
* Fix bug when setting IDs to positioners for petal 0 [PR `#117`_]
* Add Front illuminated focal plane "bright image" analysis code [PR `#116`_]
* Add robustness to NaN in coordinates [PR `#115`_]
* Update GFA metrology using onsky data [PR `#113`_]
* Improve xy2tp [PR `#112`_]
* Add ptl2int convenience wrapper function [PR `#111`_]

.. _`#111`: https://github.com/desihub/desimeter/pull/111
.. _`#112`: https://github.com/desihub/desimeter/pull/112
.. _`#113`: https://github.com/desihub/desimeter/pull/113
.. _`#115`: https://github.com/desihub/desimeter/pull/115
.. _`#116`: https://github.com/desihub/desimeter/pull/116
.. _`#117`: https://github.com/desihub/desimeter/pull/117
.. _`#118`: https://github.com/desihub/desimeter/pull/118
.. _`#119`: https://github.com/desihub/desimeter/pull/119

0.5.0 (2020-08-05)
------------------

* Added teststand petal 0 metrology [PR `#109`_]
* Patch pseudo-metrology for Guide3 and Guide7 GFAs [PR `#108`_]
* Updated raytracing model [PR `#107`_]
* Added tools for analysing front-illuminated images [PR `#104`_, `#105`_]
* Added fiber matching disambiguation code [PR `#103`_]
* Added tools for matching moving spots to move requests [PR `#101`_]
* Convert ptlxy to flat, improve log message detail [PR `#99`_]
* Combined posmoves from DB with FVC image data [PR `#98`_]
* Script to prepare calib posparmas for DB upload [PR `#97`_]

.. _`#97`: https://github.com/desihub/desimeter/pull/97
.. _`#98`: https://github.com/desihub/desimeter/pull/98
.. _`#99`: https://github.com/desihub/desimeter/pull/99
.. _`#101`: https://github.com/desihub/desimeter/pull/101
.. _`#103`: https://github.com/desihub/desimeter/pull/103
.. _`#104`: https://github.com/desihub/desimeter/pull/104
.. _`#105`: https://github.com/desihub/desimeter/pull/105
.. _`#107`: https://github.com/desihub/desimeter/pull/107
.. _`#108`: https://github.com/desihub/desimeter/pull/108
.. _`#109`: https://github.com/desihub/desimeter/pull/109

0.4.0 (2020-06-15)
------------------

* Positioner calibration tools [PR `#73`_, `#77`_, `#78`_, `#79`_, `#81`_,
  `#82`_, `#84`_, `#85`_, `#86`_, `#87`_, `#88`_, `#89`_, `#90`_, `#91`_,
  `#92`_, `#93`_, `#94`_, `#95`_, `#96`_]
* Codacy [PR `#75`_, `#76`_]
* Support astropy 3.0.4 and matplotlib 2.1.2 [PR `#80`_]

.. _`#73`: https://github.com/desihub/desimeter/pull/73
.. _`#75`: https://github.com/desihub/desimeter/pull/75
.. _`#76`: https://github.com/desihub/desimeter/pull/76
.. _`#77`: https://github.com/desihub/desimeter/pull/77
.. _`#78`: https://github.com/desihub/desimeter/pull/78
.. _`#79`: https://github.com/desihub/desimeter/pull/79
.. _`#80`: https://github.com/desihub/desimeter/pull/80
.. _`#81`: https://github.com/desihub/desimeter/pull/81
.. _`#82`: https://github.com/desihub/desimeter/pull/82
.. _`#84`: https://github.com/desihub/desimeter/pull/84
.. _`#85`: https://github.com/desihub/desimeter/pull/85
.. _`#86`: https://github.com/desihub/desimeter/pull/86
.. _`#87`: https://github.com/desihub/desimeter/pull/87
.. _`#88`: https://github.com/desihub/desimeter/pull/88
.. _`#89`: https://github.com/desihub/desimeter/pull/89
.. _`#90`: https://github.com/desihub/desimeter/pull/90
.. _`#91`: https://github.com/desihub/desimeter/pull/91
.. _`#92`: https://github.com/desihub/desimeter/pull/92
.. _`#93`: https://github.com/desihub/desimeter/pull/93
.. _`#94`: https://github.com/desihub/desimeter/pull/94
.. _`#95`: https://github.com/desihub/desimeter/pull/95
.. _`#96`: https://github.com/desihub/desimeter/pull/96

0.3.0 (2020-05-05)
------------------
* Connection to posmov db at LBL, documentation [PR `#70`_]
* LBL petal 1 [PR `#69`_]
* Make posparam fits more atomic [PR `#68`_]
* Added CircleCI and Coveralls [PR `#57`_] [PR `#58`_] [PR `#62`_]
* pos2ptl [PR `#56`_] [PR `#59`_]
* simplification of call to ptl2fp [PR `#55`_]
* Added `fit_posparams` and `desimeter.transform.pos2ptl` for fitting
  positioner parameters [PR `#54`_].

.. _`#70`: https://github.com/desihub/desimeter/pull/70
.. _`#69`: https://github.com/desihub/desimeter/pull/69
.. _`#68`: https://github.com/desihub/desimeter/pull/68
.. _`#63`: https://github.com/desihub/desimeter/pull/63
.. _`#62`: https://github.com/desihub/desimeter/pull/62
.. _`#59`: https://github.com/desihub/desimeter/pull/59
.. _`#58`: https://github.com/desihub/desimeter/pull/58
.. _`#57`: https://github.com/desihub/desimeter/pull/57
.. _`#56`: https://github.com/desihub/desimeter/pull/56
.. _`#55`: https://github.com/desihub/desimeter/pull/55
.. _`#54`: https://github.com/desihub/desimeter/pull/54

0.2.1 (2020-04-15)
------------------

* Simplified call in ptl2fp [PR `#55`_]

.. _`#55`: https://github.com/desihub/desimeter/pull/55

0.2.0 (2020-04-06)
------------------

* Improves fidicial pinhole matching robustness [PR `#15`_]
* Adds ability to match fibers to expected positions [PR `#17`_]
* Added tangent plane to/from focal plane transforms [PR `#21`_]
* Added ra,dec to/from tangent plane transforms [PR `#22`_, `#29`_, `#30`_]
* Fix NotImplementedError typos [PR `#28`_]
* Added GFA to/from focal plane transforms [PR `#31`_, `#46`_]
* Fit guide star coordinates [PR `#34`_]
* Tools to calculate per-fiber RA,dec from field model [PR `#35`_]
* Option to use last extension of FVC file [PR `#36`_]
* Fail more gracefully if very few fiducials are detected [PR `#38`_]
* Add model for new singlet FVC lens [PR `#40`_]
* Add ray trace model of corrector distortions including ADC [PR `#41`_]
* Update to metrology v6 [PR `#44`_]
* Account for z-offset of GFA sensors [PR `#45`_]
* Add field rotation prediction [PR `#48`_]
* Update to metrology v7 [PR `#50`_]
* Added tools to fit positioner calibration circles [PR `#53`_]

.. _`#15`: https://github.com/desihub/desimeter/pull/15
.. _`#17`: https://github.com/desihub/desimeter/pull/17
.. _`#21`: https://github.com/desihub/desimeter/pull/21
.. _`#22`: https://github.com/desihub/desimeter/pull/22
.. _`#28`: https://github.com/desihub/desimeter/pull/28
.. _`#29`: https://github.com/desihub/desimeter/pull/29
.. _`#30`: https://github.com/desihub/desimeter/pull/30
.. _`#31`: https://github.com/desihub/desimeter/pull/31
.. _`#34`: https://github.com/desihub/desimeter/pull/34
.. _`#35`: https://github.com/desihub/desimeter/pull/35
.. _`#36`: https://github.com/desihub/desimeter/pull/36
.. _`#38`: https://github.com/desihub/desimeter/pull/38
.. _`#40`: https://github.com/desihub/desimeter/pull/40
.. _`#41`: https://github.com/desihub/desimeter/pull/41
.. _`#44`: https://github.com/desihub/desimeter/pull/44
.. _`#45`: https://github.com/desihub/desimeter/pull/45
.. _`#46`: https://github.com/desihub/desimeter/pull/46
.. _`#48`: https://github.com/desihub/desimeter/pull/48
.. _`#50`: https://github.com/desihub/desimeter/pull/50
.. _`#53`: https://github.com/desihub/desimeter/pull/53

0.1.0 (2019-12-29)
------------------

* Initial release
