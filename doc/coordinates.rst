=======================================
Coordinates
=======================================

This document presents the different coordinate systems
 and the coordinate transformation routines in desimeter.

Coordinate systems
++++++++++++++++++

GFA (Guide Focus Array)
~~~~~~~~~~~~~~~~~~~~~~~

Pixel coordinates in GFA guide CCD images. Usually called `x_gfa , y_gfa` in desimeter.
Measurements are done by Aaron Meisner's code.

PTL (Petal)
~~~~~~~~~~~

Petal coordinates in mm, either as cartesian coordinates `x_ptl,y_ptl,z_ptl` or curved coordinates `s_ptl,q_ptl`. This coordinate system is attached to the metrology measurements for positioner mounting holes, FIF and GIF. See `data.README`_

FP (Focal Plane)
~~~~~~~~~~~~~~~~

Coordinate system as close as possible to the 'CS5' system described in ...
PTL and FP are related by a translation and a rotation. `x_fp,y_fp`  are also in mm. It is not clear yet whether the petals are moving one with respect to the others.

TP (Tangent Plane)
~~~~~~~~~~~~~~~~~~


Polar coordinates (theta,phi) projected on a plane with an arbitrary but fixed projection::

  x_tan = sin(theta) cos(phi)
  y_tan = sin(theta) sin(phi)

For a telescope pointing at `HA_tel , Dec_tel`,
`phi=0` for a star with `HA>HA_tel` and `Dec=Dec_tel`, and `phi=90 deg` for a star with `HA=HA_tel` and `Dec>Dec_tel`.

The angles (theta,phi) are the direction of photons hitting the telescope (so including refraction but excluding optical distortions) in a frame attached to the focal plane (the frame is rotating with the hexapod).

ICRS (International Celestial Reference System)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`RA,Dec` sky coordinates of targets (same system as GAIA stars used for the astrometry).



.. _`data.README`: https://github.com/desihub/desimeter/blob/master/py/desimeter/data/README.rst
.. _`DESI-5421`: https://desi.lbl.gov/DocDB/cgi-bin/private/ShowDocument?docid=5421
