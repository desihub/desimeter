=======================================
Coordinates
=======================================

This document presents the different coordinate systems and the coordinate transformation routines in desimeter.

Coordinate systems
++++++++++++++++++

GFA (Guide Focus Array)
~~~~~~~~~~~~~~~~~~~~~~~

Pixel coordinates in GFA guide CCD images. Usually called `x_gfa , y_gfa` in desimeter.
Measurements are done by Aaron Meisner's code.

FVC (Fiber View Camera)
~~~~~~~~~~~~~~~~~~~~~~~

Pixel coordinates in the Fiber View Camera CCD images. Usually called `x_fvc , y_fvc` in desimeter.
In this code, the coordinates of the first pixel is `(0,0)`. This is different from other codes. 

PTL (Petal)
~~~~~~~~~~~

Petal coordinates in mm, either as cartesian coordinates `x_ptl,y_ptl,z_ptl` or curved coordinates `s_ptl,q_ptl`. This coordinate system is attached to the metrology measurements for positioner mounting holes, FIF and GIF. See `data.README`_ .

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


Coordinate transformations
++++++++++++++++++++++++++

GFA -> FP
~~~~~~~~~

Based on the metrology compilation in `DESI-5421`_ + rigid transformation from PTL to FP.
The fit of the metrology is performed on the fly for the first call of the function.

Example::

  from desimeter.transform.gfa2fp import gfa2fp,fp2gfa
  x_fp,y_fp   = gfa2fp(petal_loc,x_gfa,y_gfa)
  x_gfa,y_gfa = fp2gfa(petal_loc,x_fp,y_fp)

FVC -> FP
~~~~~~~~~

Based on a fit of the measured pixel coordinates of fiducials and their metrology.
In most and probably all cases, the transformation is used in the same procedure
 as the fit (see the script `desi_fvc_proc`). Several implementations of this
 transform exist in `transform/fvc2fp`.

 Example::

  from desimeter.transform.fvc2fp.zb import FVCFP_ZhaoBurge
  tx = FVCFP_ZhaoBurge()
  tx.fit(spots, update_spots=True)
  tx.write_jsonfile(filename)
  x_fp,y_fp = tx.fvc2fp(x_fvc,y_fvc)
  
FP -> TAN
~~~~~~~~~

From the focal surface instruments to the tangent plane. This contains all of the optics
distortion. This transformation is minimal for now. It is based on the echo22 design model
(the default `echo22.py` transformations are "linked" in `transform/tan2fp/__init__.py`).
We have to do much better than this.

 Example::

   from desimeter.transform.tan2fp import fp2tan,tan2fp
   x_tan,y_tan = fp2tan(x_fp,y_fp)

TAN -> ICRS
~~~~~~~~~~~

This transformation includes precession, aberration, refraction, polar mis-alignment, hexapod rotation
angle.

 Example::

   from desimeter.transform.radec2tan import radec2tan,tan2radec
   ra,dec = radec2tan(x_tan,y_tan,telescope_ra,telescop_dec,mjd,lst,precession=True,aberration=True)

FieldModel
~~~~~~~~~~

This is not a transformation from one system to another but rather the combination of the transformations, including a correction based on the location of GFA guide stars. This "field model" is fit to the data
 using the script `desi_fit_guide_star_coordinates`.

Example::

  import json
  from desimeter.fieldmodel import FieldModel
  file=open("model.json")
  fm = FieldModel.fromjson(file.read())
  ...
  x_fp,y_fp = fm.all_gfa2fp(x_gfa,y_gfa,petal)
  ra,dec=fm.fp2radec(x_fp,y_fp)

  
.. _`data.README`: https://github.com/desihub/desimeter/blob/master/py/desimeter/data/README.rst
.. _`DESI-5421`: https://desi.lbl.gov/DocDB/cgi-bin/private/ShowDocument?docid=5421


