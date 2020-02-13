=======================================
Coordinates
=======================================

This document presents the different coordinate systems and the coordinate transformation routines in desimeter.

Coordinate systems
++++++++++++++++++

GFA (Guide Focus Array)
~~~~~~~~~~~~~~~~~~~~~~~

Pixel coordinates in GFA guide CCD images. Usually called ``x_gfa , y_gfa`` in desimeter.
Measurements are done by Aaron Meisner's code.

FVC (Fiber View Camera)
~~~~~~~~~~~~~~~~~~~~~~~

Pixel coordinates in the Fiber View Camera CCD images. Usually called ``x_fvc , y_fvc`` in desimeter.
In this code, the coordinates of the first pixel is ``(0,0)``. This is different from other codes. 

PTL (Petal)
~~~~~~~~~~~

Petal coordinates in mm, either as cartesian coordinates ``x_ptl,y_ptl,z_ptl`` or curved coordinates ``s_ptl,q_ptl``. This coordinate system is attached to the metrology measurements for positioner mounting holes, FIF and GIF. See `data.README`_ .

FP (Focal Plane)
~~~~~~~~~~~~~~~~

Coordinate system as close as possible to the 'CS5' system described in ...
PTL and FP are related by a translation and a rotation. ``x_fp,y_fp``  are also in mm. It is not clear yet whether the petals are moving one with respect to the others.

TP (Tangent Plane)
~~~~~~~~~~~~~~~~~~


Polar coordinates (theta,phi) projected on a plane with an arbitrary but fixed projection::

  x_tan = sin(theta) cos(phi)
  y_tan = sin(theta) sin(phi)

For a telescope pointing at ``HA_tel , Dec_tel``,
``phi=0`` for a star with ``HA>HA_tel`` and ``Dec=Dec_tel``, and ``phi=90 deg`` for a star with ``HA=HA_tel`` and ``Dec>Dec_tel``.

The angles (theta,phi) are the direction of photons hitting the telescope (so including refraction but excluding optical distortions) in a frame attached to the focal plane (the frame is rotating with the hexapod).

ICRS (International Celestial Reference System)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``RA,Dec`` sky coordinates of targets (same system as GAIA stars used for the astrometry).


Coordinate transformations
++++++++++++++++++++++++++

PTL -> FP
~~~~~~~~~

Solid transformation from PTL to FP composed of a translation and a rotation.
The transformation coefficients are saved in the yaml file `petal-alignments.yaml`_ . In the current
version, they are the same as the ones used by PlateMaker. They were dowloaded from the online
data base using the script `load_petal_alignments_from_db`_ .


GFA -> FP
~~~~~~~~~

Based on the metrology compilation in `DESI-5421`_ + the transformation from PTL to FP.
The fit of the metrology is performed on the fly for the first call of the function.

Example::

  from desimeter.transform.gfa2fp import gfa2fp,fp2gfa
  x_fp,y_fp   = gfa2fp(petal_loc,x_gfa,y_gfa)
  x_gfa,y_gfa = fp2gfa(petal_loc,x_fp,y_fp)

FVC -> FP
~~~~~~~~~

Based on a fit of the measured pixel coordinates of fiducials and their metrology.
In most and probably all cases, the transformation is used in the same procedure
as the fit (see the script `desi_fvc_proc`_ ). Two implementations of this
transform exist in `transform.fvc2fp`_ , but we used by default the Zhao-Burge transform
implemented in `transform.fvc2fp.zb.py`_ .

Example::

  from desimeter.transform.fvc2fp.zb import FVCFP_ZhaoBurge
  tx = FVCFP_ZhaoBurge()
  tx.fit(spots, update_spots=True)
  tx.write_jsonfile(filename)
  x_fp,y_fp = tx.fvc2fp(x_fvc,y_fvc)
  
FP -> TAN
~~~~~~~~~

From the focal surface instruments to the tangent plane. We have two implementations
for now.

* `transform.tan2fp.echo22.py`_ It is based on the echo22 design model. Does not account for ADC rotation.
* `transform.tan2fp.raytracefit.py`_ It is based on a ray tracing model fitted fitted ZB polynomials. It does include ADC rotation and is now the default. More information on this in `raytrace.rst`_ . 
 
Example::

   from desimeter.transform.tan2fp import fp2tan,tan2fp
   x_tan,y_tan = fp2tan(x_fp,y_fp,adc1=60.,adc2=120.)
   x_fp,y_fp = fp2tan(x_tan,y_tan,adc1=60.,adc2=120.)
   
   
TAN -> ICRS
~~~~~~~~~~~

This transformation includes precession, aberration, refraction, polar mis-alignment, hexapod rotation angle.
The code is in `transform.radec2tan.py`_ .

Example::

   from desimeter.transform.radec2tan import radec2tan,tan2radec
   ra,dec = radec2tan(x_tan,y_tan,telescope_ra,telescop_dec,mjd,lst,precession=True,aberration=True)

Field Model
~~~~~~~~~~~

This is not a transformation from one system to another but rather the combination of the transformations,
including a correction based on the location of GFA guide stars. This "field model" is fit to the data
using the script `desi_fit_guide_star_coordinates`_ . The routines are in `fieldmodel.py`_ .

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
.. _`petal-alignments.yaml`: https://github.com/desihub/desimeter/blob/master/py/desimeter/data/petal-alignments.yaml
.. _`load_petal_alignments_from_db`: https://github.com/desihub/desimeter/blob/master/bin/load_petal_alignments_from_db
.. _`desi_fvc_proc`: https://github.com/desihub/desimeter/blob/master/bin/desi_fvc_proc
.. _`transform.fvc2fp`: https://github.com/desihub/desimeter/tree/master/py/desimeter/transform/fvc2fp
.. _`transform.fvc2fp.zb.py`: https://github.com/desihub/desimeter/blob/master/py/desimeter/transform/fvc2fp/zb.py
.. _`transform.tan2fp.echo22.py`: https://github.com/desihub/desimeter/blob/master/py/desimeter/transform/tan2fp/echo22.py
.. _`transform.tan2fp.raytracefit.py`: https://github.com/desihub/desimeter/blob/master/py/desimeter/transform/tan2fp/raytracefit.py
.. _`raytrace.rst`: https://github.com/desihub/desimeter/blob/master/doc/raytrace.rst
.. _`transform.radec2tan.py`: https://github.com/desihub/desimeter/blob/master/py/desimeter/transform/radec2tan.py
.. _`desi_fit_guide_star_coordinates`: https://github.com/desihub/desimeter/blob/master/bin/desi_fit_guide_star_coordinates
.. _`fieldmodel.py`: https://github.com/desihub/desimeter/blob/master/py/desimeter/fieldmodel.py

