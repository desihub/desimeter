========
desimeter
========

Introduction
------------

This package contains scripts and packages for coordinates transformations in DESI. It's all in python.
It comprises

* spots detection (based on a 2D FFT convolution of the input FVC image by Gaussian).
* spots centroid (simultaneous chi2 fit of the amplitude and position of a Gaussian integrated in the pixels).
* spots match to existing fiducials.
* fit of transformation from FVC pixels to focal plane X Y coordinates.
* metrology data, with all the routines to convert the engineering data in DocDB to the files used to fit the transformation, including a patch that corrects for missing or erroneous metrology, see py/desimeter/data/README.rst for more information on this.

Example
------------

Here is an example::

    desi_fvc_proc.py -i fvc.20191125125304.fits  -o out.csv

Plot of the residuals with respect to the metrology::

    plot_fvc_residuals.py -i out.csv

Plot of the metrology data ::

    plot_metrology.py
