========
desicoord
========

Introduction
------------

This package contains scripts and packages for coordinates transformations in DESI.

Example
------------

Here is an example::

    desi_fvc_proc.py -i fvc.20191125125304.fits  -o out.csv

Plot of the residuals with respect to the metrology::

    plot_fvc_residuals.py -i out.csv

Plot of the metrology data ::

    plot_metrology.py
