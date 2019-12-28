========
desimeter
========

Introduction
------------

This package contains scripts and data for coordinates transformations in DESI. It's all in python, and runs on a laptop.

It comprises

* spots detection (based on a 2D FFT convolution of the input FVC image by Gaussian).
* spots centroid (simultaneous chi2 fit of the amplitude and position of a Gaussian integrated in the pixels).
* spots match to existing fiducials.
* fit of transformation from FVC pixels to focal plane X Y coordinates.
* metrology data, with all the routines to convert the engineering data in DocDB to the files used to fit the transformation, including a patch that corrects for missing or erroneous metrology, see py/desimeter/data/README.rst for more information on this.

Script Examples
---------------

Here is an example::

    desi_fvc_proc.py -i fvc.20191125125304.fits  -o out.csv

Plot of the residuals with respect to the metrology::

    plot_fvc_residuals.py -i out.csv

Plot of the metrology data ::

    plot_metrology.py

Code transform examples
-----------------------

Loading spots and fitting a FVC -> FP transform::

    from astropy.table import Table
    from desimeter.transform.fvc2fp.poly2d import FVCFP_Polynomial
    spots = Table.read('spots.csv')
    tx = FVCFP_Polynomial()
    tx.fit(spots)

Save that transform for later use::

    tx.write_jsonfile('fvc2fp.json')

Read it back in and do transforms between FVC and FP::

    t2 = FVCFP_Polynomial.read_jsonfile('fvc2fp.json')

    import numpy as np
    xpix, ypix = np.random.uniform(1000,5000, size=(2,50))
    xfp, yfp = t2.fvc2fp(xpix, ypix)

    xpix2, ypix2 = t2.fp2fvc(xfp, yfp)
    dr = np.sqrt((xpix2-xpix)**2 + (ypix2-ypix)**2)
    print(np.median(dr))

Note: class names and module organization will change.

