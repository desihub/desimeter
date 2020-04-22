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
* matching positioner spots based upon their expected locations.
* sky transformations (precession, aberration, refraction, polar mis-alignment)
* fit of GFA guide stars to adjust telescope pointing, field rotation while absorbing a scaling factor (potentially accounting for thermal dilatation)
* prediction of fiber sky RA,Dec coordinates

See `desimeter/doc/coordinates.rst` for a definition of the coordinate systems.


Automated build & testing
-------------------------

.. image:: https://img.shields.io/circleci/project/github/desihub/desimeter.svg
  :alt: Build
  :target: https://circleci.com/gh/desihub/desimeter

.. image:: https://coveralls.io/repos/github/desihub/desimeter/badge.svg?branch=master
  :alt: Coverage
  :target: https://coveralls.io/github/desihub/desimeter?branch=master
  
Script Examples
---------------

These examples use a fiber view camera (FVC) image taken during PSF stability
tests, available at NERSC at::
     /global/project/projectdirs/desi/cmx/psfstability/thetaphi-20191113/data/fvc/fvc.20191113120837.fits

They also work on FVC images written by ICS as part of platemaker acquisition
sequences by using the images in the ``F0000``, ``F0001``, etc FITS extensions.

Here is an example fitting FVC spots and matching to fiducials::

    desi_fvc_proc -i fvc.20191113120837.fits  -o spots.csv

Plot of the residuals with respect to the metrology::

    plot_fvc_residuals -i out.csv

Plot of the metrology data ::

    plot_metrology

Fit of guide stars (at KPNO), saving the transformation in a json file ::
  
     desi_fit_guide_star_coordinates -i /n/home/datasystems/users/ameisner/reduced/realtime/20200119/00042312/gfa-00042312_catalog.fits --fits-header /exposures/desi/20200119/00042312/gfa-00042312.fits.fz -o fp-00042312.json  

Then using it to predict the fiber RA,Dec from a FVC image ::

     desi_fvc_proc -i /exposures/desi/20200119/00042311/fvc-00042311.fits.fz --field-model fp-00042312.json --expected /exposures/desi/20200119/00042311/coordinates-00042311.fits -o fvc-00042311.csv

Code examples
-------------

Detect spots and match fiducials::

    import fitsio
    from desimeter.detectspots import detectspots
    from desimeter.findfiducials import findfiducials
    image = fitsio.read('fvc.20191113120837.fits')

    #- Detects spots: XPIX, YPIX, XERR, YERR, COUNTS
    spots = detectspots(image)
    print(spots.colnames)

    #- matches to fiducials and adds LOCATION and PINHOLE_ID columns
    spots = findfiducials(spots)
    print(spots.colnames)
    spots.write('spots.csv', overwrite=True)

Load spots and fit a FVC -> FP transform::

    from astropy.table import Table
    from desimeter.transform import fvc2fp
    spots = Table.read('spots.csv')

    #- Fit transform; adds X_FP, Y_FP, X_FP_METRO, Y_FP_METRO columns
    tx = fvc2fp.fit(spots, update_spots=True)
    print(spots.colnames)

Check consistency with fiducial metrology (PINHOLE_ID>0)::

    ii = spots['PINHOLE_ID'] > 0
    dr2 = (spots['X_FP'] - spots['X_FP_METRO'])**2 + \
         (spots['Y_FP'] - spots['Y_FP'])**2
    rms = np.sqrt(np.mean(dr2[ii]))
    print('1D RMS = {:.1f} um'.format(1000*rms))

Save that transform for later use::

    tx.write_jsonfile('fvc2fp.json')

Read it back in and do transforms between FVC and FP::

    t2 = fvc2fp.read_jsonfile('fvc2fp.json')

    import numpy as np
    xpix, ypix = np.random.uniform(1000,5000, size=(2,50))
    xfp, yfp = t2.fvc2fp(xpix, ypix)

    xpix2, ypix2 = t2.fp2fvc(xfp, yfp)
    dr = np.sqrt((xpix2-xpix)**2 + (ypix2-ypix)**2)
    print(np.median(dr))

Dependencies
------------

desimeter requires numpy, scipy, astropy, fitsio, and matplotlib.

Python 3.6 or greater is required.

It purposefully does *not* require desiutil, desimodel, or any other
offline desidata packages to facilitate integration with the DESI online
environment and to minimize getting started overhead for non-desidata users.

Similarly, it does *not* use the ICS ops database or any online code to
facilitate offline development and studies, e.g. on a laptop.

Installation
------------

Get a copy of the code::

    git clone https://github.com/desihub/desimeter

If you want to use desimeter but don't intend to actively develop it::

    cd desimeter
    python setup.py install

For developers, we recommend adding `desimeter/py` to `$PYTHONPATH`
and `desimeter/bin` to `$PATH` instead of installing desimeter.

Other Notes
-----------

desimeter is a work in progress and we expect that class names and module
organization will change.

