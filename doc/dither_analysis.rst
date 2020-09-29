Example dither analysis
=================================

Copy, link or download files
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

From Eddie's dither analysis::
 https://portal.nersc.gov/project/cosmo/temp/schlafly/dither/dither20200315-63224-B.fits

Choose one exposure in the file (column 'expid' in first HDU), here the first one: 55670

From the exposure directory::
 /project/projectdirs/desi/spectro/data/20200315/00055670/guide-00055670.fits.fz
 /project/projectdirs/desi/spectro/data/20200315/00055670/fvc-00055670.fits.fz
 /project/projectdirs/desi/spectro/data/20200315/00055670/fiberassign-063225.fits
 /project/projectdirs/desi/spectro/data/20200315/00055670/coordinates-00055670.fits

From Aaron's pipeline::
 /project/projectdirs/desi/users/ameisner/GFA/reduced/v0005_guide/20200315/00055670/guide-00055670_catalog-00001.fits


Fit the astrometric solution
++++++++++++++++++++++++++++

Fit pointing offsets, the field rotation and a field compression (3 numbers: sxx,syy,sxy) using the sky coordinates and the CCD pixel coordinates of guide stars::

 desi_fit_guide_star_coordinates -i guide-00055670_catalog-00001.fits --fits-header guide-00055670.fits.fz -o fm-00055670.json

``fm-00055670.json`` is the "field model" file with the astometric solution::

 {"name": "Field model", "version": "1", "ra": 106.33481916831745, "dec": 31.889476531165776, "mjd": 58904.25116974065, "lst": 133.24085389263928, "hexrot_deg": 0.0, "adc1": 70.412607, "adc2": 95.190051, "sxx": 0.9999611908551399, "syy": 0.9999892189227471, "sxy": -4.799877494061412e-05, "fieldrot_zp_deg": 359.9706671491963, "fieldrot_deg": 0.14414265500548332, "expid": 55670, "nstars": 341, "rms_arcsec": 0.2659517231204168}

One can see the field dilatation is only of 4e-5.

Fit FVC image
+++++++++++++

Extract fiber coordinates from FVC image, convert to focal plane coordinates using fiducials, and use the field model derived above to convert the focal plane coordinates to sky coordinates::

 desi_fvc_proc -i fvc-00055670.fits.fz --extname F0002 --field-model fm-00055670.json --expected coordinates-00055670.fits -o fvc-00055670.csv

Display the result of the match between the measured fiducial locations and the metrology::

 plot_fvc_residuals -i fvc-00055670.csv


Compare with dither analysis
++++++++++++++++++++++++++++

Compare the desimeter prediction for the fiber to target RA,Dec offsets with the
results from the dither analysis::

 plot_dither_vs_desimeter  dither20200315-63224-B.fits fvc-00055670.csv fiberassign-063225.fits
