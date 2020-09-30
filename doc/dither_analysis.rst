Example dither analysis
=================================

Copy, link or download files
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

From Eddie's dither analysis

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

``fm-00055670.json`` is the "field model" file with the astrometric solution::

 {"name": "Field model", "version": "1", "ra": 169.97172680058398, "dec": 60.001396953233744, "mjd": 58924.473182622765, "lst": 233.09726528124884, "hexrot_deg": 5.555555555555556e-05, "adc1": 51.980022, "adc2": 138.720115, "sxx": 0.9999537892790107, "syy": 0.9999774313698416, "sxy": -4.109184104764672e-05, "fieldrot_zp_deg": 359.946666865039, "fieldrot_deg": 0.06743416213752557, "expid": 55670, "nstars": 79, "rms_arcsec": 0.1640910115685826}

One can see the field dilatation is only of 5e-5. The astrometric residuals to 79 guide stars are of 0.16 arcsec.

Fit FVC image
+++++++++++++

Extract fiber coordinates from FVC image, convert to focal plane coordinates using fiducials, and use the field model derived above to convert the focal plane coordinates to sky coordinates:

 desi_fvc_proc -i fvc-00055670.fits.fz --extname F0002 --field-model fm-00055670.json --expected coordinates-00055670.fits -o fvc-00055670.csv --zbfit

(note the --zbfit option to absorb some of the distortions)

Display the result of the match between the measured fiducial locations and the metrology::

 plot_fvc_residuals -i fvc-00055670.csv


Compare with dither analysis
++++++++++++++++++++++++++++

Compare the desimeter prediction for the fiber to target RA,Dec offsets with the
results from the dither analysis::

 plot_dither_vs_desimeter  dither20200315-63224-B.fits fvc-00055670.csv fiberassign-063225.fits 55670


Here the desimeter to dither residuals are of 0.42 arcsec. After correcting for a rotation, shift, and dilatation,
the residuals are of 0.26 arcsec::

  number of valid targets from dither analysis = 1758
  residual scale along HA       = 1-7.8e-05
  residual scale along Dec      = 1-4.1e-06
  residual pointing offset dHA  = -0.11 arcsec
  residual pointing offset dDec = -0.15 arcsec
  residual rotation             = -9.19 arcsec
     desimeter(guide+fvc)-target rms 2D = 0.44 arcsec
                   dither-target rms 2D = 0.23 arcsec
     desimeter(guide+fvc)-dither rms 2D = 0.42 arcsec
  transformed(desimeter) -dither rms 2D = 0.26 arcsec
