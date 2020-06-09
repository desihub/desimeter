=======================================
Extracting posmoves from DB and fitting
positioner parameters
=======================================

1) Read posmoves from the data base and write them as cvs ASCI tables.
   (one file per positioner)

Example::

  get_posmoves --host beyonce.lbl.gov --port 5432 --password reader --petal-ids 1 --exposure-ids 2107 --outdir . --pos-id M02155

2) Fit the positioner calibration

Example::

  fit_posparams -i M02155.csv -o .

3) Plot the moves along with the fitted positioner parameters

Example::

  plot_posmoves -i M02155.csv -p M02155_paramfits.csv

=======================================
Deprecated method
to match posmov db entries and FVC data
=======================================

1) Make a list of FVC images.

Example::

  find /global/cfs/cdirs/desi/engineering/focalplane/logs/kpno -name \*fvc\*.fits\* > /global/cfs/cdirs/desi/engineering/fvc/desimeter/v0.2.0/fvc-images.txt
  find /global/cfs/cdirs/desi/spectro/data -name \*fvc\*.fits\* >> /global/cfs/cdirs/desi/engineering/fvc/desimeter/v0.2.0/fvc-images.txt

2) Write and run desimeter

Run::

  write_fvc_image_processing_scripts -i  /global/cfs/cdirs/desi/engineering/fvc/desimeter/v0.2.0/fvc-images.txt --jobdir ./jobs --outdir /global/cfs/cdirs/desi/engineering/fvc/desimeter/v0.2.0


Then execute the jobs written in the job directory.
This script writes the file ``/global/cfs/cdirs/desi/engineering/fvc/desimeter/v0.2.0/fvc-desimeter-and-dates.csv`` which contains a list of desimeter FVC CSV tables with time stamps::

  DATE,FILE
  2019-10-07T20:31:49,/global/cfs/cdirs/desi/engineering/fvc/desimeter/v0.2.0/fvc-00017943-F0000-2019-10-07T20:31:49.csv
  2019-10-07T20:31:58,/global/cfs/cdirs/desi/engineering/fvc/desimeter/v0.2.0/fvc-00017943-F0001-2019-10-07T20:31:58.csv

3) Match to the posmove db::

   get_posmov_fvc_data --fvc-file-list /global/cfs/cdirs/desi/engineering/fvc/desimeter/v0.2.0/fvc-desimeter-and-dates.csv --outdir /global/cfs/cdirs/desi/engineering/fvc/desimeter/v0.2.0/positioners-3

This scripts matches each pos move db entry to FVC coordinates of positioner, for each positioner of each petal.
