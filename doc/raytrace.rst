=======================================
Ray trace runs
=======================================

The table ``desimeter/py/desimeter/data/raytrace-tan2fp-4957-v17.csv`` has been obtained by running::

  cd desimeter/py/desimeter/transform/tan2fp/raytrace
  ./run.py

The code in this directory has been copied from `DESI-4957`_, currently v17
  
The ``tan2fp`` transform parameters saved in ``desimeter/py/desimeter/data/raytrace-tan2fp.json`` have been obtained by running::

  fit_raytrace_tan2fp -i raytrace-tan2fp-4957-v17.csv -o raytrace-tan2fp.json --plot

.. _`DESI-4957`: https://desi.lbl.gov/DocDB/cgi-bin/private/ShowDocument?docid=4957

