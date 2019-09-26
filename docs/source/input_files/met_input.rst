.. _met_input:

Meteorological Input File
-------------------------

SUEWS is designed to run using commonly measured meteorological
variables.

-  Required inputs must be continuous – i.e. **gap fill** any missing
   data.
-  Temporal information (i.e., ``iy``, ``id``, ``it`` and ``imin``
   should be in local time.
-  The table below gives the must-use (MU) and optional (O) additional
   input variables.
-  If an optional input variable is not available or will not be used by
   the model, enter ‘-999.0’ for this column.
-  Since v2017a forcing files no longer need to end with two rows
   containing ‘-9’ in the first column.

-  One single meteorological file can be used for all grids
   (**MultipleMetFiles=0** in `RunControl.nml`, no
   grid number in file name) if appropriate for the study area, or
-  separate met files can be used for each grid if data are available
   (**MultipleMetFiles=1** in `RunControl.nml`,
   filename includes grid number).

-  The meteorological forcing file names should be appended with the
   temporal resolution in minutes (SS_YYYY_data_tt.txt, or
   SSss_YYYY_data_tt.txt for multiple grids).

-  Separate met forcing files should be provided for each year.
-  Files do not need to start/end at the start/end of the year, but they
   must contain a whole number of days.
-  The meteorological input file should match the information given in
   `SUEWS_SiteSelect.txt`.
-  If a *partial year* is used that specific year must be given in
   SUEWS_SiteSelect.txt.
-  If *multiple years* are used, all years should be included in
   SUEWS_SiteSelect.txt.
-  If a *whole year* (e.g. 2011) is intended to be modelled using and
   hourly resolution dataset, the number of lines in the met data file
   should be 8760 and begin and end with::

     iy     id  it  imin
     2011   1   1   0 …
     …
     2012   1   0   0 …


.. _SSss_YYYY_data_tt.txt:

SSss_YYYY_data_tt.txt
~~~~~~~~~~~~~~~~~~~~~

Main meteorological data file.

.. csv-table::
  :file: SSss_YYYY_data_tt.csv
  :header-rows: 1
  :widths: auto
