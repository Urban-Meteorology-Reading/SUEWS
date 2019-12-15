.. _output_files:

Output files
============

Runtime diagnostic information
------------------------------

.. _problems.txt:

Error messages: problems.txt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If there are problems running the program serious error messages will be
written to problems.txt.

-  Serious problems will usually cause the program to stop after writing
   the error message. If this is the case, the last line of problems.txt
   will contain a non-zero number (the error code).
-  If the program runs successfully, problems.txt file ends with::

    Run completed.
    0

SUEWS has a large number of error messages included to try to capture
common errors to help the user determine what the problem is. If you
encounter an error that does not provide an error message please capture
the details so we can hopefully provide better error messages in future.

See `Troubleshooting` section for help solving
problems. If the file paths are not correct the program will return an
error when run (see `Preparing_to_run_the_model`).

.. _warnings.txt:

Warning messages: warnings.txt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  If the program encounters a more minor issue it will not stop but a
   warning may be written to warnings.txt. It is advisable to check the
   warnings to ensure there is not a more serious problem.
-  The warnings.txt file can be large (over several GBs) given warning
   messages are written out during a large scale simulation, you can use
   :code:`tail`/:code:`head` to view the ending/starting part without opening
   the whole file on Unix-like systems (Linux/mac OS), which may slow
   down your system.
-  To prevent warnings.txt from being written, set :option:`SuppressWarnings`
   to 1 in `RunControl.nml`.
-  Warning messages are usually written with a grid number, timestamp
   and error count. If the problem occurs in the initial stages (i.e.
   before grid numbers and timestamps are assigned, these are printed as
   00000).

Summary of model parameters: SS_FileChoices.txt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For each run, the model parameters specified in the input files are
written out to the file SS_FileChoices.txt.

Model output files
------------------

SSss_YYYY_SUEWS_TT.txt
~~~~~~~~~~~~~~~~~~~~~~

SUEWS produces the main output file (SSss_YYYY_SUEWS_tt.txt) with time
resolution (TT min) set by :option:`ResolutionFilesOut` in `RunControl.nml`.

Before these main data files are written out, SUEWS provides a summary
of the column names, units and variables included in the file
Ss_YYYY_TT_OutputFormat.txt (one file per run).

The variables included in the main output file are determined according
to :option:`WriteOutOption` set in :ref:`RunControl.nml`.


.. csv-table::
  :file: SSss_YYYY_SUEWS_TT.csv
  :header-rows: 1
  :widths: auto


SSss_DailyState.txt
~~~~~~~~~~~~~~~~~~~

Contains information about the state of the surface and soil and
vegetation parameters at a time resolution of one day. One file is
written for each grid so it may contain multiple years.

.. csv-table::
  :file: SSss_DailyState.csv
  :header-rows: 1
  :widths: auto

.. _initialconditionsssss_yyyy.nml:

InitialConditionsSSss_YYYY.nml
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

At the end of the model run (or the end of each year in the model run) a
new InitialConditions file is written out (to the input folder) for each
grid, see `Initial_Conditions`

SSss_YYYY_snow_TT.txt
~~~~~~~~~~~~~~~~~~~~~

SUEWS produces a separate output file for snow (when :option:`snowUse` = 1 in
`RunControl.nml`) with details for each surface type.

File format of SSss_YYYY_snow_TT.txt

.. csv-table::
  :file: SSss_YYYY_snow_TT.csv
  :header-rows: 1
  :widths: auto

SSss_YYYY_RSL_TT.txt
~~~~~~~~~~~~~~~~~~~~~

SUEWS produces a separate output file for wind, temperature and humidity
profiles in the roughness sublayer at 30 levels:
levels 1 and 30 are positioned at 0.1 and 3.0 ``Zh`` (i.e., canopy height)
with other levels evenly distributed in between.

File format of SSss_YYYY_RSL_TT.txt:

.. csv-table::
  :file: SSss_YYYY_RSL_TT.csv
  :header-rows: 1
  :widths: auto

SSss_YYYY_BL_TT.txt
~~~~~~~~~~~~~~~~~~~~

Meteorological variables modelled by CBL portion of the model are output
in to this file created for each day with time step (see section CBL
Input).

.. csv-table::
  :file: SSss_YYYY_BL_TT.csv
  :header-rows: 1
  :widths: auto


.. SOLWEIG is fully removed since 2019a

.. SOLWEIGpoiOut.txt
.. ~~~~~~~~~~~~~~~~~

.. Calculated variables from POI, point of interest (row, col) stated in
.. `SOLWEIGinput.nml`.

.. SOLWEIG model output file format: SOLWEIGpoiOUT.txt


.. .. csv-table::
..   :file: SOLWEIGpoiOut.csv
..   :header-rows: 1
..   :widths: auto



SSss_YYYY_ESTM_TT.txt
~~~~~~~~~~~~~~~~~~~~~

If the ESTM model option is run, the following output file is created.
**Note: First time steps of storage output could give NaN values during
the initial converging phase.**

ESTM output file format

.. csv-table::
  :file: SSss_YYYY_ESTM_TT.csv
  :header-rows: 1
  :widths: auto
