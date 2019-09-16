.. _SUEWS_Irrigation.txt:

SUEWS_Irrigation.txt
~~~~~~~~~~~~~~~~~~~~

SUEWS includes a simple model for external water use if observed data
are not available. The model calculates daily water use from the mean
daily air temperature, number of days since rain and fraction of
irrigated area using automatic/manual irrigation. The sub-daily pattern
of water use is modelled according to the daily cycles specified in
`SUEWS_Profiles.txt`.

Alternatively, if available, the external water use can be provided in
the met forcing file (and set `WaterUseMethod` = 1 in
`RunControl.nml`) by filling the `Wuh` columns with valid values.

.. DON'T manually modify the csv file below
.. as it is always automatically regenrated by each build:
.. edit the item descriptions in file `Input_Options.rst`

.. csv-table::
  :file: csv-table/SUEWS_Irrigation.csv
  :header-rows: 1
  :widths: 5 25 5 65

.. only:: html

    An example `SUEWS_Irrigation.txt` can be found below:

    .. literalinclude:: sample-table/SUEWS_Irrigation.txt

.. only:: latex

    An example `SUEWS_Irrigation.txt` can be found in the online version.
