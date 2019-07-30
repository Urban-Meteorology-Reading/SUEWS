.. _SUEWS_AnthropogenicHeat.txt:

SUEWS_AnthropogenicHeat.txt
~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUEWS_AnthropogenicHeatFlux.txt provides the parameters needed to model
the anthropogenic heat flux using either the method of Järvi et al.
(2011) based on heating and cooling degree days (`EmissionsMethod` = 2
in `RunControl.nml`) or the method of Loridan et
al. (2011) based on air temperature (`EmissionsMethod` = 1 in
`RunControl.nml`). The sub-daily variation in
anthropogenic heat flux is modelled according to the daily cycles
specified in SUEWS_Profiles.txt. 

Alternatively, if available, the anthropogenic heat flux can be provided in the met forcing file (and set `EmissionsMethod` = 0 in `RunControl.nml`) by filling the `qf` column with valid values.

.. DON'T manually modify the csv file below
.. as it is always automatically regenrated by each build

.. csv-table::
  :file: csv-table/SUEWS_AnthropogenicHeat.csv
  :header-rows: 1
  :widths: 5 25 5 65

.. only:: html

    An example `SUEWS_AnthropogenicHeat.txt` can be found below:

    .. literalinclude:: sample-table/SUEWS_AnthropogenicHeat.txt

.. only:: latex

    An example `SUEWS_AnthropogenicHeat.txt` can be found in the online version.
