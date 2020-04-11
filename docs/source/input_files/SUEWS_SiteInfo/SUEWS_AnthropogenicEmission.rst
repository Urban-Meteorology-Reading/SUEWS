.. SUEWS_AnthropogenicEmission.txt:

SUEWS_AnthropogenicEmission.txt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note::
    this file used to be named as ``SUEWS_AnthropogenicHeat.txt``
    and is changed to this name in v2019a.


SUEWS_AnthropogenicEmission.txt provides the parameters needed to model
the anthropogenic heat flux using either the method of Järvi et al.
(2011) based on heating and cooling degree days (`EmissionsMethod` = 2
in `RunControl.nml`) or the method of Loridan et
al. (2011) based on air temperature (`EmissionsMethod` = 1 in
`RunControl.nml`).


The sub-daily variation in
anthropogenic heat flux is modelled according to the daily cycles
specified in SUEWS_Profiles.txt.

Alternatively, if available, the anthropogenic heat flux can be provided in the met forcing file (and set `EmissionsMethod` = 0 in `RunControl.nml`) by filling the :option:`qf` column with valid values.

.. DON'T manually modify the csv file below
.. as it is always automatically regenrated by each build:
.. edit the item descriptions in file `Input_Options.rst`

.. csv-table::
  :file: csv-table/SUEWS_AnthropogenicEmission.csv
  :header-rows: 1
  :widths: 5 25 5 65

.. only:: html

    An example `SUEWS_AnthropogenicEmission.txt` can be found below:

    .. literalinclude:: sample-table/SUEWS_AnthropogenicEmission.txt

.. only:: latex

    An example `SUEWS_AnthropogenicEmission.txt` can be found in the online version.
