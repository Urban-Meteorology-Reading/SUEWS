.. _SUEWS_Conductance.txt:

SUEWS_Conductance.txt
~~~~~~~~~~~~~~~~~~~~~

SUEWS_Conductance.txt contains the parameters needed for the Jarvis
(1976) [Ja76]_ surface conductance model used in the modelling of evaporation in
SUEWS. These values should **not** be changed independently of each
other. The suggested values below have been derived using datasets for
Los Angeles and Vancouver (see Järvi et al. (2011) [J11]_) and should be
used with `gsModel` = 1. An alternative formulation ( `gsModel` =2) uses
slightly different functional forms and different coefficients (with
different units).

.. DON'T manually modify the csv file below
.. as it is always automatically regenrated by each build

.. csv-table::
  :file: csv-table/SUEWS_Conductance.csv
  :header-rows: 1
  :widths: 5 25 5 65

.. only:: html

    An example `SUEWS_Conductance.txt` can be found below:

    .. literalinclude:: sample-table/SUEWS_Conductance.txt

.. only:: latex

    An example `SUEWS_Conductance.txt` can be found online
