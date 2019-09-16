.. _SUEWS_Profiles.txt:

SUEWS_Profiles.txt
~~~~~~~~~~~~~~~~~~

SUEWS_Profiles.txt specifies the daily cycle of variables related to
human behaviour (energy use, water use and snow clearing). Different
profiles can be specified for weekdays and weekends. The profiles are
provided at hourly resolution here; the model will then interpolate the
hourly energy and water use profiles to the resolution of the model time
step and normalize the values provided. Thus it does not matter whether
columns 2-25 add up to, say 1, 24, or another number, because the model
will handle this. Currently, the snow clearing profiles are not
interpolated as these are effectively a switch (0 or 1).

If the anthropogenic heat flux and water use are specified in the met
forcing file, the energy and water use profiles are not used.

Profiles are specified for the following

-  Anthropogenic heat flux (weekday and weekend)
-  Water use (weekday and weekend; manual and automatic irrigation)
-  Snow removal (weekday and weekend)
-  Human activity (weekday and weekend).


.. DON'T manually modify the csv file below
.. as it is always automatically regenrated by each build:
.. edit the item descriptions in file `Input_Options.rst`

.. csv-table::
  :file: csv-table/SUEWS_Profiles.csv
  :header-rows: 1
  :widths: 5 25 5 65

.. only:: html

    An example `SUEWS_Profiles.txt` can be found below:

    .. literalinclude:: sample-table/SUEWS_Profiles.txt

.. only:: latex

    An example `SUEWS_Profiles.txt` can be found in the online version.
