.. _SUEWSWUDAPT_Beijing:

Urban Energy Balance - SUEWS, WUDAPT and WATCH
==============================================

Introduction
------------

In this tutorial you will use UMEP to answer a research question. The
research question can be answered using SUEWS within the UMEP plugin. We
will go through each of the steps necessary to prepare, run and analyse
the model output.

Research question:
**What is the difference in sensible heat flux between a summer with a
large number of heat wave days and a low number of heat wave days in
Beijing, China?**

Getting started
---------------

Before getting started with this tutorial, make sure you have followed
these steps:

#. Install `QGIS <http://umep-docs.readthedocs.io/en/latest/Getting_Started.html>`__,
#. Install `UMEP <http://umep-docs.readthedocs.io/en/latest/Getting_Started.html>`__
#. Make sure the following python packages are installed: numpy,
   matplotlib and pandas. Start the WATCH-plugin (*UMEP > Pee-Processor > Meteorological Data > Download data (WATCH)*). If a message pops up that libaraies are missing, follow this `link <http://umep-docs.readthedocs.io/en/latest/Getting_Started.html#adding-missing-python-libraries-and-other-osgeo-functionalities>`__,.
#. Download and load the LCZ map from Beijing, available from the `WUDAPT portal <http://www.wudapt.org/>`__.

Overview
--------

.. figure::  /images/750px-Flowchart_beijing2.png
   :alt:  Flowchart of the workflow
   :width: 500px

   Flowchart of the workflow

#. **Extreme finder:** Find years with heat waves and years without
   waves
#. **Download data (WATCH):** Download data for the selected years.
#. **Vector grid:** Make a vector grid at the resolution for which you
   would like to run LQF and SUEWS.
#. **LCZ converter:** Convert Beijings LCZ map into morphometric
   parameters and Land cover fractions.
#. **Spatial downloader:** Downloads maps of population count and
   density.
#. **LQF:** Compute anthropogenic heat fluxes.
#. **SUEWS:** Calculates the energy balance.
#. **Benchmark:** Compares results of different runs.
#. **SUEWS analyser:** Displays the output for one single run.

Step 1: Extreme finder
----------------------

.. figure::  /images/600px-Extremefinder.png
   :alt: Screenshot of extreme finder
   :width: 100%

   Screenshot of the *Extreme Finder* tool

#. Open extreme finder at *UMEP > Processer > Outdoor Themal Comfort >
   ExtremeFinder*.
#. Manually enter coordinates of the location you are interested in, or
   click *fetch coordinates from map canvas*.
#. Select a period over which you would like to identify heat waves at
   *Start Date* and *End Date*.
#. Select a place and name for the output file.
#. Click *generate*.

Plots of maximum temperature, number of heat wave days, and box plots of
the maximum temperature of heat wave days over the selected year should
appear.

We choose 2009 as a year with a heat wave, due to the long heat wave
event at the end of June/start July. On the other hand, 2006 was
selected as a non-heat wave year.

Step 2: Download WATCH data
---------------------------

.. figure::  /images/600px-Watch.png
   :alt: Screenshot of Download Data (WATCH)

   Screenshot of Download Data (WATCH)

#. Open the WATCH data downloader at *UMEP > Pre-Processer >
   Meteorological Data > Download Data (WATCH).*
#. Click on *Fetch coordinates from map canvas* and click in the centre
   of the LCZ map, this will make the chosen coordinates show up on in
   the Latitude and Longitude boxes.
#. Specify the hours offset from UTC, for Beijing this is 8.
#. Specify the terrain height of the chosen coordinates. In the centre
   of Beijing this is about 50 meters.
#. Under *Compressed source data* specify a folder where the data should
   be downloaded.
#. Leave *Path to the AH results (optional)* blank for now.
#. Specify the dates for which meteorological data should be downloaded
#. At *Extract data to* specify the meteorological data text file.
#. Click *Generate*

This will take some minutes if you are downloading a year. Finally, you
should have a text file with the meteorological data and netcdf files of
each of the individual variables (Incoming shortwave and longwave
radiation, pressure, rain, temperature and humidity).

Step 3: Vector grid
-------------------

.. figure::  /images/450px-Vector.png
   :alt: vector.png

   vector.png

#. Open vector grid at *Vector > Research Tools > Vector grid*.
#. Select the extend of your LCZ map by clicking the ... next to *Grid
   extent (xmin, xmax, ymin, ymax)* and select *Use layer/canvas
   extent*.
#. Select the LCZ layer.
#. Specify the desired grid spacing, depending on the projection this
   will either be in meters or in degrees!
#. Make sure the output is in polygons, not lines.
#. Save the grid to a new layer.

Step 4: LCZ converter
---------------------

.. figure::  /images/600px-LCZdialog1.png
   :alt: Screenshot of LCZ converter

   Screenshot of LCZ converter

#. Open the LCZ converter at *UMEP > Pre-Processer > Spatial data > LCZ
   converter*.
#. Select the LCZ raster layer at '' LCZ raster''.
#. Select the vector grid you have just created in step 3 at *Vector
   grid* and select the ID field of the polygon grid at *ID field*.
#. By clicking *Adjust default parameters* you can edit the table. This
   table specifies the pervious, trees, grass, etc. fractions for each
   of the LCZ classes. For more information about each of the classes
   see `LCZConverter <http://umep-docs.readthedocs.io/en/latest/pre-processor/Spatial%20Data%20LCZ%20Converter.html>`__.
   If you choose to edit the table, make sure all fractions add up to
   1.0.
#. If you are unsure about the exact fractions for each of the LCZ click
   the tab *Pervious distribution*. Select *Same for all LCZ's*
#. Now you can select your best estimate about the distribution of the
   pervious surface fractions for urban and the tree distribution for
   rural. In addition, also specify the expected height of the trees.
#. Once you are satisfied click *Update Table*.
#. Select add results to polygon.
#. Add a file prefix if desired.
#. Finally select an output folder where you would like to receive the
   text files and click *Run*.

This should generate 3 text files, one with the land cover fractions,
one with morphometric parameters for buildings and one for trees for
each grid cell of the polygon grid.

Step 5: Spatial downloader
--------------------------

.. figure:: /images/600px-Spatialdownloader.png
   :alt: Spatialdownloader

   Spatialdownloader

In order to run LQF you will need
population counts for each of the grid cells you are modelling.

#. Open de spatial downloader at *UMEP > Pre-Processer > Spatial data >
   Spatial Data Downloader*.
#. Select *population density* and select the *GPWv4: UN-Adjusted
   Population Density* closest to the year you intend to model.
#. Make sure your canvas is zoomed out to the entire LCZ map and click
   *Use canvas extent*
#. Now click *Get data*.

You should get a raster of population density. These raster values will
need to be added to your vector grid by following `these
instructions <http://umep-docs.readthedocs.io/en/latest/OtherManuals/LQF_Manual.html#appendix-a-converting-a-population-raster-to-a-vector-shapefile-using-qgis>`__.
Finally the population densities need to be converted into population
counts:

#. Right-mouse click on your vector grid and click *Open Attribute
   Table*.
#. Click the abacus shaped symbol this is the Field calculator.
#. Under *Output field name* write "Pop, the *Output field type* should
   be “Decimal number (real)”, and the *Output Precision* can be set to
   2.
#. In the expression dialog box write popdens*$area/1000000, here
   popdens is the name of your population density field, $area
   is the size of the area for each grid cell and the 1 000 000 is to
   convert the data from m\ :sup:`2` to ha.
#. Click *OK* and you should have a new field called “Pop”.

Step 6: LQF
-----------


.. figure::  /images/LQf.png
   :alt: Screenshot of LQf

   Screenshot of LQf

Before running LQf, you will need to prepare some of
the data required to run it.

#. Convert the hourly temperatures you downloaded in step 2 into daily
   averaged temperatures in Excel, or a programming language of your
   choice.
#. Save the daily mean temperatures as a csv file with the first column
   the day of the year and the second column the temperature. The header
   and the data should look like:
   ::
     Data,T_Celsius
     StartDate,2006-1-1
     EndDate,2006-12-31
     Timezone,Asia/Shanghai
     1,-0.255517391
     2,-0.303882609
     3,-2.570373913
     4,-7.982847826
     5,-7.119765217
     6,-0.255517391

Prepare the
`database.nml <http:/www.urban-climate.net/umep/LQF_Manual#Data_sources_file>`__
and the
`modelparm.nml <http:/www.urban-climate.net/umep/LQF_Manual#Parameters_file>`__
as written in the manual. Make sure you change timezone,
use_uk_holidays, use_custom_holidays and custom_holidays in
modelparm.nml.

#. Open LQf at: *UMEP > Processer > Urban Energy Balance > Anthropogenic
   Heat LQf (LUCY)*.
#. Select the locations of the created modelparm.nml and database.nml
   and the folder you would to save the output to.
#. The extra disaggregation of the input data is optional, but the user
   could specify the land cover fractions generated in step 4 and the
   vector grid.
#. Click *Prepare input data using Data sources*. This may take a
   minute.
#. Once this first process has finished specify the *Date range* *Start
   date* and *End date*.
#. Click *Run model* and the model will take a while to run. If you are
   simulating an entire year this process may take a few hours.

Finally you should have csv files of anthropogenic heat fluxes for each
hour in the date range and for each grid cell of the vector grid.

Step 2b: Download WATCH data
----------------------------


.. figure::  /images/600px-Watch.png
   :alt: Screenshot of LQf

   Screenshot of Download Data (WATCH）

Inconsistant

#. Open the WATCH data downloader at *UMEP > Pre-Processer >
   Meteorological Data > Download Data (WATCH).*
#. Click on *Fetch coordinates from map canvas* and click in the centre
   of the LCZ map, this will make the chosen coordinates show up on in
   the Latitude and Longitude boxes.
#. Specify the hours offset from UTC, for Beijing this is 8.
#. Specify the terrain height of the chosen coordinates. In the centre
   of Beijing this is about 50 meters.
#. Under *Compressed source data* specify a folder where the data should
   be downloaded.
#. **Under *Path to the AH results (optional)* specify the folder the
   results from step 6 are saved to.**
#. Specify the dates for which meteorological data should be downloaded
#. At *Extract data to* specify the meteorological data text file.
#. Click *Generate*

This will take some minutes if you are downloading a year. Finally, you
should have a text file with the meteorological data and netcdf files of
each of the individual variables (Incoming shortwave and longwave
radiation, pressure, rain, temperature and humidity).

Step 7: SUEWS
-------------



.. figure::  /images/600px-Suews_sc.png
   :alt: Screenshot of SUEWS advance

   Screenshot of SUEWS advance

Before running SUEWS, you will need to
prepare some of the data required to run it.

#. Open SUEWS prepare at: *UMEP > Pre-Processer > SUEWS prepare*.
#. Under *vector polygon grid* specify your created vector grid and the
   *ID field*.
#. Select the location of the *Meteorological file* from step 2, the
   *Building morphology*, *Tree morphology* and *Land cover fractions*
   from step 4 and the population **density** from step 5.
#. Enter the start and end of day light savings time and the UTC offset.
#. Specify the *Leaf cycle* = winter when initialising in January.
   Unless the user has better information initialise the *Soil moisture
   state* at 100 %.
#. Select an output folder where the initial data to run SUEWS should be
   saved and press *Generate*.

When using LQf output as input for SUEWS there will be different
meteorological file for each grid cell of the vector grid. these files
all need to be renamed to the following format be7_2006_data_60.txt,
where 7 is the number of the grid cell. Rename all meteorological files
to this format.

#. Open SUEWS at *UMEP > Processer > Urban Energy Balance > Urban Energy
   Balance (SUEWS/BLUEWS, advanced).*
#. For the *Anthropogenic heat flux*, select option “[0] Observed data”.
   Feel free to change any other options.
#. Specify the *Temporal resolution of forcing data (minutes)* to be 60
   minutes.
#. Specify the *Temporal resolution of output (minutes)* to be 60
   minutes.
#. Select the *Input folder* specified in SUEWS prepare and select an
   output folder for the SUEWS output to be saved in. Finally, click
   *Run*.

This process will take several hours dependent how many grid cells are
used in the simulation. If the simulation is successful the output
folder will contain txt files with SUEWS output for each of the grid
cells in the vector grid.

Step 8: Benchmarking
--------------------

.. figure::  /images/600px-Benchmark.png
   :alt: Screenshot of LCZ converter

   Screenshot of LCZ converter

This system allows for comparison of runs
with observed data.

#. Open the Benchmarking system at: *UMEP > Post-Processer >
   Benchmarking system.*
#. Import data from the observation in the *Import of base dataset*.
   Specify the the number of rows in the header and the column
   separator. Note the names of the variables in the observational
   dataset should be the same as that of the SUEWS output
#. Select the *First comparison dataset (reference)* by pressing *Import
   data*' and importing the default SUEWS run at the location of the
   observations loaded in the previous step.
#. It is possible to load another dataset by checking *Add another
   comparison dataset* and selecting another SUEWS run or a different
   grid cell.
#. Once you have selected the appropriate datasets create a PDF at
   *Specify output PDF* and pressing *Run*.

This should generate a PDF of statistics for each variable with the
overall performance score, the mean absolute error, mean bias error and
the root mean squared error.

Step 9: SUEWS analyser
----------------------


.. figure::  /images/1200px-Suewsana.png
   :alt: Screenshot of SUEWS analyser
   :width: 100%

   Screenshot of *SUEWS Analyser* tool

This system allows for plotting of SUEWS
output.

#. Open the SUEWS analyser at: *UMEP > Post-Processer > Urban Energy
   Balance > SUEWS analyser.*
#. At *SUEWS RunControl namelist* select SUEWS RunControl.nml created
   with SUEWS prepare in step 7
#. Plot some of the basic data, such as the radiation and energy balance
   and the soil moisture and precipitation by selecting a grid cell on
   the upper left hand side. In addition select the *Year to
   investigate*. Check the box *Plot basic data* and press *Plot*. This
   should create a plot comparable to that in the screenshot above.

Plotting for example the mean daytime sensible heat flux during June
could be done as follows:

#. On the right-hand side of the dialog select “Sensible heat flux” at
   *Variable to analyse*.
#. Select the *Year to investigate* and the days of the year. June is
   DOY 152 - 181.
#. Select *Average* and *Only daytime*.
#. Finally select the vector grid created in step 3 and ID field at
   *Vector polygon grid used in the SUEWS model*, check *Add result to
   polygon grid* and click *Generate*.
#. This should generate an additional field in the attribute table of
   your vector grid. If it does not show up in the attribute table,
   reopen the vector grid.
#. In order to visualise the mean June daytime sensible heat flux, right
   click the vector grid in the layer panel and select *Properties*.
#. Go to *Style* and select *Graduated* in the top box. And select the
   QH column. Under *Color ramp* select the colour bar you prefer and
   click *Ok*
