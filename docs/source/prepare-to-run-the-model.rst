.. _Preparing_to_run_the_model:

Preparing to run the model
==========================

The following is to help with the model setup. Note that there are also
starting `tutorials`_  for the
version of SUEWS in `UMEP`_.
The version there is the same (i.e. the executable) as the
standalone version so you can swap to that later once you have some
familiarity.

Preparatory reading
-------------------

Read the manual and relevant papers (and references therein):

-  Järvi L, Grimmond CSB & Christen A (2011) The Surface Urban Energy
   and Water Balance Scheme (SUEWS): Evaluation in Los Angeles and
   Vancouver. J. Hydrol. 411, 219-237.
   `doi:10.1016/j.jhydrol.2011.10.00 <http://www.sciencedirect.com/science/article/pii/S0022169411006937>`__
-  Järvi L, Grimmond CSB, Taka M, Nordbo A, Setälä H & Strachan IB
   (2014) Development of the Surface Urban Energy and Water balance
   Scheme (SUEWS) for cold climate cities. Geosci. Model Dev. 7,
   1691-1711.
   `doi:10.5194/gmd-7-1691-2014 <http://www.geosci-model-dev.net/7/1691/2014/>`__
-  Ward HC, Kotthaus S, Järvi L and Grimmond CSB (2016) Surface Urban
   Energy and Water Balance Scheme (SUEWS): development and evaluation
   at two UK sites. Urban Climate 18, 1-32.
   `doi:10.1016/j.uclim.2016.05.001 <http://www.sciencedirect.com/science/article/pii/S2212095516300256/>`__

`See other publications with example applications <Recent_publications>`

Decide what type of model run you are interested in
---------------------------------------------------

.. list-table::
   :widths: 50 50
   :header-rows: 1

   * -
     - Available in this release
   * - LUMPS
     - Yes – not standalone
   * - SUEWS at a point or for an individual area
     - Yes
   * - SUEWS for multiple grids or areas
     - Yes
   * - SUEWS with Boundary Layer (BL)
     - Yes
   * - SUEWS with snow
     - Yes
   * - SUEWS with SOLWEIG
     - No
   * - SUEWS with SOLWEIG and BL
     - No

Download the program and example data files
-------------------------------------------

Visit the `website <https://urban-meteorology-reading.github.io/SUEWS>`_ to receive a link to download the program and example
data files. Select the appropriate compiled version of the model to download. For windows there is an installation version which will put the programs and all the files into the appropriate place. There is also a version linked to QGIS:
`UMEP`_.

Note, as the definition of long double precision varies between
computers (e.g. Mac vs Windows) slightly different results may occur in
the output files.

Test/example files are given for the London KCL site, 2011 data (denoted :code:`Kc11`)

In the following, :code:`SS` is the site code (e.g. :code:`Kc`), :code:`ss` the grid ID, :code:`YYYY` the year and :code:`tt` the time interval.

.. list-table::
   :widths: 33 33 33
   :header-rows: 1

   * - Filename
     - Description
     - Input/output
   * - SSss_data.txt
     - Meteorological input
     - Input file (60-min)
   * - SSss_YYYY_data_5.txt
     - Meteorological input
     - Input file (5-min)
   * - InitialConditionsSSss
     - Initial conditions
     - Input - _YYYY.nml(+) file
   * - SUEWS_SiteInfo_SSss.x
     - Spreadsheet
     - Input lsm containing all other input information
   * - RunControl.nml
     - Sets model run
     - Input (located in options main directory)
   * - SS_Filechoices.txt
     - Summary of model run
     - Output  options
   * - SSss_YYYY_5.txt
     - (Optional) 5-min
     - Output resolution output file
   * - SSss_YYYY_60.txt
     - 60-min resolution
     - Output output file
   * - SSss_DailyState.txt
     - Daily state variables
     - Output (all years in one file)





(+) There is a second file InitialConditionsSSss_YYYY_EndOfRun.nml or
InitialConditionsSSss_YYYY+1.nml in the input directory. At the end of
the run, and at the end of each year of the run, these files are written
out so that this information could be used to initialize further model
runs.

Run the model for example data
------------------------------

Before running the model with your own data, check that you get the same results as the test run
example files provided. Copy the example output files elsewhere so you can compare the
results. When you run the program it will write over the supplied files.

To run the model you can use **Command Prompt** (in the directory where
the programme is located type the model name) or just double click the
executable file.

Please see `Troubleshooting` if you have problems running the model.

Preparation of data
-------------------

The information required to run SUEWS for your site consists of:

#. Continuous *meteorological forcing data* for the entire period to be modelled without gaps. If you need help preparing the data you can use some of the `UMEP`_ tools.
#. Knowledge of the *surface and soil conditions immediately prior to the first model timestep*. If these initial conditions are unknown, model spinup can help; i.e. run the model and use the output at the end of the run to infer the conditions at the start of the main run).
#. The *location of the site* (latitude, longitude, altitude).
#. Information about the *characteristics of the surface*, including
   land cover, heights of buildings and trees, radiative characteristics
   (e.g. albedo, emissivity), drainage characteristics, soil
   characteristics, snow characteristics, phenological characteristics
   (e.g. seasonal cycle of LAI). For guidance on how to derive parameters related to
   LAI, albedo, surface conductance and surface roughness, the reader is referred to this `link <https://suews-parameters-docs.readthedocs.io/>`_.
#. Information about *human behaviour*, including energy use and water
   use (e.g. for irrigation or street cleaning) and snow clearing (if
   applicable). The anthropogenic energy use and water use may be
   provided as a time series in the meteorological forcing file if these
   data are available or modelled based on parameters provided to the
   model, including population density, hourly and weekly profiles of
   energy and water use, information about the proportion of properties
   using irrigation and the type of irrigation (automatic or manual).

It is particularly important to ensure the following input information
is appropriate and representative of the site:

-  Fractions of different land cover types and (less so) heights of
   buildings [W16]_
-  Accurate meteorological forcing data, particularly precipitation and
   incoming shortwave radiation [Ko17]_
-  Initial soil moisture conditions [Best2014]_
-  Anthropogenic heat flux parameters, particularly if there are
   considerable energy emissions from transport, buildings, metabolism,
   etc [W16]_
-  External water use (if irrigation or street cleaning occurs)
-  Snow clearing (if running the snow option)
-  Surface conductance parameterisation [J11]_ [W16]_

SUEWS can be run either for an individual area or for multiple areas.
There is no requirement for the areas to be of any particular shape but
here we refer to them as model 'grids'.

Preparation of site characteristics and model parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The area to be modelled is described by a set of characteristics that
are specified in the `SUEWS_SiteSelect.txt`
file. Each row corresponds to one model grid for one year (i.e. running
a single grid over three years would require three rows; running two
grids over two years would require four rows). Characteristics are often
selected by a code for a particular set of conditions. For example, a
specific soil type (links to `SUEWS_Soil.txt`) or
characteristics of deciduous trees in a particular region (links to
`SUEWS_Veg.txt`). The intent is to build a library of
characteristics for different types of urban areas. The codes are
specified by the user, must be integer values and must be unique within
the first column of each input file, otherwise the model will return an
error. (Note in `SUEWS_SiteSelect.txt` the first column is labelled 'Grid' and can contain repeat values for different years.) See `Input_files` for details. Note `UMEP`_ maybe helpful for components of this.

Land cover
^^^^^^^^^^

For each grid, the land cover must be classified using the following
surface types:

.. list-table::
   :widths: 25 25 50
   :header-rows: 1
   :stub-columns: 1

   * - Classification
     - Surface type
     - File where characteristics are specified
   * - Non-vegetated
     - Paved surfaces
     - `SUEWS_NonVeg.txt`
   * -
     - Building
     - `SUEWS_NonVeg.txt`
   * -
     - Bare soil
     - `SUEWS_NonVeg.txt`
   * - Vegetation
     - Evergreen trees
     - `SUEWS_Veg.txt`
   * -
     - Deciduous trees
     - `SUEWS_Veg.txt`
   * -
     - Grass
     - `SUEWS_Veg.txt`
   * - Water
     - Water
     - `SUEWS_Water.txt`
   * - Snow
     - Snow
     - `SUEWS_Snow.txt`


The surface cover fractions (i.e. proportion of the grid taken up by
each surface) must be specified in
`SUEWS_SiteSelect.txt`. The surface cover
fractions are **critical**, so make certain that the different surface
cover fractions are appropriate for your site.

For some locations, land cover information may be already available
(e.g. from various remote sensing resources). If not, websites like Bing
Maps and Google Maps allow you to see aerial images of your site and can
be used to estimate the relative proportion of each land cover type. If
detailed spatial datasets are available,
`UMEP`_ allows for a direct link
to a GIS environment using QGIS.

.. _anthropogenic-heat-flux-qf-1:

Anthropogenic heat flux (|QF|)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can either model |QF| within SUEWS or provide it as an input.

-  To model it population density is needed as an input for LUMPS and
   SUEWS to calculate |QF|.
-  If you have no information about the population of the site we
   recommend that you use the `LUCY`_ model [lucy]_  [lucy2]_ to estimate the
   anthropogenic heat flux which can then be provided as input SUEWS
   along with the meteorological forcing data.

Alternatively, you can use the updated version of LUCY called
`LQF`_, which is included in
`UMEP`_.

Other information
^^^^^^^^^^^^^^^^^

The surface cover fractions and population density can have a major
impact on the model output. However, it is important to consider the
suitability of all parameters for your site. Using inappropriate
parameters may result in the model returning an error or, worse,
generating output that is simply not representative of your site. Please
read the section on `input_files`. Recommended or
reasonable ranges of values are suggested for some parameters, along
with important considerations for how to select appropriate values for
your site.

.. _data_entry:

Data Entry
^^^^^^^^^^

To create the series of input text files describing the characteristics
of your site, there are three options:

#. Data can be entered directly into the input text files. The example
   (.txt) files provide a template to create your own files which can be
   edited with :ref:`A_text_editor` directly.
#. Data can be entered into the spreadsheet **SUEWS_SiteInfo.xlsm** and
   the input text files generated by running the macro.
#. Use `UMEP`_.

**To run the xlsm macro:** Enter the data for your site into the xlsm
spreadsheet **SUEWS_SiteInfo.xlsm** and then use the macro to create the
text files which will appear the same directory.

If there is a problem

-  Make sure none of the text files to be generated are open.
-  It is recommended to close the spreadsheet before running the actual
   model code.

Note that in all txt files:

-  The first two rows are headers. The first row is the column number;
   the second row is the column name.
-  The names and order of the columns should not be altered from the
   templates, as these are checked by the model and errors will be
   returned if particular columns cannot be found.
-  Since v2017a it is no longer necessary for the meteorological forcing
   data to have two rows with -9 in column 1 as their last two rows.
-  “!” indicates a comment, so any text following "!" on the same line
   will not be read by the model.
-  If data are unavailable or not required, enter the value -999 in the
   correct place in the input file.
-  Ensure the units are correct for all input information. See `Input_files` for a description of parameters.

In addition to these text files, the following files are also needed to
run the model.

Preparation of the RunControl file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the RunControl.nml file the site name (:code:`SS`) and directories for the
model input and output are given. This means **before running** the
model (even the with the example datasets) you must either

#. open the RunControl.nml file and edit the input and output file paths
   and the site name (with a :ref:`A_text_editor`) so that
   they are correct for your setup, or
#. create the directories specified in the RunControl.nml file

From the given site identification the model identifies the input files
and generates the output files. For example if you specify::

    FileOutputPath = “C:\FolderName\SUEWSOutput\” 

and use site code SS the model creates an output file::

    C:\FolderName\SUEWSOutput\SSss_YYYY_TT.txt 

.. note:: remember to add the last backslash in windows and slash in Linux/Mac


If the file paths are not correct the program will return an error when
run and write the error to the `problems.txt` file.

Preparation of the Meteorological forcing data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The model time-step is specified in `RunControl.nml`
(5 min is highly recommended). If meteorological forcing data are not
available at this resolution, SUEWS has the option to downscale (e.g.
hourly) data to the time-step required. See details about the
`SSss_YYYY_data_tt.txt` to learn more
about choices of data input. Each grid can have its own meteorological
forcing file, or a single file can be used for all grids. The forcing
data should be representative of the local-scale, i.e. collected (or
derived) above the height of the roughness elements (buildings and
trees).

Preparation of the InitialConditions file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Information about the surface state and meteorological conditions just
before the start of the run are provided in the Initial Conditions file.
At the very start of the run, each grid can have its own Initial
Conditions file, or a single file can be used for all grids. For details
see `Initial_Conditions`.

Run the model for your site
---------------------------

To run the model you can use **Command Prompt** (in the directory where
the programme is located type the model name) or just double click the
executable file.

Please see `Troubleshooting` if you have problems
running the model.

Analyse the output
------------------

It is a good idea to perform initial checks that the model output looks
reasonable.

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Characteristic
     - Things to check
   * - Leaf area index
     - Does the phenology look appropriate?
        * what does the seasonal cycle of `leaf area index (LAI) <http://glossary.ametsoc.org/wiki/Leaf_area_index>`__ look like?
        * Are the leaves on the trees at approximately the right time of the year?
   * - Kdown
     - Is the timing of diurnal cycles correct for the incoming solar radiation?
        * Although Kdown is a required input, it is also included in the output file. It is a good idea to check that the timing of Kdown in the output file is appropriate, as problems can indicate errors with the timestamp, incorrect time settings or problems with the disaggregation. In particular, make sure the sign of the longitude is specified correctly in `SUEWS_SiteSelect.txt`.
        * Checking solar angles (zenith and azimuth) can also be a useful check that the timing is correct.
   * - Albedo
     -
      Is the bulk albedo correct?
        * This is critical because a small error has an impact on all the fluxes (energy and hydrology).
        * If you have measurements of outgoing shortwave radiation compare these with the modelled values.
        * How do the values compare to literature values for your area?



Summary of files
----------------

The table below lists the files required to run SUEWS and the output
files produced. SS is the two-letter code (specified in RunControl)
representing the site name, ss is the grid identification (integer
values between 0 and 2,147,483,647 (largest 4-byte integer)) and YYYY is
the year. TT is the resolution of the input/output file and tt is the
model time-step.

The last column indicates whether the files are needed/produced once per
run (1/run), or once per day (1/day), for each year (1/year) or for each
grid (1/grid)::

    [B] indicates files used with the CBL part of SUEWS (BLUEWS) and therefore are only needed/produced if this option is selected
    [E] indicates files associated with ESTM storage heat flux models and therefore are only needed/produced if this option is selected

Get in contact
--------------
For issues met in using SUEWS, we recommend the following ways to get in contact with the developers and the SUEWS community:

#. Report issues on `our GitHub page <https://github.com/Urban-Meteorology-Reading/Urban-Meteorology-Reading.github.io/issues>`_.

#. Ask for help by joining `the Email-list for SUEWS <https://www.lists.reading.ac.uk/mailman/listinfo/met-suews>`_.


.. _`tutorials`: http://umep-docs.readthedocs.io/en/latest/Tutorials/Tutorials.html
.. _`UMEP`: http://umep-docs.readthedocs.io/en/latest/index.html
.. _`LQF`: http://umep-docs.readthedocs.io/en/latest/OtherManuals/LQF_Manual.html
.. _`LUCY`: https://urban-meteorology-reading.github.io/LUCY
