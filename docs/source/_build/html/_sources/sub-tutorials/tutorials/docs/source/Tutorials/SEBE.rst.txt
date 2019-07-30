.. _SEBE:

Solar Energy - Introduction to SEBE
===================================

Introduction
------------

Lindberg et al.’s (2015) solar radiation model SEBE (Solar Energy on
Building Envelopes) allows estimates of solar irradiance on ground
surfaces, building roofs and walls. It uses a shadow casting algorithm
with a digital surface model (**DSM**) and the solar position to
generate pixel-level information of shadow or sunlit areas. The shadow
casting algorithm (Ratti and Richens 1999) has been developed to
incorporate walls (Lindberg et al. 2015). This is of interest for a
broad range of applications, for example solar energy potential and
thermal comfort.

SEBE is incorporated in `UMEP(Urban Multi-scale Environmental
Predictor) <http://umep-docs.readthedocs.io>`__, a plugin for
`QGIS <http://www.qgis.org>`__. As SEBE was initially developed to
estimate solar energy potential on building roofs, the Digital Surface
Models (DSMs) used need to include roof structures, such as tilted
roofs, chimneys etc. Methods to produce accurate ground and building
DSMs for SEBE include the use of LiDAR technology and 3D roof structure
objects in vector format.

In this exercise you will apply the model in `Gothenburg,
Sweden <https://en.wikipedia.org/wiki/Gothenburg>`__ and `Covent
Garden <https://en.wikipedia.org/wiki/Covent_Garden>`__, London to
investigate solar energy potential, how it changes between buildings,
with seasons, with the effects of vegetation etc.

Initial Practical steps
-----------------------

-  If **QGIS** is not on your computer you will `need to install it <http://umep-docs.readthedocs.io/en/latest/Getting_Started.html>`__
-  Then install the `**UMEP** plugin <http://umep-docs.readthedocs.io/en/latest/Getting_Started.html>`__

-  Start the QGIS software
-  *Windows:* If not visible on the desktop use the **Start** button to
   find the software (i.e. Find QGIS under your applications)
-  Select **QGIS 2.18 Desktop** (or the latest version installed). Do not use QGIS3 at this point.

When you open it on the top toolbar you will see **UMEP**.

.. figure:: /images/SEBE_Interfacelocation.png
   :alt:  None
   :width: 100%

   Location of SEBE in UMEP

-  Read through the section in the online
   manual BEFORE using the model, so you are familiar with it’s operation and
   terminology used.

Data for Tutorial
~~~~~~~~~~~~~~~~~

.. figure:: /images/SEBE_Gothenburg.png
   :alt:  None
   :width: 356px

   Central Gothenburg study area (red square).
   The Open layers plugin in QGIS was used to generate
   this snapshot.

Geodata and meteorological data for **Gothenburg, Sweden**.

-  Data are projected in SWEREF99 1200 (EPSG:3007) the national
   coordinate system of Sweden.

Data requreiments:
S: Spatial, M: Meteorological,

.. list-table:: Input data and parameters
   :widths: 20 20 10 50

   * - **Name**
     - **Definition**
     - **Type**
     - **Description**
   * - krbig_dsm.asc
     - Ground and building DSM
     - S
     - Raster dataset: derived from a 3D vector roof structure dataset and a digital elevation model (DEM)
   * - krbig_cdsm.asc
     - Vegetation canopy DSM
     - S
     - Raster dataset: drived from a LiDAR dataset
   * - kr_buildings.shp
     - Building footprint polygon layer
     - S
     - Vector dataset
   * - GBG_typicalweatheryear_1977.txt
     - Meteorological forcing data
     - M
     - Meteorological data, hourly time resolution for 1977 Gothenburg, Sweden.


:download:`Datasets in Gothenburg, Sweden </data/Goteborg_SWEREF99_1200.zip>`

:download:`Datasets in London Covent Garden </data/DataCoventGarden.zip>`

`Google map link to Covent Garden <https://www.google.co.uk/maps/@51.5117012,-0.1231273,356m/data=!3m1!1e3>`__


Steps
-----

#. Start with the Gothenbrug data. If the data are zipped - unzip the data first.
#. Examine the geodata by adding the layers to your project.
#. Use *Layer > Add Layer > Add Raster
   Layer* to open the .asc raster files
   and *Layer > Add Layer > Add Vector Layer*. The Vector layer is a
   shape file which consists of multiple files. It is the
   **kr_building.shp** that should be used to load the vector layer into
   QGIS.
#. You will need to indicate the co-ordinate system
   (`CRS <https://docs.qgis.org/2.18/en/docs/gentle_gis_introduction/coordinate_reference_systems.html>`__)
   that is associated with these data. If you look at the lower right
   hand side you can see the CRS used in the current QGIS project.

   -  You can use the filter to find this then.
   -  Select SWEREF99 1200 as CRS and the files will load into the map
      canvas.
   -  Do this for all of the geodata files.

#. Open the **meteorological file** in a text editor or in a spreadsheet
   such as MS excel or LibreOffice (Open office).

   -  Data file is formatted for the UMEP plugin (in general) and the
      SEBE plugin (in particular).
   -  First four columns are *time related*.
   -  Columns of interest are **kdown, kdiff and kdir**. These are
      related to shortwave radiation and give global, diffuse and direct
      radiation, respectively.
   -  The meteorological file should be at least a year long, but
      preferably multi-year.
   -  One option is to use a `**typical meteorological
      year** <https://en.wikipedia.org/wiki/Typical_meteorological_year>`__
      as you will do in this tutorial

Variables included in the **meteorological data file**. No. indicates
the column the file is in. Use indicates if it is **R – required** or
*O- optional* (in this application) or **N- Not used in this
application**. All columns must be present but can be filled with
numbers to indicate they are not in use (e.g. -999).

+------+------+-------------+-----------------+
| No.  | USE  | Column name | Description     |
+======+======+=============+=================+
| 1    | R    | iy          | Year [YYYY]     |
+------+------+-------------+-----------------+
| 2    | R    | id          | Day of year     |
|      |      |             | [DOY]           |
+------+------+-------------+-----------------+
| 3    | R    | it          | Hour [H]        |
+------+------+-------------+-----------------+
| 4    | R    | imin        | Minute [M]      |
+------+------+-------------+-----------------+
| 5    | N    | qn          | Net all-wave    |
|      |      |             | radiation [W    |
|      |      |             | m\ :sup:`-2`]   |
+------+------+-------------+-----------------+
| 6    | N    | qh          | Sensible heat   |
|      |      |             | flux [W         |
|      |      |             | m\ :sup:`-2`]   |
+------+------+-------------+-----------------+
| 7    | N    | qe          | Latent heat     |
|      |      |             | flux [W         |
|      |      |             | m\ :sup:`-2`]   |
+------+------+-------------+-----------------+
| 8    | N    | qs          | Storage heat    |
|      |      |             | flux [W         |
|      |      |             | m\ :sup:`-2`]   |
+------+------+-------------+-----------------+
| 9    | N    | qf          | Anthropogenic   |
|      |      |             | heat flux [W    |
|      |      |             | m\ :sup:`-2`]   |
+------+------+-------------+-----------------+
| 10   | N    | U           | Wind speed [m   |
|      |      |             | s\ :sup:`-1`]   |
+------+------+-------------+-----------------+
| 11   | O    | RH          | Relative        |
|      |      |             | Humidity [%]    |
+------+------+-------------+-----------------+
| 12   | O    | Tair        | Air temperature |
|      |      |             | [°C]            |
+------+------+-------------+-----------------+
| 13   | N    | pres        | Barometric      |
|      |      |             | pressure [kPa]  |
+------+------+-------------+-----------------+
| 14   | N    | rain        | Rainfall [mm]   |
+------+------+-------------+-----------------+
| 15   | R    | kdown       | Incoming        |
|      |      |             | shortwave       |
|      |      |             | radiation [W    |
|      |      |             | m\ :sup:`-2`]   |
|      |      |             | Must be >= 0 W  |
|      |      |             | m\ :sup:`-2`.   |
+------+------+-------------+-----------------+
| 16   | N    | snow        | Snow [mm]       |
+------+------+-------------+-----------------+
| 17   | N    | ldown       | Incoming        |
|      |      |             | longwave        |
|      |      |             | radiation [W    |
|      |      |             | m\ :sup:`-2`]   |
+------+------+-------------+-----------------+
| 18   | N    | fcld        | Cloud fraction  |
|      |      |             | [tenths]        |
+------+------+-------------+-----------------+
| 19   | N    | Wuh         | External water  |
|      |      |             | use[m\ :sup:`3`]|
+------+------+-------------+-----------------+
| 20   | N    | xsmd        | Observed soil   |
|      |      |             | moisture [m3    |
|      |      |             | m\ :sup:`-3` or |
|      |      |             | kg              |
|      |      |             | kg\ :sup:`-1`]  |
+------+------+-------------+-----------------+
| 21   | N    | lai         | Observed leaf   |
|      |      |             | area index [m2  |
|      |      |             | m\ :sup:`-2`]   |
+------+------+-------------+-----------------+
| 22   | O    | kdiff       | Diffuse         |
|      |      |             | radiation [W    |
|      |      |             | m\ :sup:`-2`]   |
+------+------+-------------+-----------------+
| 23   | O    | kdir        | Direct          |
|      |      |             | radiation [W    |
|      |      |             | m\ :sup:`-2`]   |
+------+------+-------------+-----------------+
| 24   | N    | wdir        | Wind direction  |
|      |      |             | [°]             |
+------+------+-------------+-----------------+



Preparing data for SEBE
-----------------------

SEBE plugin: located at *UMEP -> Processor -> Solar Energy -> Solar
Energy on Building Envelopes (SEBE)* in the menu bar.

.. figure:: /images/SEBE_SEBE1.png
   :alt: SEBE1.png
   :width: 514px

   The interface for SEBE in UMEP

#. *Top frame*: for input data for the SEBE calculations.

   -  Critical is the **building and ground**
      `DSM <http://umep-docs.readthedocs.io/en/latest/Abbreviations.html>`__
      for the calculations in SEBE.
   -  Optionally **vegetation** (trees and bushes) can be included as
      they can shadow buildings, walls and roofs reducing the potential
      solar energy production
   -  Two vegetation DSMs are required when the Use vegetation DSMs is
      ticked:

      + One to describe the top of the vegetation (Vegetation Canopy DSM).

      + One to describe the bottom, underneath the canopies (Vegetation Trunk Zone DSM).

      As Trunk Zone DSMs are very rare, an option to create this from the
      canopy DSM is available. You can set the amount of light (shortwave radiation) that is
      transmitted through the vegetation.

#. Two raster datasets, height and wall aspect, are needed to calculate
   irradiance on building walls.

   -  The average albedo (one value is used for all surfaces) can be
      changed.

#. The
   `UTC <https://en.wikipedia.org/wiki/Coordinated_Universal_Time>`__
   offset is needed to accurately estimate the sun position, positive
   numbers for easterly position and negative for westerly. For example,
   Gothenburg is located in CET which is UTC +1.
#. Meteorological file needs to be specified.
#. Wall data are created with the `UMEP plugin - **Wall Height and
   Aspect** <http://umep-docs.readthedocs.io/en/latest/pre-processor/Urban%20Geometry%20Wall%20Height%20and%20Aspect.html>`__:

   -  This uses a 3 by 3 pixels kernel minimum filter where the four
      cardinal points (N, W, S,E) are investigated. The pixels just
      ‘inside’ the buildings are identified and give values to indicate
      they are a building edge. The aspect algorithm originates from a
      linear filtering technique (Goodwin et al. 2009). It identifies
      the linear features plus (a new addition) the aspect of the
      identified line. Other more accurate techniques include using a
      vector building layer and spatially relating this to the wall
      pixels.

#. UMEP -> Pre-Processor -> Urban Geometry -> Wall Height and Aspect.
#. Close the SEBE plugin and open the Wall and Height and Aspect plugin
#. Use your ground and building DSM as input
#. Tick the option to Calculate wall aspect.
#. Create a folder in your Documents folder called e.g. SEBETutorial
#. Use this to save the result.
#. Name your new raster datasets aspect and height, respectively.
#. Tick: Add result to project and click OK.

Running the model
-----------------

Now you have all data ready to run the model.

.. figure:: /images/SEBE_SEBEnoVeg.png
   :alt:  Settings for running SEBE without vegetation.

   Example of settings for running SEBE without vegetation.

#. First run the model *without* including vegetation.

   -  Open the SEBE-plugin again
   -  Make the setting according to the figure to the LHS
   -  Save your results in a subfolder (**NoVeg**) of **SEBETutorial**.
   -  The model takes some time to calculate irradiance on all the
      surfaces.
   -  The result added to your map canvas is the horizontal radiation,
      i.e. irradiance on the ground and roofs.

#. Run the model again but this time also use the vegetation DSM.

   -  Save your result in a subfolder called **Veg**.

Irradiance on building envelopes (alternatively see the tips below – currrently better)
---------------------------------------------------------------------------------------

To determine the irradiance on building walls:

#. Open the SunAnalyser located at *UMEP > Post-Processor > Solar
   Radiation > SEBE (Visualisation)*.

   -  This can be used to visualize the irradiance on both roofs and
      walls.

#. Choose the input folder where you saved your result for one of the
   runs.
#. Mark an area with the tool (Area of Visualisation) on the map canvas
   by click first once
#. Drag to produce an area
#. Click again to finish.
#. Click Visualise. Now you should be able to see the results in 3D.

**3D Visualisation for Mac currently not working properly**

Use the **Profile tool**, which is a plugin for QGIS, to see the range of values along a transect.

#. Plugins > Profile tool > Terrain profile.

   -  Draw a line across the screen on the area of interest. Double
      click and you will see the profile drawn. Make certain you use the
      correct layer (see Tips).

#. If this is not installed you will need to install it from official
   QGIS-plugin reporistory (Plugins > Manage and Install Plugins).

Solar Energy Potential
----------------------

In order to obtain the solar energy potential for a specific building:

#. The actual area of the roof needs to be considered.
#. Determine the area of each pixel (A\ :sub:`P`): e.g. 1 m\ :sup:`2`
#. As some roofs are tilting the area may be larger for some pixels. The
   actual area (*A*\ :sub:`A`) can be computed from:

      *A*\ :sub:`A` = *A*\ :sub:`P` / *cos(S*\ :sub:`i`)

      where the slope (*S*\ :sub:`i`) of the raster pixel should be in radians (1 deg = pi/180 rad).

**To make a slope raster:**

#. *Raster > Terrain analysis > Slope*. If the tool is missing, Go to
   *Manage and Install Plugins* and activate (*Raster Terrain Analysis
   Plugin*)
#. Use the DSM for elevation layer
#. Create the slope z factor =1 - area

.. figure:: /images/SEBE_Slope.png
   :alt: None

   The Slope tool in QGIS

Use the raster menu: *Raster> Raster Calculator*.

#. To determine the area after you have removed the wall area from the
   buildings.
#. Enter the equation indicated.
#. To visualize where to place solar panels the amount of energy
   received needs to be cost effective. As irradiance below 900 kWh is
   considered to be too low for solar energy production (*Per Jonsson
   personal communication Tyréns Consultancy*), pixel cells lower than
   900 can be filtered out (Figure LHS). Transparency – allows you to
   make visible above a threshold of interest.

   -  Right-click on the Energyyearroof-layer and go to **Properties**
      and then **Transparency**.
   -  Add a custom transparency (green cross) where values between 0 and
      900 are set to 100% transparency.

.. figure:: /images/SEBE_RasterCalculator.png
   :alt: None

   The RasterCalculator in QGIS

Irradiance map with values less than 900 kWh filtered out
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To estimate solar potential on building roofs we can use the Zonal
statistics tool:

#. Raster > Zonal statistics.

   -  Use the roof area raster layer (**energyPerm2_slope65_RoofArea**)
      created before and use **kr_building.shp** as the polygon layer to
      calculate as your zone layer. Make sure that you calculate sum
      statistics.

#. On your building layer – Right click Open Attribute Table
#. Or use the identifier to click a building (polygon) of interest to
   see the statistics you have just calculated

   Note that we will not consider the performance of the solar panels.

   .. figure:: /images/SEBE_GOT_Irradiance.png
      :alt: None

      Irradiance map on building roofs in Gothenburg

Covent Garden data set
----------------------

A second GIS data set is available for the Covent Garden area in London

#. Close the Gothenburg data (it may be easiest to completely close QGIS
   and re-open).
#. Download from
   `1 <https://drive.google.com/open?id=0B7D8dqiua0uzWWhwWmU4c1lnTG8>`__
#. Add the Covent Garden data
#. Extract the data to a directory
#. Load the Raster data (DEM, DSM) files (as you did before)
#. Shadows

   -  `UMEP -> Processor -> Solar Radiation -> Daily Shadow
      pattern <http://umep-docs.readthedocs.io/en/latest/processor/Solar%20Radiation%20Daily%20Shadow%20Pattern.html>`__
   -  Allows you to calculate the shadows for a particular time of day
      and `Day of
      Year <http://disc.sci.gsfc.nasa.gov/julian_calendar.html>`__.

Questions for you to explore with UMEP:SEBE
-------------------------------------------

#. Use the Gothenburg dataset consider the impact of vegetation.

   -  What are the main differences between the two model runs with
      respect to ground and roof surfaces?
   -  To what extent are the building roofs affected by vegetation?

#. Consider the differences between London and Gothenburg. You can run
   the model for different times of the year by modifying the
   meteorological data so the file only has the period of interest.
#. For Covent Garden, determine the solar energy potential for a
   specific building within the model domain. Work in groups to consider
   different areas. What would be the impact of having a smaller/larger
   area domain modelled for this building? Identify the possibilities of
   solar energy production for that building.


References
----------

-  Goodwin NR, Coops NC, Tooke TR, Christen A, Voogt JA 2009:
   Characterizing urban surface cover and structure with airborne lidar
   technology. `Can J Remote Sens
   35:297–309 <http://pubs.casi.ca/doi/abs/10.5589/m09-015?journalCode=cjrs>`__
-  Lindberg F, Jonsson P, Honjo T, Wästberg D 2015: Solar energy on
   building envelopes - 3D modelling in a 2D environment. `Solar Energy.
   115,
   369–378 <http://www.sciencedirect.com/science/article/pii/S0038092X15001164>`__
-  Ratti CF, Richens P 1999: Urban texture analysis with image
   processing techniques Proc CAADFutures99, Atlanta, GA

**Authors of this document**: Lindberg and Grimmond (2015, 2016)

*Contributors to the material covered*

-  University of Gothenburg: Fredrik Lindberg
-  University of Reading: Sue Grimmond
-  Background work also comes from: UK (Ratti & Richens 1999), Sweden
   (Lindberg et al. 2015), Canada (Goodwin et al. 2009)

In the `repository <https://bitbucket.org/fredrik_ucg/umep/>`__ of UMEP you can find the code and report bugs and other suggestions on future improvments.

Tips
----

**Meteorological** file in UMEP has a special format. If you have data
in another format there is a `UMEP plugin that can convert your
meteorological data into the UMEP
format <http://umep-docs.readthedocs.io/en/latest/pre-processor/Meteorological%20Data%20MetPreprocessor.html>`__.

-  Plugin is found at *UMEP > Pre-Processor > Meteorological data >Prepare Existing data*.

Plugin to **visualize data** in 3D: called
`Qgis2Threejs <https://media.readthedocs.org/pdf/qgis2threejs/docs-release/qgis2threejs.pdf>`__.

-  Available for download from the official repository Plugins -> Manage
   and Install Plugins.

.. figure:: /images/SEBE_CoventGarden.png
   :alt: None

   3D visualisation with Qgis2Threejs over Convent Garden

TIFF (TIF) and ASC are **raster data file formats** In the left Hand
Side there is a list of layers.

-  The layer that is checked at the top of the list is the layer that is
   seen, If you want to see another layer you can either:

#*Un-tick the layers above the one you are interested in and/or

#*Move the layer you are interested in to the top of the list by
dragging it.

You can save all of you work for different areas as a project – so you
can return to it as whole.

-  Project > Save as

You can change the *shading etc*. on different layers.

-  Right Click on the Layer name Properties > Style > Singlebandpseudo
   color
-  Choose the color band you would like.
-  Classify
-  Numerous things can be modified from this point.

`UMEP repository <https://bitbucket.org/fredrik_ucg/umep/>`__.
