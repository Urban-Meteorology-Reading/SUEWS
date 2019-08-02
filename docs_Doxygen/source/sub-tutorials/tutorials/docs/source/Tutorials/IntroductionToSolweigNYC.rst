.. _IntroductionToSOLWEIGNYC:

Thermal Comfort - Introduction to SOLWEIG (NYC)
===============================================

Introduction
------------

In this tutorial you will use a model **SOlar and LongWave Environmental
Irradiance Geometry model (SOLWEIG)** to estimate the mean radiant
temperature (T\ :sub:`mrt`).

SOLWEIG is a model that simulates spatial variations of 3D radiation
fluxes and the T\ :sub:`mrt` in complex urban settings. It is also able
to model spatial variations of shadow patterns. T\ :sub:`mrt` is one of
the key meteorological variables governing human energy balance and the
thermal comfort of people. It is derived from summing all the radiative
(shortwave and longwave) fluxes (both direct and reflected) to which the
human body is exposed. In SOLWEIG, T\ :sub:`mrt` is derived by modelling
shortwave and longwave radiation fluxes in six directions (upward,
downward and from the four cardinal points) and angular factors.

The model requires **meteorological** forcing data (global shortwave
radiation (K\ :sub:`down`), air temperature (T\ :sub:`a`), relative humidity (RH)),
urban geometry (DSMs), and geographic information (latitude, longitude
and elevation). To determine T\ :sub:`mrt`, continuous maps of sky view
factors are required. Both vegetation and ground cover information can
be added to increase the accuracy of the model output. Below 
a schematic flowchart of SOLWEIG in shown. The `full
manual <http://umep-docs.readthedocs.io>`__ provides more
detail.

.. figure:: /images/SOLWEIG_flowchart.png
   :alt:  Overview of SOLWEIG

   Overview of SOLWEIG

Objectives
----------

To introduce SOLWEIG and how to run the model within `UMEP (Urban
Multi-scale Environmental Predictor) <http://umep-docs.readthedocs.io>`__. 

Help with Abbreviations can be found `here <http://umep-docs.readthedocs.io/en/latest/Abbreviations.html>`__.

Steps
~~~~~

#. Generation of the different kinds of input data that are needed to
   run the model
#. How to run the model
#. How to examine the model output
#. Add additional information (vegetation and ground cover) to improve
   the model outcome and to examine the effect of climate sensitive
   design

Initial Practical steps
-----------------------

UMEP is a python plugin used in conjunction with
`QGIS <http://www.qgis.org>`__. To install the software and the UMEP
plugin see the `getting
started <http://umep-docs.readthedocs.io/en/latest/Getting_Started.html>`__
section in the UMEP manual.

As UMEP is under constant development, some documentation may be missing
and/or there may be instability. Please report any issues or suggestions
to our `repository <https://bitbucket.org/fredrik_ucg/umep/>`__.

Data for this exercise
~~~~~~~~~~~~~~~~~~~~~~

The UMEP tutorial datasets around CCNY campus can be downloaded from our here repository
`here <https://github.com/Urban-Meteorology-Reading/Urban-Meteorology-Reading.github.io/raw/master/other%20files/CCNY_ESPG26918.zip>`__

-  Download, extract and add the raster layers (DSM, CDSM, DEM and land
   cover) into a new QGIS session (see below).

   -  Create a new project
   -  Examine the geodata by adding the layers (*ccny_dsm_1m*,
      *ccny_cdsm_1m*, *ccny_dem_1m* and *ccny_lc_1m*) to your project (***Layer
      > Add Layer > Add Raster Layer**).

-  Coordinate system of the grids is NAD83 / UTM zone 18N (EPSG:26918). If you
   look at the lower right hand side you can see the CRS used in the
   current QGIS project
-  Examine the different datasets before you move on.

-  To add a legend to the **land cover** raster you can load
   **landcoverstyle.qml** found in the test dataset. Right click on the
   land cover (*Properties -> Style (lower left) -> Load Style*).

The origin of the datasets are shown below
   
.. list-table:: Table 1. Spatial data used for this tuorial
   :widths: 10 10 40 40

   * - **Geodata**
     - **Year**
     - **Source**
     - **Description**
   * - Digital surface model (DSM)
     - 2013 (Lidar), 2016 (building polygons)
     - United States Geological Survey (USGS). New York CMGP Sandy 0.7m NPS Lidar and NYC Open Data Portal. `link <https://data.cityofnewyork.us>`__
     - A raster grid including both buildings and ground given in meter above sea level.
   * - Digital elevation model (DEM)
     - 2013
     - United States Geological Survey (USGS). New York CMGP Sandy 0.7m NPS Lidar. `link <https://data.cityofnewyork.us>`__	
     - A raster grid including only ground heights given in meter above sea level.
   * - Digital canopy model (CDSM)
     - 2013 (August)
     - United States Geological Survey (USGS). New York CMGP Sandy 0.7m NPS Lidar. `link <https://coast.noaa.gov/htdata/lidar1_z/geoid12b/data/4920/>`__
     - A vegetation raster grid where vegetation heights is given in meter above ground level. Vegetation lower than 2.5 meter Pixels with no vegetation should be zero.
   * - Land cover (UMEP formatted)
     - 2010
     - New York City Landcover 2010 (3ft version). University of Vermont Spatial Analysis Laboratory and New York City Urban Field Station. `link <https://opendata.cityofnewyork.us/>`__
     - A raster grid including: 1. Paved surfaces, 2. Building surfaces, 3. Evergreen trees and shrubs, 4. Deciduous trees and shrubs, 5. Grass surfaces, 6. Bare soil, 7. Open water			

     
   
SOLWEIG Model Inputs
--------------------

Details of the model inputs and outputs are provided in the `SOLWEIG
manual <http://umep-docs.readthedocs.io/en/latest/OtherManuals/SOLWEIG.html>`__. As this tutorial is
concerned with a **simple application** only the most critical
parameters are used. Many other parameters can be modified to more
appropriate values if applicable. The table below provides an overview
of the parameters that can be modified in the Simple application of
SOLWEIG.

Data requreiments:
R: required, O: Optional, N : not needed, 
S: Spatial, M: Meteorological, 

.. list-table:: Input data and parameters
   :widths: 30 30 5 5 30

   * - **Data**
     - **Definition**
     - **Use**
     - **Type**
     - **Description**
   * - Ground and building digital surface model (DSM)
     - High resolution surface model of ground and building heights
     - R
     - S
     - Given in metres above sea level (m asl)
   * - Digital elevation model (DEM) 
     - High resolution surface model of the ground 
     - R\* 
     - S 
     - R\* if land cover is absent to identify buildings. Given in m asl. Must be same resolution as the DSM.
   * - Digital canopy surface model (CDSM) 
     - High resolution surface model of 3D vegetation 
     - O 
     - S
     - Given in metres above ground level (m agl). Must be same resolution as the DSM.
   * - Digital trunk zone surface model (TDSM) 
     - High resolution surface model of trunk zone heights (underneath tree canopy) 
     - O 
     - S 
     - Given in m agl. Must be same resolution as the DSM.
   * - Land (ground) cover information (LC) 
     - High resolution surface model of ground cover 
     - O 
     - S 
     - Must be same resolution as the DSM. Five different ground covers are currently available (building, paved, grass, bare soil and water)
   * - UMEP formatted meteorological data 
     - Meteorological data from one nearby observation station, preferably at 1-2 m above ground. 
     - R 
     - M 
     - Any time resolution can be given.
   * - Latitude (°) 
     - Solar related calculations 
     - R 
     - O
     - Obtained from the ground and building DSM coordinate system
   * - Longitude (°) 
     - Solar related calculations 
     - R
     - O
     - Obtained from the ground and building DSM coordinate system
   * - `UTC (h) <https://en.wikipedia.org/wiki/Coordinated_Universal_Time>`__
     - Time zone 
     - R
     - O 
     - Influences solar related calculations. Set in the interface of the model.
   * - Human exposure parameters 
     - Absorption of radiation and posture 
     - R 
     - O 
     - Set in the interface of the model.
   * - Environmental parameters
     - e.g. albedos and emissivites of surrounding urban fabrics 
     - R 
     - O 
     - Set in the interface of the model.
 

Meterological input data should be in UMEP format. You can use the
`Meterological Preprocessor <http://umep-docs.readthedocs.io/en/latest/pre-processor/Meteorological%20Data%20MetPreprocessor.html>`__
to prepare your input data. There is also a possibility to use a single point in time in the plugin. 

Requred meteorological data is: 

#. Air temperature (°C)
#. Relative humidity (%)
#. Incoming shortwave radiation (W m\ :sup:`2`)

The model performance will increase if also diffure and direct beam solar radiation is 
available but the mdoel can also calculate these variables. 


How to Run SOLWEIG from the UMEP-plugin
---------------------------------------

#. Open SOLWEIG from *UMEP -> Processor -> Outdoor Thermal Comfort ->
   Mean radiant temperature (SOLWEIG)*.

   -  Some additional information about the plugin is found in the lower
      left window. You will make use of a test dataset from observations
      for Gothenburg, Sweden.

    .. figure:: /images/SOLWEIG_Interface.png
       :alt:  None
       :width: 100%

       Dialog for the SOLWEIG model (click on image for larger image)

#. To be able to run the model some additional spatial datasets needs to
   be created.

   -  Close the SOLWEIG plugin and open *UMEP -> Pre-Processor -> Urban
      geometry -> Sky View Factor*.
   -  To run SOLWEIG various sky view factor (SVF) maps for both
      vegetation and buildings must be created (see `Lindberg and
      Grimmond
      (2011) <http://link.springer.com/article/10.1007/s00704-010-0382-8>`__
      for details).
   -  You can create all SVFs needed (vegetation and buildings) at the
      same time. Use the settings as shown below. Use an appropriate
      output folder for your computer. 
	  
    .. figure:: /images/SOLWEIG_Svf_solweig.png
       :alt:  None
       :width: 487px
       
       Settings for the SkyViewFactorCalculator.
	   
   -  When the calculation is done, map will appear in the map canvas.
      This is the 'total' SVF i.e., including both buildings and
      vegetation. Examine the dataset.
   -  Where are the highest and lowest values found?
   -  Look in your output folder and find a zip-file containing all the
      necessary SVF maps needed to run the SOLWEIG-model.

#. Another preprocessing plugin needed is to create the building wall
   heights and aspect. Open *UMEP -> Pre-Processor -> Urban geometry ->
   Wall height and aspect* and use the settings as shown below. QGIS scales loaded rasters by a *cumulative count out* approach (98%). As the height and aspect layers are filled with zeros where no wall are present it might appear as there is no walls identified. Rescale your results to see the wall identified (*Layer Properties > Style*).
   
    .. figure:: /images/SOLWEIG_wallgeight_solweig.png
       :alt:  None
       :width: 505px
       
       Settings for the Wall height and aspect plugin.

#. Re-open the SOLWEIG plugin and use the settings as shown below. **NOTE** Remember to change UTC time to *-5* (figure not updated). 
   You will use the GUI to set one point in time (i.e. a summer hour in
   Gothenburg, Sweden) hence, no input meteorological file is needed for
   now. No information on vegetation and ground cover is added for this
   first try. Click **Run**. 
   
    .. figure:: /images/SOLWEIG_Tmrt1_solweig.png
       :alt:  None
       :width: 100%
       
       The settings for your first SOLWEIG run (click on image for larger image).
	   
#. Examine the output (Average T\ :sub:`mrt` (°C). What is the main
   driver to the spatial variations in T\ :sub:`mrt`?
#. Add 3D vegetation information by ticking in *Use vegetation scheme
   (Lindberg, Grimmond 2011)* and add the vegetation canopy dataset as the *Vegetation
   Canopy DSM*. As no TDSM exists we estimate the it by using 25% of the
   canopy height. Leave the tranmissivity as 3%. Tick in *Save generated
   Trunk Zone DSM* (a tif file, **TDSM.tif**, will be generated in the
   specified output folder and used in a later section: **Climate
   sensitive planning**). Also tick in *Save generated building grid* as
   this will be needed later in this tutorial. Leave the other setting
   as before (Step 4) except for changing your output directory
   Otherwise, results from your first run will be overwritten. Run the
   model again and compare the result with your first run.
#. Add your last spatial dataset, the **land cover** grid by ticking in
   *Use land cover scheme (Lindberg et al. 2016)*. Run and compare the
   result again with the previous runs.

Using meteorolgical data and POIs
---------------------------------

SOLWEIG is also able to run a continuous dataset of meteorological data.
You will make use of a single summer day as well as a winter day for
Gothenburg, Sweden. The GUI is also able to derive full model output
(all calculated variables) from certain points of interest (POIs).

#. First you need to create a point vector layer to store the POIs. Go
   to *Layer -> Create Layer -> New Shape file*. Choose *Point* as
   *Type* and add a new text field called **name**. Name the new layer
   **POI_Kr.shp**. Specify the coordinate system to be same as the other dataset.
#. Now you should add two points within the study area. To add points to
   the layer it has to be editable and Add Feature should be activated.

    .. figure:: /images/SOLWEIG_Addpoint.png
       :alt:  None
       :width: 411px
       
       Setting to add points 
   
   Two points should be added and the attributes should be id=\ **1** and
   name=\ **point1** for the right point and id=\ **2** and
   name=\ **point2** for the left point. Put them out at any location within the domain
   
   When you are finished, save layer edits. Close the editing by pressing Toggle editing (the pencil).
#. Now open the SOLWEIG plugin. Use both the vegetation and land cover
   schemes as before. This time, tick in *Include POI(s)*, select your
   point layer and use the ID attribute as *ID field*.
#. Tick in *Use continuous meteorological dataset* and choose
   **NYC_2010_Data_60_doy190.txt** as *Input meteorological file*. This is a warm and clear day in 2010 (9 July).
   Also, tick in to save T\ :sub:`mrt` as *Output maps*. Run the model again.

Examine your output with SOLWEIG Analyzer
-----------------------------------------

To perform a first set of analysis of your result you can make use of
the SOLWEIG Analyzer plug-in.

#. Open the Analyzer located in *UMEP -> Post-Processor -> Outdoor
   Thermal Comfort -> SOLWEIG Analyzer*. Here you can analyze both data
   from your POIs as well as perform statistical analysis based on saved
   output maps. Start by locating your output folder in the top section
   (*Load Model Result*). 
   
    .. figure:: /images/SOLWEIG_SOLWEIGAnalyzer.png
       :alt:  None
       :width: 100%
       
       Dialog for the SOLWEIG Analyzer plug-in

#. Firstly you will compare differences in T\ :sub:`mrt` for the two
   locations (point1 and point2). This can done using the left frame
   (*Point of Interest data*). Specify *point1* as *POI* and *Mean
   Radiant Temperature* in the two top scroll down lists. Then tick in
   *Include another POI/variable* and chose *point2* and *Mean Radiant
   Temperature* below. Click *Plot*. What explains the differences?
#. Now lets us move on to analyse the output maps generated from our
   last model run. In the right frame, specify *Mean Radiant
   Temperature* as *Variable to visualize*. Start by clicking *Show
   Animation*. Now the output maps of T\ :sub:`mrt` generated before are
   displayed in a sequence.
#. Next step is to generate some statistical maps from the last model
   run. Specify *Mean Radiant Temperature* as *Variable to visualize*
   and tick in to *Exclude building pixels*. Choose the building grid
   that you saved earlier in this tutorial. If it is not in the
   drop-down list you need to add this layer (**buildings**) to your
   project. Tick in *T*\ :sub:`mrt`: *Percent of time above threshold
   (degC)* and specify 55.0 as threshold. Specify an output folder and
   tick also in *Add analysis to map canvas* before you generate the
   result. The resulting map show the time that a pixel has been above
   55 degC based on the whole analysis time i.e. 24 hours. This type of
   maps can be used to identify areas prone to e.g. heat stress

Climate sensitive planning
--------------------------

Vegetation is one effective measure to reduce areas prone to heat
related health issues. In this section you make use of the Tree
Generator plugin to see the effect of adding more vegetation into our
study area. The municipality in Gothenburg have identified a "hot spot"
south of the german church and they want to see the effect of planting
three new trees in that area.

The Tree Generator
~~~~~~~~~~~~~~~~~~

The Tree Generator plugin make use of a point vector file including the
necessary attributes to generate/add/remove vegetation suitable for
either mean radiant temperature modelling with SOLWEIG or urban energy
balance modelling with SUEWS.

#. Create a point vector shape file named (**Trees.shp**) as described
   in the previous section adding five attributes (*id, ttype, trunk,
   totheight, diameter*). The attributes should all be decimal (float)
   numbers (see table below). Locate the new tree on the laws south west of Shepards Hall 
   The values for all three vegetation units should
   be **ttype=2, trunk=4, totheight=15, diameter=10**. 
   
#. Add your created trunk zone dsm (TDSM.tif) that was created
   previously (located in your output directory).
#. Open the TreeGenerator (UMEP -> PreProcessor -> TreeGenerator) and
   use the settings as shown in figure below. 

    .. figure:: /images/SOLWEIG_Treegeneratorsolweig.png
       :alt:  None
       :width: 574px
       
       The settings for the Tree Generator

#. As the vegetation DSMs have been changed, the SVFs has to be
   recalculated. This time use the two generated vegetation DSMs.
#. Now re-run SOLWEIG using the same settings as before but now use the
   new vegetation surface models as well as the new SVFs generated in
   the previous step.
#. Generate a new, updated threshold map based on the new results and
   compare the differences.

The table below show the input variables needed for each tree point.

+-----------------------+-----------------------+-----------------------+
| Attribute name        | Name                  | Description           |
+=======================+=======================+=======================+
| ttype                 | Tree type             | Two shapes are        |
|                       |                       | available:            |
|                       |                       |                       |
|                       |                       | -  conifer = 1 and    |
|                       |                       | -  deciduous = 2.     |
|                       |                       | -  To remove          |
|                       |                       |    vegetation set     |
|                       |                       |    ttype = 0.         |
+-----------------------+-----------------------+-----------------------+
| trunk                 | Trunk zone height (m  | Height of the trunk   |
|                       | agl)                  | zone.                 |
+-----------------------+-----------------------+-----------------------+
| totheight             | Total tree height (m  | Maximum height of the |
|                       | agl)                  | vegetation unit       |
+-----------------------+-----------------------+-----------------------+
| diameter              | Canopy diameter (m)   | Circular diameter of  |
|                       |                       | the vegetation unit   |
+-----------------------+-----------------------+-----------------------+



