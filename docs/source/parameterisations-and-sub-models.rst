.. _physics_schemes:

Parameterisations and sub-models within SUEWS
=============================================

Net all-wave radiation, Q\*
---------------------------

There are several options for modelling or using observed radiation
components depending on the data available. As a minimum, SUEWS requires
incoming shortwave radiation to be provided.

#. Observed net all-wave radiation can be provided as input instead of
   being calculated by the model.
#. Observed incoming shortwave and incoming longwave components can be
   provided as input, instead of incoming longwave being calculated by
   the model.
#. Other data can be provided as input, such as cloud fraction (see
   options in `RunControl.nml`).
#. **NARP** (Net All-wave Radiation Parameterization, Offerle et al.
   2003 [O2003]_ , Loridan et al. 2011 [L2011]_ ) scheme calculates outgoing
   shortwave and incoming and outgoing longwave radiation components
   based on incoming shortwave radiation, temperature, relative humidity
   and surface characteristics (albedo, emissivity).



Anthropogenic heat flux, Q\ :sub:`F`
------------------------------------

#. Two simple anthropogenic heat flux sub-models exist within SUEWS:

   -  Järvi et al. (2011) [J11]_ approach, based on heating and cooling
      degree days and population density (allows distinction between
      weekdays and weekends).
   -  Loridan et al. (2011) [L2011]_ approach, based on a linear piece-wise
      relation with air temperature.

#. Pre-calculated values can be supplied with the meteorological forcing
   data, either derived from knowledge of the study site, or obtained
   from other models, for example:

   -  **LUCY** (Allen et al. 2011 [lucy]_, Lindberg et al. 2013 [lucy2]_). A
      new version has been now included in UMEP. To distinguish it is
      referred to as
      `LQF`_
   -  **GreaterQF** (Iamarino et al. 2011 [I11]_). A new version has been
      now included in UMEP. To distinguish it is referred to as
      `GQF`_

Storage heat flux, ΔQ\ :sub:`S`
-------------------------------

#. Three sub-models are available to estimate the storage heat flux:

   -  **OHM** (Objective Hysteresis Model, Grimmond et al. 1991 [G91OHM]_,
      Grimmond & Oke 1999a [GO99QS]_, 2002 [GO2002]_). Storage heat heat flux is
      calculated using empirically-fitted relations with net all-wave
      radiation and the rate of change in net all-wave radiation.
   -  **AnOHM** (Analytical Objective Hysteresis Model, Sun et al.
      2017 [AnOHM17]_). OHM approach using analytically-derived coefficients.
      |NotRecmd|
   -  **ESTM** (Element Surface Temperature Method, Offerle et al.
      2005 [OGF2005]_). Heat transfer through urban facets (roof, wall, road,
      interior) is calculated from surface temperature measurements and
      knowledge of material properties. |NotRecmd|

#. Alternatively, 'observed' storage heat flux can be supplied with the
   meteorological forcing data.

Turbulent heat fluxes, Q\ :sub:`H` and Q\ :sub:`E`
--------------------------------------------------

#. **LUMPS** (Local-scale Urban Meteorological Parameterization Scheme,
   Grimmond & Oke 2002 [GO2002]_) provides a simple means of estimating
   sensible and latent heat fluxes based on the proportion of vegetation
   in the study area.
#. **SUEWS** adopts a more biophysical approach to calculate the latent
   heat flux; the sensible heat flux is then calculated as the residual
   of the energy balance. The initial estimate of stability is based on
   the LUMPS calculations of sensible and latent heat flux. Future
   versions will have alternative sensible heat and storage heat flux
   options.

Sensible and latent heat fluxes from both LUMPS and SUEWS are provided in the `output_files`.
Whether the turbulent heat fluxes are calculated using LUMPS or SUEWS can have a major impact on the results.
For SUEWS, an appropriate surface conductance parameterisation is also critical [J11]_ [W16]_.
For more details see `Differences_between_SUEWS_LUMPS_and_FRAISE` .

Water balance
-------------

The running water balance at each time step is based on the urban water
balance model of Grimmond et al. (1986) [G86]_ and urban
evaporation-interception scheme of Grimmond and Oke (1991) [G91]_.

-  Precipitation is a required variable in the meteorological forcing
   file.
-  Irrigation can be modelled [J11]_ or observed values can be provided
   if data are available.
-  Drainage equations and coefficients to use must be specified in the
   input files.
-  Soil moisture can be calculated by the model.
-  Runoff is permitted:

   -  between surface types within each model grid
   -  between model grids (|NotAvail|)
   -  to deep soil
   -  to pipes.

Snowmelt
--------

The snowmelt model is described in Järvi et al. (2014) [Leena2014]_.
Changes since v2016a:
1) previously all surface states could freeze in 1-h time step, now the freezing surface state is
calculated similarly as melt water and can freeze within the snow pack.
2) Snowmelt-related coefficients have also slightly changed (see
`SUEWS_Snow.txt`).

Convective boundary layer
-------------------------

A convective boundary layer (CBL) slab model (Cleugh and Grimmond
2001 [CG2001]_) calculates the CBL height, temperature and humidity during
daytime (Onomura et al. 2015 [Shiho2015]_).

.. SOLWEIG is fully removed since 2019a

.. Thermal comfort
.. ---------------

.. **SOLWEIG** (Solar and longwave environmental irradiance geometry model,
.. Lindberg et al. 2008 [FL2008]_, Lindberg and Grimmond 2011 [FL2011]_) is a 2D
.. radiation model to estimate mean radiant temperature.

.. .. figure:: /assets/img/Bluews_2.jpg
..     :alt:  Overview of scales. Source: Onomura et al. (2015) [Shiho2015]_

..     Overview of scales. Source: Onomura et al. (2015) [Shiho2015]_

Surface Diagnostics
-------------------

A `MOST <https://en.wikipedia.org/wiki/Monin–Obukhov_similarity_theory>`_-based surface diagnostics module is implemented in 2017b for calculating the surface level diagnostics, including:

  * T2: air temperature at 2 m agl
  * Q2: air specific humidity at 2 m agl
  * U10: wind speed at 10 m agl

The details for formulation of these diagnostics can be found in equations 2.54, 2.55 and 2.56 in Brutsaert (2005) [B05]_


.. _LQF: http://umep-docs.readthedocs.io/en/latest/OtherManuals/LQF_Manual.html
.. _GQF: http://umep-docs.readthedocs.io/en/latest/OtherManuals/GQF_Manual.html

.. _rsl_mod:

Wind, Temperature and Humidity Profiles in the Roughness Sublayer
----------------------------------------------------------------------------
Wind, temperature and humidity profiles are derived at 30 levels in the surface layer.
In order to account for the roughness sublayer and canopy layer,
we follow Harman and Finnigan (2007) [HF07]_,
Harman and Finnigan (2008) [HF08]_, and Theeuwes et al. (2019) [T19]_.

The 30 levels have a step of 0.1 times the canopy height ``zh``
(should still output zh somewhere) ``dz = 0.1 * zh``.
However. if 3 x canopy height is less the 10 m steps of 0.3333 m are used:

.. code-block:: fortran

   IF ((3.*Zh) < 10.) THEN
   dz = 1./3.
   zarray = (/(I, I=1, nz)/)*dz...

Here ``nz = 30``.

.. note::

   All the diagnostic profiles (wind speed, temperature and humidity) are calculated
   from the forcing data down into the canopy.
   Therefore it is assumed that the forcing temperature and humidity
   are above the blending height.
