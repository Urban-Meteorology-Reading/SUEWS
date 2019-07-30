.. _introduction:

Introduction
============





Surface Urban Energy and Water Balance Scheme (**SUEWS**) (Järvi et al.
2011 [J11]_, Ward et al. 2016 [W16]_) is able to simulate the urban
radiation, energy and water balances using only commonly measured
meteorological variables and information about the surface cover. SUEWS
utilizes an evaporation-interception approach (Grimmond et al.
1991 [G91]_), similar to that used in forests, to model evaporation from
urban surfaces.


.. figure:: /assets/img/SUEWS_Overview_s.png
	:alt: Overview of SUEWS

	Overview of SUEWS




The model uses seven surface types: paved, buildings, evergreen
trees/shrubs, deciduous trees/shrubs, grass, bare soil and water. The
surface state for each surface type at each time step is calculated from
the running water balance of the canopy where the evaporation is
calculated from the Penman-Monteith equation. The soil moisture below
each surface type (excluding water) is taken into account.

Horizontal movement of water above and below ground level is allowed.
The user can specify the model time-step, but 5 min is strongly
recommended. The main output file is provided at a resolution of 60 min
by default. The model provides the radiation and energy balance
components, surface and soil wetness, surface and soil runoff and the
drainage for each surface. Timestamps refer to the end of the averaging
period.

Model applicability: SUEWS is a neighbourhood-scale or local-scale
model.

.. figure:: /assets/img/SUEWS_SurfaceWaterBalance_v2_xxs.jpg
	:alt: The seven surface types considered in SUEWS

	The seven surface types considered in SUEWS
