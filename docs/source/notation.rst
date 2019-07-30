Notation
========


.. glossary::

   *λF*
     Frontal area  index

   ΔQS
     Storage heat flux

   BLUEWS
     Boundary Layer  part of SUEWS

        .. figure:: /assets/img/Bluews_1.jpg
          :alt: Relation between BLUEWS and SUEWS

          Relation between BLUEWS and SUEWS

   CDD
     Cooling degree days

   GDD
     Growing degree days

   HDD
     Heating degree days

   Bldgs
     Building  surface

   CBL
     Convective  boundary layer

   DEM
     Digital   Elevation Model

   DSM
     Digital surface  model

   DTM
     Digital Terrain Model

   DecTr
     Deciduous trees and shrubs

   EveTr
     Evergreen trees and shrubs

   ESTM
     Element Surface Temperature Method (Offerle et al.,2005 [OGF2005]_)

   Grass
     Grass surface

   BSoil
     Unmanaged land and/or bare soil

   Runoff
      The water that drains freely off the impervious surface

   SoilStore
      The water stored in the underlying soil that infiltrates from the pervious surface

   L↓
     Incoming longwave radiation

   LAI
     Leaf area index

   LUMPS
     Local-scale Urban Meteorological Parameterization Scheme (Loridan et al. 2011 [L2011]_)

   MU
     Parameters which must be supplied and must be specific for the site/grid being run.

   MD
     Parameters which must be supplied and must be specific for the site/grid being run (but default values may be ok if these values are not known specifically for the site).

   O
     Parameters that are optional, depending on the model settings in `RunControl.nml`. Set any parameters that are not used/not known to ‘-999’.

   L
     Codes that are used to link between the input files. These codes are required but their values are completely arbitrary, providing that they link the input files in the correct way. The user should choose these codes, bearing in mind that the codes they match up with in column 1 of the corresponding input file must be unique within that file. Codes must be integers. Note that the codes must match up with column 1 of the corresponding input file, even if those parameters are not used (in which case set all columns except column 1 to ‘-999’ in the corresponding input file), otherwise the model run will fail.

   NARP
     Net All-wave Radiation Parameterization (Offerle et al. 2003 [O2003]_, Loridan et al. 2011 [L2011]_)

   OHM
     Objective Hysteresis Model (Grimmond et al. 1991 [G91OHM]_, Grimmond & Oke 1999a [GO99QS]_, 2002 [GO2002]_)

   Paved
     Paved surface

   |Qstar|
     Net all-wave radiation

   QE
     Latent heat flux

   QF
     Anthropogenic  heat flux

   QH
     Sensible heat  flux

   SOLWEIG
     The solar and longwave environmental irradiance geometry model (Lindberg et al. 2008 [FL2008]_,   Lindberg and Grimmond 2011 [FL2011]_)

   SVF
     Sky view factor

   *θ*
     Potential temperature

   tt
     Time step of data

   UMEP
     `Urban Multi-scale Environmental Predictor`_

   Water
     Water surface

   WATCH
     The WATCH project has produced a large number of data sets which should be of considerable use in regional and global studies of climate and water. see `WATCH webpage <http://www.eu-watch.org/data_availability>`__

   zi
     Convective boundary layer height

.. _Urban Multi-scale Environmental Predictor: http://umep-docs.readthedocs.io
