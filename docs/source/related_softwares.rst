.. _suews_related_softwares:

SUEWS-related Software
================================

.. _supy:

SuPy
----

`SuPy <https://supy.readthedocs.io/en/latest/>`_ is a Python-enhanced urban climate model with `SUEWS`_ as its computation core.

The scientific rigour in SuPy results is thus gurranteed by SUEWS (see :ref:`SUEWS publications <Recent_publications>` and :ref:`Parameterisations and sub-models within SUEWS`).

Meanwhile, the data analysis ability of SuPy is greatly enhanced by `the Python-based SciPy Stack <https://scipy.org>`_, notably `numpy`_ and `pandas`_.


.. _SUEWS: https://suews-docs.readthedocs.io/en/latest/
.. _numpy: https://www.numpy.org
.. _pandas: http://pandas.pydata.org/


- **How to get SuPy?**

  SuPy is available on all major platforms (macOS, Windows, Linux) for Python 3.5+
  via `PyPI <https://pypi.org/project/supy/>`_:

  .. code-block:: shell

    python3 -m pip install supy --upgrade

- **How to use SuPy?**

    * Please follow :ref:`Quickstart of SuPy` and :ref:`other tutorials <tutorial_index>`.

    * Please see :ref:`SuPy API <supy:api>` for usage details of SuPy functions.


.. _suews_umep:

SUEWS and UMEP
--------------


SUEWS can be run as a standalone model but also can be used within
`UMEP <http://umep-docs.readthedocs.io/en/latest/UMEP_Manual>`_. There are numerous
tools included within UMEP to help a user get started. The `SUEWS (Simple)`_
within UMEP is a fast way to start using SUEWS.

The version of SUEWS within UMEP is the complete model. Thus all options
that are listed in this manual are available to the user. In the UMEP
`SUEWS (Simple)`_ runs all options are set to values to allow intial exploration of the
model behaviour.


- Pre-Processor
	- Meteorological Data
		- `Prepare Existing Data`_
			Transforms meteorological data into UMEP format
		- `Download data (WATCH)`_
			Prepare meteorological dataset from `WATCH`


	- Spatial Data
		- `Spatial Data Downloader`_
			Plugin for retrieving geodata from online services suitable for various UMEP related tools
			- `LCZ Converter`_
			Conversion from Local Climate Zones (LCZs) in the WUDAPT database into SUEWS input data

	- Urban land cover
		- `Land Cover Reclassifier`_
			Reclassifies a grid into UMEP format land cover grid. Land surface models
		- `Land Cover Fraction (Point)`_
			Land cover fractions estimates from a land cover grid based on a specific point in space
		- `Land Cover Fraction (Grid)`_
			Land cover fractions estimates from a land cover grid based on a polygon grid

	- Urban Morphology
		- `Morphometric Calculator (Point)`_
			Morphometric parameters from a DSM based on a specific point in space
		- `Morphometric Calculator (Grid)`_
			Morphometric parameters estimated from a DSM based on a polygon grid
		- `Source Area Model (Point)`_
			Source area calculated from a DSM based on a specific point in space.

	- SUEWS input data
		- `SUEWS Prepare`_
			Preprocessing and preparing input data for the SUEWS model

- Processor
	- Anthropogenic Heat (|QF|)
		- `LQF`_
			Spatial variations anthropogenic heat release for urban areas
		- `GQF`_
			Anthropogenic Heat (|QF|).

	- Urban Energy Balance
		- `SUEWS (Simple)`_
			Urban Energy and Water Balance.
		- `SUEWS (Advanced)`_
			Urban Energy and Water Balance.

- Post-Processor
	- Urban Energy Balance
		- `SUEWS analyser`_
			Plugin for plotting and statistical analysis of model results from SUEWS simple and SUEWS advanced
	- Benchmark
		- `Benchmark System`_
			For statistical analysis of model results, such as SUEWS

.. _Prepare Existing Data: http://umep-docs.readthedocs.io/en/latest/pre-processor/Meteorological%20Data%20MetPreprocessor.html

.. _Download data (WATCH): http://umep-docs.readthedocs.io/en/latest/pre-processor/Meteorological%20Data%20Download%20data%20(WATCH).html

.. _Spatial Data Downloader: http://umep-docs.readthedocs.io/en/latest/pre-processor/Spatial%20Data%20Spatial%20Data%20Downloader.html

.. _LCZ Converter: http://umep-docs.readthedocs.io/en/latest/pre-processor/Spatial%20Data%20LCZ%20Converter.html

.. _Land Cover Reclassifier: http://umep-docs.readthedocs.io/en/latest/pre-processor/Urban%20Land%20Cover%20Land%20Cover%20Reclassifier.html

.. _Land Cover Fraction (Point): http://umep-docs.readthedocs.io/en/latest/pre-processor/Urban%20Land%20Cover%20Land%20Cover%20Fraction%20(Point).html

.. _Land Cover Fraction (Grid): http://umep-docs.readthedocs.io/en/latest/pre-processor/Urban%20Land%20Cover%20Land%20Cover%20Fraction%20(Grid).html

.. _Morphometric Calculator (Point): http://umep-docs.readthedocs.io/en/latest/pre-processor/Urban%20Morphology%20Morphometric%20Calculator%20(Point).html

.. _Morphometric Calculator (Grid): http://umep-docs.readthedocs.io/en/latest/pre-processor/Urban%20Morphology%20Morphometric%20Calculator%20(Grid).html

.. _Source Area Model (Point): http://umep-docs.readthedocs.io/en/latest/pre-processor/Urban%20Morphology%20Source%20Area%20(Point).html

.. _SUEWS Prepare: http://umep-docs.readthedocs.io/en/latest/pre-processor/SUEWS%20Prepare.html

.. _LQF: http://umep-docs.readthedocs.io/en/latest/processor/Urban%20Energy%20Balance%20LQ.html

.. _GQF: http://umep-docs.readthedocs.io/en/latest/processor/Urban%20Energy%20Balance%20GQ.html

.. _SUEWS (Simple): http://umep-docs.readthedocs.io/en/latest/processor/Urban%20Energy%20Balance%20Urban%20Energy%20Balance%20(SUEWS,%20simple).html

.. _SUEWS (Advanced): http://umep-docs.readthedocs.io/en/latest/processor/Urban%20Energy%20Balance%20Urban%20Energy%20Balance%20(SUEWS.BLUEWS,%20advanced).html

.. _SUEWS analyser: http://umep-docs.readthedocs.io/en/latest/post_processor/Urban%20Energy%20Balance%20SUEWS%20Analyser.html

.. _Benchmark System: http://umep-docs.readthedocs.io/en/latest/post_processor/Benchmark%20System.html



.. _Differences_between_SUEWS_LUMPS_and_FRAISE:


Differences between SUEWS, LUMPS and FRAISE
--------------------------------------------------------


The largest difference between LUMPS and SUEWS is that the latter
simulates the urban water balance in detail while LUMPS takes a simpler
approach for the sensible and latent heat fluxes and the water balance
(“water bucket”). The calculation of evaporation/latent heat in SUEWS is
more biophysically based. Due to its simplicity, LUMPS requires less
parameters in order to run. SUEWS gives turbulent heat fluxes calculated
with both models as an output.

Similarities and differences between LUMPS and SUEWS.

+--------------------+----------------------+-----------------------+
|                    | LUMPS                | SUEWS                 |
+====================+======================+=======================+
| Net all-wave       | Input or NARP        | Input or NARP         |
| radiation (Q*)     |                      |                       |
+--------------------+----------------------+-----------------------+
| Storage heat flux  | Input or from OHM    | Input or from OHM     |
| (ΔQS)              |                      |                       |
+--------------------+----------------------+-----------------------+
| Anthropogenic heat | Input or calculated  | Input or calculated   |
| flux (QF)          |                      |                       |
+--------------------+----------------------+-----------------------+
| Latent heat (QE)   | DeBruin and Holtslag | Penman-Monteith       |
|                    | (1982)               | equation2             |
+--------------------+----------------------+-----------------------+
| Sensible heat flux | DeBruin and Holtslag | Residual from         |
| (QH)               | (1982)               | available energy      |
|                    |                      | minus QE              |
+--------------------+----------------------+-----------------------+
| Water balance      | No water balance     | Running water balance |
|                    | included             | of canopy and water   |
|                    |                      | balance of soil       |
+--------------------+----------------------+-----------------------+
| Soil moisture      | Not considered       | Modelled              |
+--------------------+----------------------+-----------------------+
| Surface wetness    | Simple water bucket  | Running water balance |
|                    | model                |                       |
+--------------------+----------------------+-----------------------+
| Irrigation         | Only fraction of     | Input or calculated   |
|                    | surface area that is | with a simple model   |
|                    | irrigated            |                       |
+--------------------+----------------------+-----------------------+
| Surface cover      | Buildings, paved,    | Buildings, paved,     |
|                    | vegetation           | coniferous and        |
|                    |                      | deciduous             |
|                    |                      | trees/shrubs,         |
|                    |                      | irrigated and         |
|                    |                      | unirrigated grass     |
+--------------------+----------------------+-----------------------+

FRAISE Flux Ratio – Active Index Surface Exchange
-------------------------------------------------

FRAISE provides an estimate of mean midday (±3 h around solar noon)
energy partitioning from information on the surface characteristics and
estimates of the mean midday incoming radiative energy and anthropogenic
heat release. Please refer to Loridan and Grimmond (2012) [LG2012]_ for
further details.

+----------------+----------------+-----------------+-----------------+
| Topic          | FRAISE         | LUMPS           | SUEWS           |
+================+================+=================+=================+
| **Complexity** | Simplest:      |                 | More complex:   |
|                | FRAISE         |                 | SUEWS           |
+----------------+----------------+-----------------+-----------------+
| **Software     | R code         | Windows exe     | Windows exe     |
| provided:**    |                | (written in     | (written in     |
|                |                | Fortran)        | Fortran) -      |
|                |                |                 | other versions  |
|                |                |                 | available       |
+----------------+----------------+-----------------+-----------------+
| Applicable     | Midday (within | hourly          | 5               |
| period:        | 3 h of solar   |                 | min-hourly-annu |
|                | noon)          |                 | al              |
+----------------+----------------+-----------------+-----------------+
| Unique         | Calculates     | Radiation and   | Radiation,      |
| features:      | active surface | energy balances | energy and      |
|                |  and fluxes    |                 | water balance   |
|                |                |                 | (includes       |
|                |                |                 | LUMPS)          |
+----------------+----------------+-----------------+-----------------+
