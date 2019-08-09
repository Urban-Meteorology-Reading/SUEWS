.. _coding_guideline:

Coding Guidelines
-------------------------

If you are interested in contributing to the code please contact Sue
Grimmond.

Coding
******

#. Core physics and calculatoin schemes of SUEWS are written in Fortran 90

#. Code is hosted in GitHub as private repository

#. Variables

   -  Names should be defined at least in one place in the code –
      ideally when defined
   -  Implicit None should be used in all subroutines
   -  Variable name should include units. e.g. Temp\_C, Temp\_K
   -  Output variable attributes should be provided in the TYPE
      structure defined in the ctrl_output module as follows:

       ::

           : TYPE varAttr
           : CHARACTER(len = 15) :: header ! short name in headers
           : CHARACTER(len = 12) :: unit   ! unit
           : CHARACTER(len = 14) :: fmt    ! output format
           : CHARACTER(len = 50) :: longNm ! long name for detailed description
           : CHARACTER(len = 1)  :: aggreg ! aggregation method
           : CHARACTER(len = 10) :: group  ! group: datetime, default, ESTM, Snow, etc.
           : INTEGER             :: level  ! output priority level: 0 for highest (defualt output)
           : END TYPE varAttr

#. Code should be written generally
#. Data set for testing should be provided
#. Demonstration that the model performance has improved when new code
   has been added or that any deterioration is warranted.
#. Additional requirements for modelling need to be indicated in the
   manual
#. All code should be commented in the program (with initials of who
   made the changes – name specified somewhere and institution)
#. The references used in the code and in the equations will be
   collected to a webpage
#. Current developments that are being actively worked on


Testing
*******

#. The testing of SUEWS is done using Python 3
#. The following tests are done for each release of SUEWS:

  #. Working status of `all physics schemes <scheme_options>`
  #. Year-grid looping logic
  #. Identity of output results with internal test dataset

Please use pre-defined `make test` option to check if your code can pass all tests or not.
If not, the correctness of added code should be justified with caution.



Preparation of SUEWS Manual
***************************

#. The SUEWS manual is written in `reStructuredText (aka rst) <http://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`_ with a `Sphinx <http://www.sphinx-doc.org/>`_ flavour
#. The SUEWS manual is hosted by `readthedocs.org <https://www.readthedocs.org>`_
#. CSV tables used in following pages are automatically generated from the *Description* field in `Input_Options` by each build, so **DON'T** manually edit them as your edits will be swiped automatically:

  * `SUEWS_AnthropogenicEmission.txt`
  * `SUEWS_BiogenCO2.txt`
  * `SUEWS_Conductance.txt`
  * `SUEWS_Irrigation.txt`
  * `SUEWS_NonVeg.txt`
  * `SUEWS_OHMCoefficients.txt`
  * `SUEWS_Profiles.txt`
  * `SUEWS_SiteSelect.txt`
  * `SUEWS_Snow.txt`
  * `SUEWS_Soil.txt`
  * `SUEWS_Veg.txt`
  * `SUEWS_Water.txt`
  * `SUEWS_WithinGridWaterDist.txt`

F2PY tips
*********

This includes several **DON'T**'s
that have never been mentioned by F2PY docs:

1. DON'T mix comments as lines into argument list of Fortran subroutines/functions:

  DONT:

  .. code-block:: fortran

      subroutine(&
      ! DONT DO this
      args&
      )

  OK:

  .. code-block:: fortran

      subroutine(&
      args& ! OK this way
      )

2. DON'T end a subroutine as ``ENDSUBROUTINE``.
Instead, leave a space in between
to form ``END SUBROUTINE``.
Otherwise, the subroutines won't be correctly
parsed and picked up by F2PY.