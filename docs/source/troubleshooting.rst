.. _Troubleshooting:

Troubleshooting
===============

How to report an issue of this manual?
--------------------------------------

    Please submit your issue via `our GitHub page. <https://github.com/Urban-Meteorology-Reading/SUEWS/issues>`_


How to join your email-list?
----------------------------

    Please join our email-list `here. <https://www.lists.reading.ac.uk/mailman/listinfo/met-suews>`_


How to create a directory?
--------------------------

    Please search the web using this phrase if you do not know how to
    create a folder or directory

How to unzip a file
-------------------

    Please search the web using this phrase if you do not know how to
    unzip a file


.. _A_text_editor:

A text editor
-------------

    A program to edit plain text files. If you search on the web
    using the phrase ‘text editor’ you will find numerous programs.
    These include for example, NotePad, EditPad, Text Pad etc

Command prompt
--------------

    From Start select run –type cmd – this will open a window. Change
    directory to the location of where you stored your files. The
    following website may be helpful if you do not know what a command
    prompt is: http://dosprompt.info/

Day of year [DOY]
-----------------

    January 1st is day 1, February 1st is day 32. If you search on the
    web using the phrase ‘day of year calendar’ you will find tables
    that allow rapid conversions. Remember that after February 28th DOY
    will be different between leap years and non-leap years.

ESTM output
-----------

First time steps of storage output could give NaN values during the
initial converging phase.

First things to Check if the program seems to have problems
-----------------------------------------------------------

-  Check the problems.txt file.
-  Check file options – in RunControl.nml.
-  Look in the output directory for the SS_FileChoices.txt. This allows
   you to check all options that were used in the run. You may want to
   compare it with the original version supplied with the model.
-  Note there can not be missing time steps in the data. If you need
   help with this you may want to checkout `UMEP`_

A pop-up saying “file path not found"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This means the program cannot find the file paths defined in
RunControl.nml file. Possible solutions:

-  Check that you have created the folder that you specified in
   RunControl.nml.
-  Check does the output directory exist?
-  Check that you have a single or double quotes around the
   FileInputPath, FileOutputPath and FileCode

====“%sat_vap_press.f temp=0.0000 pressure dectime”==== Temperature is
zero in the calculation of water vapour pressure parameterization.

-  You don’t need to worry if the temperature should be (is) 0°C.
-  If it should not be 0°C this suggests that there is a problem with
   the data.

%T changed to fit limits
~~~~~~~~~~~~~~~~~~~~~~~~

-  [TL =0.1]/ [TL =39.9] You may want to change the coefficients for
   surface resistance. If you have data from these temperatures, we
   would happily determine them.

%Iteration loop stopped for too stable conditions.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  [zL]/[USTAR] This warning indicates that the atmospheric stability
   gets above 2. In these conditions `MO
   theory <http://glossary.ametsoc.org/wiki/Monin-obukhov_similarity_theory>`__
   is not necessarily valid. The iteration loop to calculate the
   `Obukhov length <http://glossary.ametsoc.org/wiki/Obukhov_length>`__
   and `friction
   velocity <http://glossary.ametsoc.org/wiki/Friction_velocity>`__ is
   stopped so that stability does not get too high values. This is
   something you do not need to worry as it does not mean wrong input
   data.

“Reference to undefined variable, array element or function result”
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Parameter(s) missing from input files.

See also the error messages provided in problems.txt and warnings.txt

Email list
~~~~~~~~~~

-  SUEWS email list

`https://www.lists.reading.ac.uk/mailman/listinfo/met-suews <https://www.lists.reading.ac.uk/mailman/listinfo/met-suews>`__

-  UMEP email list

`https://www.lists.reading.ac.uk/mailman/listinfo/met-umep <https://www.lists.reading.ac.uk/mailman/listinfo/met-umep>`__


.. _`UMEP`: http://umep-docs.readthedocs.io/en/latest/index.html
