.. _input_converter:

SUEWS input converter
********************************

SUEWS input converter is a Python 3 script to convert input files between different versions based on pre-defined rules.

How to use
----------------
`Download the converter script and rule.csv <download_converter>` below, and specify these arguments in the script:

#. :code:`fromVer`: which version to convert from.
#. :code:`toVer`: which version to convert to.
#. :code:`fromDir`: where the input files are located.
#. :code:`toDir`: where the converted files are produced.


.. _download_converter:

Downloads
----------------

- SUEWS input converter in python

  :download:`SUEWS_TableConverter.py <SUEWS_TableConverter.py>`


- Rules for conversions between different SUEWS versions

  :download:`rules.csv <rules.csv>`

Description of rules
--------------------
The converter currently picks up the following types of actions:

#. Add: New entries or files to be added with default values.
#. Rename: Entries to be renamed from one version to another.
#. Delete: Entries to be deleted from one version to another.

.. note::

	For entries introduced in a version via a new file, the new file will be created to hold the new entries without extra delaration for new files.

The current available rules are listed below:

.. csv-table::
  :class: longtable
  :file: rules.csv
  :header-rows: 1
  :widths: auto
