.. _index_page:

SUEWS: Surface Urban Energy and Water Balance Scheme
----------------------------------------------------

.. image:: http://readthedocs.org/projects/suews/badge/?version=latest
    :target: https://suews.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status


- **How to get SUEWS?**

  - **Latest release:**

     The **latest formal** release of SUEWS is `new_latest` and can be downloaded via `our Zenodo repository`_ (a sample input dataset is included in the release archive).

  - **Previous releases:**

     Previous releases can be downloaded via `our GitHub page`_.



- **How to use SUEWS?**

  - **For existing users:**

  .. epigraph::

    Overview of changes in this version, see :ref:`new_latest`.
    If these changes impact your existing simulations, please see appropriate parts of the manual. It may be necessary to adapt some of your input files for for the current version.

    .. tip::

        A helper python script, `SUEWS table converter <input_converter>`, is provided to help facilitate the conversion of input files between different SUEWS versions.

    Additionally, the manuals for previous versions can be accessed in respective sections under `version_history`.


  - **For new users:**

  .. epigraph::

    Before performing SUEWS simulations, new users should read the overview :ref:`introduction`, then follow the steps in `Preparing_to_run_the_model` to prepare `input files <input_files>` for SUEWS.

    Note there are tutorials learning about running SUEWS available :ref:`the tutorial. <tutorials_suews>`


- **How has SUEWS been used?**

.. epigraph::

  The scientific details and application examples of SUEWS can be found in `Recent_publications`.

.. _cite_suews:

- **How to cite SUEWS?**

  .. tip::

      Visit the repositories below for different citation styles.

  - Software:

       Sun Ting, Järvi Leena, Grimmond Sue, Lindberg Fredrik, Li Zhenkun, Tang Yihao, Ward Helen: (2019, February 21). SUEWS: Surface Urban Energy and Water Balance Scheme (Version 2018c). Zenodo. |doi_software|


  - Manual:

       Sun Ting, Järvi Leena, Grimmond Sue, Lindberg Fredrik, Li Zhenkun, Tang Yihao, Ward Helen: (2019, February 21). SUEWS Documentation (Version 2018c). Zenodo. |doi_docs|




- **How to support SUEWS?**

  #. `Cite SUEWS <cite_suews>` appropriately in your work.
  #. Contribute to the `development <Development_Suggestions_Support>`.
  #. Report issues via the `GitHub page`_.
  #. Provide `suggestions and feedback <Development_Suggestions_Support>`.


.. _our GitHub page: https://urban-meteorology-reading.github.io/SUEWS
.. _our Zenodo repository: https://zenodo.org/record/2274254
.. _this form: `dowload form`_
.. _dowload form: http://micromet.reading.ac.uk/software/

.. |doi_software| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.2574410.svg
    :target: https://doi.org/10.5281/zenodo.2574410

.. |doi_docs| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.2574997.svg
    :target: https://doi.org/10.5281/zenodo.2574997

.. toctree::
   :maxdepth: 2
   :numbered:
   :hidden:

   introduction
   parameterisations-and-sub-models
   prepare-to-run-the-model
   input_files/input_files
   output_files/output_files
   troubleshooting
   recent-publications
   related_softwares
   sub-tutorials/tutorials
   development/development
   benchmark/benchmark_report
   version-history/version-history
   acknowledgement
   notation
   references
   api
