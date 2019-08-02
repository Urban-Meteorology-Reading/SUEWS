# -*- coding: utf-8 -*-
#
# Configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/stable/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))
import os
import platform
import sys
# from datetime import datetime
from pathlib import Path

import pandas as pd
import nbsphinx


# -- processing code --------------------------------------------------------
# load all csv as a whole df
def load_df_csv(path_csv):
    if not path_csv.exists():
        print(str(path_csv), 'not existing!')
        sys.exit()

    # get list of all csv file names
    list_csv = list(path_csv.glob('*csv'))
    df_csv = pd.concat(
        {csv.stem: pd.read_csv(csv, skipinitialspace=True, quotechar='"')
            for csv in list_csv},
        sort=False)
    return df_csv


# retrieve description from rst files
def load_df_opt_desc(file_options):
    ser_opts = pd.read_csv(
        file_options,
        sep='\n', skipinitialspace=True)
    ser_opts = ser_opts.iloc[:, 0]
    ind_opt = ser_opts.index[ser_opts.str.contains('.. option::')]
    ser_opt_name = ser_opts[ind_opt].str.replace('.. option::', '').str.strip()
    ser_opt_desc = ser_opts[ind_opt + 2].str.strip()
    df_opt_desc = pd.DataFrame(
        {'desc': ser_opt_desc.values}, index=ser_opt_name.rename('option'))
    return df_opt_desc


# generate dataframe for a specific SUEWS table `csv_suews`
def gen_df_suews(df_csv, df_opt_desc, csv_suews):
    print('\t', csv_suews + '.csv')
    df_csv_suews = df_csv.loc[csv_suews].dropna(axis=1).copy()
    df_csv_suews.loc[:, 'No.'] = df_csv_suews.loc[:, 'No.'].astype(int)
    for ind, row in df_csv_suews.iterrows():
        var = row.loc['Column Name'].strip('`')
        if var in df_opt_desc.index:
            df_csv_suews.at[
                ind, 'Description'] = df_opt_desc.loc[var].values[0]
    return df_csv_suews


# save all re-generated CSV files to `path_csv`
def gen_csv_suews(path_csv):
    print('re-generating summary tables ...')
    # load all csv as a whole df
    df_csv = load_df_csv(path_csv)

    # retrieve description from rst files
    file_options = path_csv.parent / 'Input_Options.rst'
    df_opt_desc = load_df_opt_desc(file_options)

    list_csv_suews = df_csv.index.levels[0].to_series().filter(like='SUEWS')
    for csv_suews in list_csv_suews:
        df_csv_suews = gen_df_suews(df_csv, df_opt_desc, csv_suews)
        df_csv_suews.to_csv(path_csv / (csv_suews + '.csv'), index=False)

    return list_csv_suews


# -- Project information ----------------------------------------------------
project = u'SUEWS'
doc_name = u'SUEWS Documentation'
# today = datetime.today()
copyright = '2018 – 2019' + \
    ', micromet@University of Reading, led by Prof Sue Grimmond'
author = u'micromet@University of Reading, led by Prof Sue Grimmond'


# determine latest version and release
path_source = Path('.').resolve()
list_ver = sorted(
    [x.stem for x in list((path_source / 'version-history').glob('v*rst'))
     if 'version' not in x.stem])


# The short X.Y version
version = list_ver[-1]
# The full version, including alpha/beta/rc tags
release = list_ver[-1]

path_csv = path_source / 'input_files/SUEWS_SiteInfo/csv-table'
gen_csv_suews(path_csv)
# -- General configuration ---------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.mathjax',
    'sphinx.ext.ifconfig',
    'sphinx.ext.viewcode',
    # 'rinoh.frontend.sphinx',
    'sphinx.ext.autosectionlabel',
    # 'sphinxfortran.fortran_autodoc',
    # 'sphinxfortran.fortran_domain',
    'sphinx.ext.githubpages',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx.ext.extlinks',
    'recommonmark',
    'nbsphinx',
    'sphinx.ext.mathjax',
    'breathe',
    'exhale'

]

breathe_projects = {
    "SUEWS": "./doxygenoutput/xml"
}
breathe_default_project = "SUEWS"

exhale_args = {
    # These arguments are required
    "containmentFolder":     "./api",
    "rootFileName":          "library_root.rst",
    "rootFileTitle":         "API",
    "doxygenStripFromPath":  "..",
    # Suggested optional arguments
    "createTreeView":        True,
    # TIP: if using the sphinx-bootstrap-theme, you need
    # "treeViewIsBootstrap": True,
    "exhaleExecutesDoxygen": True,
    "exhaleUseDoxyfile" :    True,
    #"exhaleDoxygenStdin":    '''INPUT = ../../../SUEWS-SourceCode\n
    #                            GENERATE_HTML  = YES
    #                            '''
}

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
source_suffix = ['.rst', '.md']
# source_suffix = '.rst'

# fortran source code for `fortran_autodoc` and `fortran_domain`
fortran_src = [
    os.path.abspath('../fortran-src'),
]
fortran_ext = ['f90', 'F90', 'f95', 'F95']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The master toctree document.
master_doc = 'index'
master_doc_latex = 'index_latex'

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = None

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path .

exclude_patterns = ['_build', '**.ipynb_checkpoints']
# tags.add('html')
# if tags.has('html'):
#     exclude_patterns = ['references.rst']

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# default interpretation of `role` markups
default_role = 'any'

# some text replacement defintions
rst_prolog = """
.. |km^-1| replace:: km\ :sup:`-1`
.. |mm^-1| replace:: mm\ :sup:`-1`
.. |m^-1| replace:: m\ :sup:`-1`
.. |m^-2| replace:: m\ :sup:`-2`
.. |m^-3| replace:: m\ :sup:`-3`
.. |m^3| replace:: m\ :sup:`3`
.. |s^-1| replace:: s\ :sup:`-1`
.. |kg^-1| replace:: kg\ :sup:`-1`
.. |K^-1| replace:: K\ :sup:`-1`
.. |W^-1| replace:: W\ :sup:`-1`
.. |h^-1| replace:: h\ :sup:`-1`
.. |ha^-1| replace:: ha\ :sup:`-1`
.. |QF| replace:: Q\ :sub:`F`
.. |Qstar| replace:: Q\ :sup:`*`
.. |d^-1| replace:: d\ :sup:`-1`
.. |d^-2| replace:: d\ :sup:`-2`
.. |)^-1| replace:: )\ :sup:`-1`
.. |Recmd| replace:: **Recommended in this version.**
.. |NotRecmd| replace:: **Not recommended in this version.**
.. |NotAvail| replace:: **Not available in this version.**
.. |NotUsed| replace:: **Not used in this version.**

.. _GitHub page: https://github.com/Urban-Meteorology-Reading/SUEWS-Docs/issues/new?template=issue-report.md

.. only:: html

    .. note::

      Please report issues with the manual on the `GitHub page`_.


"""
# -- Options for HTML output -------------------------------------------------
# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'
html_theme_path = ["_themes"]

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
# html_theme_options = {}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
#
html_static_path = ['_static']
# html_context = {
#     'css_files': [
#         '_static/theme_overrides.css',  # override wide tables in RTD theme
#         ],
#      }

# Custom sidebar templates, must be a dictionary that maps document names
# to template names.
#
# The default sidebars (for documents that don't match any pattern) are
# defined by theme itself.  Builtin themes are using these templates by
# default: ``['localtoc.html', 'relations.html', 'sourcelink.html',
# 'searchbox.html']``.
#
# html_sidebars = {}
numfig = True
# html_logo = 'assets/img/SUEWS_LOGO.png'

# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'SUEWSdoc'


# -- Options for LaTeX output ------------------------------------------------
# this can be one of ['pdflatex', 'xelatex', 'lualatex', 'platex']
if platform.system() == 'Darwin':
    latex_engine = 'lualatex'
else:
    latex_engine = 'pdflatex'

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #
    'preamble': r'''
\usepackage[titles]{tocloft}
\usepackage{ragged2e}
\addto\captionsenglish{\renewcommand{\bibname}{References}}
\cftsetpnumwidth {1.25cm}\cftsetrmarg{1.5cm}
\setlength{\cftchapnumwidth}{0.75cm}
\setlength{\cftsecindent}{\cftchapnumwidth}
\setlength{\cftsecnumwidth}{1.25cm}
\newcolumntype{T}{L}
\setlength{\tymin}{40pt}
''',


    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

latex_show_pagerefs = False

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc,
     'SUEWS.tex',
     doc_name,
     author,
     'manual'),
]
# latex_logo = 'assets/img/SUEWS_LOGO.png'

# -- Options for manual page output ------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, 'suews', u'SUEWS Documentation',
     [author], 1),
]


# -- Options for Texinfo output ----------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (master_doc, 'SUEWS', u'SUEWS Documentation',
     author, 'SUEWS', 'One line description of project.',
     'Miscellaneous'),
]


# -- Options for Epub output -------------------------------------------------

# Bibliographic Dublin Core info.
epub_title = project
epub_author = author
epub_publisher = author
epub_copyright = copyright

# The unique identifier of the text. This can be a ISBN number
# or the project homepage.
#
# epub_identifier = ''

# A unique identification for the text.
#
# epub_uid = ''

# A list of files that should not be packed into the epub file.
epub_exclude_files = ['search.html']


# -- Extension configuration -------------------------------------------------
# rinoh_documents = [('index',            # top-level file (index.rst)
#                     'target',           # output (target.pdf)
#                     'Document Title',   # document title
#                     'John A. Uthor')]   # document author

# Fix for scrolling tables in the RTD-theme
# https://rackerlabs.github.io/docs-rackspace/tools/rtd-tables.html
def setup(app):
    app.add_stylesheet('theme_overrides.css')


# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'pandas': ('http://pandas.pydata.org/pandas-docs/stable/', None),
    'numpy': ('https://docs.scipy.org/doc/numpy/', None),
    'supy': ('https://supy.readthedocs.io/en/latest/', None),
}
