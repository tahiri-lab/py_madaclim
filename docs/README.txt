################ START ################
# To generate the sphinx docs I first ran
sphinx-quickstart docs --no-sep -p py-madaclim -a "Simon Lalonde" -q

################ CONFIG ################
Then change the conf.py to have the rtd_theme + other extensions
# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'py-madaclim'
copyright = '2023, Simon Lalonde'
author = 'Simon Lalonde'

version = '0.1.0'
release = '0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.coverage',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.githubpages',
    'sphinx_rtd_theme',
    'sphinx.ext.intersphinx',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

language = 'en'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

################ GENERATE AUTO DOCSTRINGS ################
sphinx-apidoc -o docs/ src/py_madaclim

################ CUSTOMIZE THE INDEX.RST ################
# Modify the index.rst file

################ CUSTOMIZE THE MODULES + PACKAGES .RST ################
# Modify as needed for a better display

################ USE AUTO-BUILD ################
# Hot reloading with sphinx-autobuild for faster dev
cd docs
sphinx-autobuild . _build/htm

################ Deploy to githubpages ################