#import os
#import sys

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'ROMpEIT'
copyright = '2024, Matthew Walker and Leandro Beltrachini'
author = 'Matthew Walker and Leandro Beltrachini'
release = 'v0.1.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

#sys.path.insert(0, os.path.abspath(".."))

extensions = [
         "sphinxcontrib.matlab",
         "sphinx.ext.autodoc",
         "sphinx_rtd_theme",
         "sphinx.ext.viewcode",
         "sphinx.ext.napoleon"
     ]

primary_domain = "mat"
matlab_src_dir = "../functions"

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
latex_engine = "xelatex"
