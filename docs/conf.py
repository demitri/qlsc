# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
from datetime import datetime
from unittest.mock import MagicMock

# import from local files rather than an installed version:
sys.path.insert(0, os.path.abspath('../source'))


# Exclude C extension since readthedocs can't build it.
# Ref: https://read-the-docs.readthedocs.io/en/latest/faq.html#i-get-import-errors-on-libraries-that-depend-on-c-modules
class Mock(MagicMock):
	@classmethod
	def __getattr__(cls, name):
		return MagicMock()
MOCK_MODULES = ['qlsc.q3c']
sys.modules.update((mod_name, Mock()) for mod_name in MOCK_MODULES)

import qlsc

# -- Project information -----------------------------------------------------

project = 'QLSC Python Module'
copyright = '2020, Demitri Muna'
#copyright = f'2015-{format(datetime.utcnow().year)} {author}, '
author = 'Demitri Muna'

# The full version, including alpha/beta/rc tags
release = qlsc.__version__


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
	'sphinx.ext.autodoc',
	'sphinx.ext.intersphinx',
	'sphinx.ext.graphviz',
	'sphinx.ext.inheritance_diagram',
	'sphinx.ext.viewcode',
	'sphinx.ext.todo',
	'sphinx.ext.napoleon', # must appear before 'sphinx-autodoc-typehints'
	'sphinx_autodoc_typehints'
]

#
# Extensions settings
#

# number of days to cache remotely downloaded 'inv' files (default = 5)
intersphinx_cache_limit = 5

intersphinx_mapping = {
	'astropy'    : ('https://docs.astropy.org/en/stable', None),
	'sqlalchemy' : ('https://docs.sqlalchemy.org/en/13/', None),
	'numpy'      : ('https://numpy.org/doc/stable', None)
}

todo_include_todos = True

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'alabaster'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
