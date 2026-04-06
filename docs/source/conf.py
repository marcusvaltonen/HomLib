# -*- coding: utf-8 -*-
import datetime

from sphinx_gallery.sorting import ExampleTitleSortKey


extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.coverage',
    'sphinx.ext.doctest',
    'sphinx.ext.viewcode',
    'sphinx.ext.extlinks',
    'sphinx.ext.napoleon',
    'sphinx_gallery.gen_gallery',
]


# Generate autodoc stubs with summaries from code
autosummary_generate = True

# Include Python objects as they appear in source files
autodoc_member_order = 'bysource'

# Default flags used by autodoc directives
autodoc_default_flags = ['members']

#autodoc_typehints = "description"

sphinx_gallery_conf = {
    # path to your examples scripts
    'examples_dirs': ['../../python/example'],
    # path where to save gallery generated examples
    'gallery_dirs': ['gallery'],
    'filename_pattern': '\.py',
    # Remove the "Download all examples" button from the top level gallery
    'download_all_examples': False,
    # Sort gallery example by file name instead of number of lines (default)
    'within_subsection_order': ExampleTitleSortKey,
    # directory where function granular galleries are stored
    'backreferences_dir': 'api/generated/backreferences',
    # Modules for which function level galleries are created.
    'doc_module': 'homlib',
    # Insert links to documentation of objects in the examples
    'reference_url': {'satoa': None},
    # Allow animations
    'matplotlib_animations': False,
}

# Always show the source code that generates a plot
plot_include_source = True
plot_formats = ['png']

# Sphinx project configuration
templates_path = ['_templates']
exclude_patterns = ['_build', '**.ipynb_checkpoints']
source_suffix = '.rst'
# The encoding of source files.
source_encoding = 'utf-8-sig'
master_doc = 'index'

# General information about the project
year = datetime.date.today().year
project = 'HomLib'
copyright = '{}, Marcus Valtonen Örnhag'.format(year)

# Version
version = 'latest'

# These enable substitutions using |variable| in the rst files
rst_epilog = """
.. |year| replace:: {year}
""".format(year=year)

html_last_updated_fmt = '%b %d, %Y'
html_title = 'HomLib'
html_short_title = 'HomLib'
#html_logo = '_static/logo.png'
#html_favicon = '_static/favicon.ico'
#html_static_path = ['_static']
html_extra_path = []
pygments_style = 'default'
add_function_parentheses = False
html_show_sourcelink = False
html_show_sphinx = True
html_show_copyright = True

# Theme config
html_theme = "sphinx_rtd_theme"
html_theme_options = {
    'logo_only': True,
    'display_version': True,
}

