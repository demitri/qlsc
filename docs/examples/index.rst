Worked Examples
===============

Each page here is a complete problem worked end to end. The code is meant to be
copied and run: every block of output shown was produced by actually running the
code above it, on an Apple M1 with Python 3.13, NumPy 2.5, and qlsc 1.1.0.

Most of the catalogs are synthetic, generated from a fixed random seed so that
your numbers will match the ones printed here. Where synthetic data is used it is
said so plainly; :doc:`bright_star_mask` is the exception and uses a real
catalog fetched from VizieR.

The figures are drawn by ``docs/generate_figures.py`` in the repository, which
uses the same seeds; none of them are hand-made.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Example
     - What it demonstrates
   * - :doc:`cone_search`
     - ``radial_query`` with keys and distances; how much faster an index is
       than scanning
   * - :doc:`bright_star_mask`
     - a real catalog, keys as a join back to source data, many small cone
       searches
   * - :doc:`survey_footprint`
     - ``polygon_query`` for convex and concave footprints; great-circle edges
       and size limits
   * - :doc:`galaxy_aperture`
     - ``ellipse_query`` for an inclined disc, and what a circular aperture
       costs you
   * - :doc:`pixel_regions`
     - ``ipix_neighbors``, region growing, and the pixel-area correction
   * - :doc:`density_map`
     - pixel counts as a sky map, and rolling counts up the quadtree

.. toctree::
   :maxdepth: 1

   cone_search
   bright_star_mask
   survey_footprint
   galaxy_aperture
   pixel_regions
   density_map
