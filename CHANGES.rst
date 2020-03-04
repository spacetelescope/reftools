2.0.0 (unreleased)
==================

* ``interpretdq`` has new option called ``origin`` to return pixel values
  in either Python 0-indexed or IRAF 1-indexed system, where applicable. [#35]
* IMPHTTAB generation now uses ``stsynphot``, not PySynphot. [#33]
* Compatibility with ``astropy`` 4.0. [#31]
* Dropped Python 2 support. Minimum supported Python version is 3.6. [#31]
* Dropped TEAL support for ``tdspysyn``. [#31]
* Dropped PyRAF support, which means ``synpysyncomp`` is no longer
  available. [#31]
* Removed unused ``wtraxyutils`` and ``test_small_dgeo`` modules. [#31]

Older versions
==============

See GitHub release notes.
