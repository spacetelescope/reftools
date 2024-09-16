2.2.0 (unreleased)
==================

- Dropped support for Python 3.6, 3.7, and 3.8. Minimum supported Python
  version is now 3.9. [#56, #71]

- Minimum supported version of ``astropy`` is now 5. [#71]

2.1.0 (2022-02-02)
==================

* Fixed FITS record handling in ``reftools.getphotpars``. [#54]
* Various packaging, infrastructure, and compatibility updates.

2.0.0 (2020-03-04)
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
