Making IMPHTTAB Tables
======================

Making IMPHTTAB tables requires a connection to CDBS and that the ``PYSYN_CDBS``
environmental keyword point to CDBS.

Use :func:`reftools.mkimphttab.create_table` and provide the name of your
output file and the base observation mode from which to build a set observation
modes for which photometry keywords will be calculated.

Examples of a base observation mode include ``"acs,hrc"``, ``"wfc3,uvis1"``,
or ``"cos"``.

Examples
--------
>>> from reftools import mkimphttab
>>> mkimphttab.create_table('acs_wfc1', 'acs,wfc1', 'WFC1',
...                         'Mar 01 2002 00:00:00', clobber=True, verbose=False)
>>> mkimphttab.create_table('acs_sbc', 'acs,sbc', 'SBC',
...                         'Mar 01 2002 00:00:00', clobber=True)

.. currentmodule:: reftools.mkimphttab

.. autofunction:: reftools.mkimphttab.create_table
.. autofunction:: reftools.mkimphttab.create_table_from_table
.. autofunction:: reftools.mkimphttab.create_nicmos_table
.. autofunction:: reftools.mkimphttab.compute_values
.. autofunction:: reftools.mkimphttab.compute_synphot_values
