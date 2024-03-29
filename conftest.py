try:
    from pytest_astropy_header.display import (PYTEST_HEADER_MODULES,
                                               TESTED_VERSIONS)
except ImportError:
    PYTEST_HEADER_MODULES = {}
    TESTED_VERSIONS = {}

try:
    from reftools import __version__
except ImportError:
    __version__ = 'unknown'

# Uncomment and customize the following lines to add/remove entries
# from the list of packages for which version numbers are displayed
# when running the tests.
PYTEST_HEADER_MODULES['Astropy'] = 'astropy'
PYTEST_HEADER_MODULES['synphot'] = 'synphot'
PYTEST_HEADER_MODULES['stsynphot'] = 'stsynphot'
PYTEST_HEADER_MODULES.pop('Pandas')
PYTEST_HEADER_MODULES.pop('h5py')

TESTED_VERSIONS['reftools'] = __version__
