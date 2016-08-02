"""The ``reftools.interpretdq`` module contains functions for
interpreting DQ flags.

Examples
--------
>>> from reftools.interpretdq import ImageDQ, DQParser

Create a class instance for HST/ACS image using pre-defined DQ definitions.
Then, display the associated metadata and translation table:

>>> acsdq = ImageDQ.from_fits('j12345678q_flt.fits', ext=('DQ', 2))
>>> acsdq.parser.metadata
<Table length=1>
   key     val
  str10    str3
---------- ----
INSTRUMENT  ACS
>>> acsdq.parser.tab
<Table length=17>
DQFLAG ...                         LONG_DESCRIPTION
uint16 ...                              str74
------ ... ---------------------------------------------------------------
     0 ...                                                      Good pixel
     1 ...                                         Lost during compression
     2 ...                                          Replaced by fill value
     4 ... Bad detector pixel or beyond aperture or HRC upper-right defect
     8 ...              Masked by aperture feature or HRC occulting finger
    16 ...                                                       Hot pixel
   ... ...                                                             ...
   512 ...     Bad pixel in reference file (FLAT, polarizer, or dust mote)
  1024 ...                                                     Charge trap
  2048 ...                                          A-to-D saturated pixel
  4096 ...       Cosmic ray and detector artifact (AstroDrizzle, CR-SPLIT)
  8192 ...                                             Cosmic ray (ACSREJ)
 16384 ...                                        Manually flagged by user
 32768 ...                                                        Not used

Interpret DQ value for a single pixel at IRAF-style coordinate
``X=1457`` and ``Y=170``. Also display the original pixel value for comparison:

>>> acsdq.interpret_pixel(1457, 170)
<Table length=3>
DQFLAG SHORT_DESCRIPTION LONG_DESCRIPTION
uint16        str9            str74
------ ----------------- ----------------
    16               HOT        Hot pixel
    32               CTE         CTE tail
  1024              TRAP      Charge trap
>>> acsdq.data[169, 1456]
1072

Interpret DQ values for all pixels. Then, extract mask for interpreted
DQ value of 16 (hot pixel) and display pixel values for that mask:

>>> acsdq.interpret_all()
Parsing DQ flag(s)...
Done!
Run time: 2.822 s
N_FLAGGED: 550656/8388608 (6.564%)
FLAG=1    : 0 (0.000%)
FLAG=2    : 0 (0.000%)
FLAG=4    : 0 (0.000%)
FLAG=8    : 0 (0.000%)
FLAG=16   : 78901 (0.941%)
FLAG=32   : 78915 (0.941%)
FLAG=64   : 343917 (4.100%)
FLAG=128  : 20373 (0.243%)
FLAG=256  : 15202 (0.181%)
FLAG=512  : 26623 (0.317%)
FLAG=1024 : 5128 (0.061%)
FLAG=2048 : 15117 (0.180%)
FLAG=4096 : 0 (0.000%)
FLAG=8192 : 0 (0.000%)
FLAG=16384: 0 (0.000%)
FLAG=32768: 0 (0.000%)
>>> hotmask = acsdq.dq_mask(16)
>>> acsdq.data[hotmask]
array([528, 528, 528, ..., 560, 560, 560], dtype=int16)

Create a class instance for HST/WFC3 DQ parser (without image).
Then, interpret a given DQ value:

>>> wfc3dq = DQParser.from_instrument('WFC3')
>>> wfc3dq.interpret_dqval(16658)
<Table length=4>
DQFLAG SHORT_DESCRIPTION                 LONG_DESCRIPTION
uint16        str9                            str69
------ ----------------- ------------------------------------------------
     2            FILLED                           Replaced by fill value
    16               HOT                                        Hot pixel
   256         SATURATED              Full-well or A-to-D saturated pixel
 16384             CRMAX Pixel has more than max CRs, ghost, or crosstalk

"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from astropy.extern.six import iteritems
from astropy.extern.six.moves import map, zip

# STDLIB
import time
import warnings

# THIRD-PARTY
import numpy as np
from astropy.io import ascii, fits
from astropy.utils.data import get_pkg_data_filename
from astropy.utils.exceptions import AstropyUserWarning


__all__ = ['DQParser', 'ImageDQ']


class DQParser(object):
    """Class to handle parsing of DQ flags.

    **Definition Table**

    A "definition table" is an ASCII table that defines
    each DQ flag and its short and long descriptions.
    It can have optional comment line(s) for metadata,
    e.g.::

        # INSTRUMENT = HSTGENERIC

    It must have three columns:

    1. ``DQFLAG`` contains the flag value (``uint16``).
    2. ``SHORT_DESCRIPTION`` (string).
    3. ``LONG_DESCRIPTION`` (string).

    Example file contents::

        # INSTRUMENT = HSTGENERIC
        DQFLAG SHORT_DESCRIPTION LONG_DESCRIPTION
        0      "OK"              "Good pixel"
        1      "LOST"            "Lost during compression"
        ...    ...               ...

    The table format must be readable by ``astropy.io.ascii``.

    Parameters
    ----------
    definition_file : str
        ASCII table that defines the DQ flags (see above).

    Attributes
    ----------
    tab : ``astropy.table.Table``
        Table object from given definition file.

    metadata : ``astropy.table.Table``
        Table object from file metadata.

    """
    def __init__(self, definition_file):
        self._dqcol = 'DQFLAG'
        self._sdcol = 'SHORT_DESCRIPTION'
        self._ldcol = 'LONG_DESCRIPTION'
        self.tab = ascii.read(
            definition_file,
            names = (self._dqcol, self._sdcol, self._ldcol),
            converters = {self._dqcol: [ascii.convert_numpy(np.uint16)],
                          self._sdcol: [ascii.convert_numpy(np.str)],
                          self._ldcol: [ascii.convert_numpy(np.str)]})
        self.metadata = ascii.read(self.tab.meta['comments'], delimiter='=',
                                   format='no_header', names=['key', 'val'])

        # Ensure table has OK flag to detect good pixel
        self._okflag = 0
        if self._okflag not in self.tab[self._dqcol]:
            self.tab.add_row([self._okflag, 'OK', 'Good pixel'])

        # Sort table in ascending order
        self.tab.sort(self._dqcol)

        # Compile a list of flags
        self._valid_flags = self.tab[self._dqcol]
        self._total_flags = int(self._valid_flags.sum())

    @classmethod
    def from_instrument(cls, instrument):
        """Use pre-defined DQ flags for the given instrument
        (HST only for now).
        The data files are distributed with this module.

        Parameters
        ----------
        instrument : {'ACS', 'WFC3'}
            Currently not all instruments are supported.
            If the instrument is not found, a generic definition
            is used.

        Returns
        -------
        cls : `DQParser`
            An instance of this class.

        """
        if instrument is None:
            instrument = ''
        else:
            instrument = instrument.upper()

        if instrument == 'ACS':
            fname = 'data/dqflags_acs.txt'
        elif instrument == 'WFC3':
            fname = 'data/dqflags_wfc3.txt'
        else:
            warnings.warn(
                '{0} is not supported, using default'.format(instrument),
                AstropyUserWarning)
            fname = 'data/dqflags_hstgen.txt'

        return cls(get_pkg_data_filename(fname))

    def interpret_array(self, data, verbose=True):
        """Interpret DQ values for an array.

        .. warning::

            If the array is large and has a lot of flagged elements,
            this can be resource intensive.

        Parameters
        ----------
        data : ndarray
            DQ values.

        verbose : bool
            Print info to screen.

        Returns
        -------
        dqs_by_flag : dict
            Dictionary mapping each interpreted DQ value to indices
            of affected array elements.

        """
        if verbose:
            print('Parsing DQ flag(s)...')
            t_beg = time.time()

        data = np.asarray(data, dtype=np.int)  # Ensure int array
        dqs_by_flag = {}

        def _one_flag(vf):
            dqs_by_flag[vf] = np.where((data & vf) != 0)

        # Skip good flag
        list(map(_one_flag, self._valid_flags[1:]))

        if verbose:
            t_end = time.time()
            nbad = np.sum(data != self._okflag)
            ntot = data.size
            pbad = 100.0 * nbad / ntot
            print('Done!\nRun time: {0:.3f} s\nN_FLAGGED: {1}/{2} '
                  '({3:.3f}%)'.format(t_end - t_beg, nbad, ntot, pbad))
            for key in sorted(dqs_by_flag):
                nbad = len(dqs_by_flag[key][0])
                pbad = 100.0 * nbad / ntot
                print('FLAG={0:<5d}: {1} ({2:.3f}%)'.format(key, nbad, pbad))

        return dqs_by_flag

    def interpret_dqval(self, dqval):
        """Interpret DQ values for a single pixel.

        Parameters
        ----------
        dqval : int
            DQ value.

        Returns
        -------
        dqs : ``astropy.table.Table``
            Table object containing a list of interpreted DQ values and
            their meanings.

        """
        dqval = int(dqval)

        # Warn of non-standard DQ flags
        unknown_flags = dqval & ~self._total_flags
        if unknown_flags:
            warnings.warn(
                'Undefined DQ flags (sum={0}) found for input value {1}. '
                'Ignoring these flags...'.format(unknown_flags, dqval),
                AstropyUserWarning)

        # Good pixel, nothing to do
        if dqval == self._okflag:
            idx = 0

        # Find all the possible DQ flags
        else:
            idx = (dqval & self._valid_flags) != 0

        return self.tab[idx]


class ImageDQ(object):
    """Class to handle DQ flags in an image.

    Parameters
    ----------
    data : ndarray
        DQ data array to interpret.

    dqparser : `DQParser` or `None`
        DQ parser for interpretation. If not given, default is used.

    Attributes
    ----------
    data : ndarray
        Same as input.

    parser : `DQParser`
        DQ parser for interpretation.

    Raises
    ------
    ValueError
        Invalid image dimension.

    """
    def __init__(self, data, dqparser=None):
        ndim = 2
        data = np.asarray(data)

        if data.ndim != ndim:
            raise ValueError('Expected ndim={0} but data has '
                             'ndim={1}'.format(ndim, data.ndim))

        if 'int' not in data.dtype.name:
            warnings.warn('Data has dtype={0}, will be converted to '
                          'int...'.format(data.dtype), AstropyUserWarning)

        if dqparser is None or not isinstance(dqparser, DQParser):
            dqparser = DQParser.from_instrument(None)

        self.data = data
        self.parser = dqparser
        self._dqs_by_flag = None

    @classmethod
    def from_fits(cls, filename, ext=('DQ', 1), inskey='INSTRUME'):
        """Use data from given FITS image and create DQ parser based
        on header value.

        Parameters
        ----------
        filename : str
            Image filename.

        ext : str or int or tuple
            Extension identifier, as accepted by ``astropy.io.fits``.

        inskey : str
            Keyword value in ``PRIMARY`` header that defines the instrument.

        Returns
        -------
        cls : `ImageDQ`
            An instance of this class.

        """
        with fits.open(filename) as pf:
            data = pf[ext].data
            try:
                instrument = pf['PRIMARY'].header[inskey]
            except Exception as e:
                warnings.warn(
                    'Failed to read {0} from PRIMARY header, using default: '
                    '{1}'.format(inskey, str(e)), AstropyUserWarning)
                dqparser = None
            else:
                dqparser = DQParser.from_instrument(instrument)
        return cls(data, dqparser=dqparser)

    def interpret_pixel(self, x, y):
        """Construct a user-friendly table object with interpreted
        DQ values for a given pixel location (in IRAF format).

        Parameters
        ----------
        x, y : int
            IRAF X and Y pixel coordinates.

        Returns
        -------
        dqs : ``astropy.table.Table``
            Table object containing a list of interpreted DQ values and
            their meanings.

        """
        return self.parser.interpret_dqval(self.data[y - 1, x - 1])

    # This is not done in init because might be resource intensive.
    def interpret_all(self, verbose=True):
        """Interpret DQ values for all flagged pixels.
        Results are stored in ``self._dqs_by_flag`` internal cache.

        .. warning::

            If the image is large and has a lot of flagged pixels,
            this can be resource intensive.

        Parameters
        ----------
        verbose : bool
            Print info to screen.

        """
        self._dqs_by_flag = self.parser.interpret_array(
            self.data, verbose=verbose)

    def _check_cache(self):
        if self._dqs_by_flag is None:
            raise ValueError('Run interpret_all() method first!')

    def dq_mask(self, dqval):
        """Generate mask for the given interpreted DQ value.

        Parameters
        ----------
        dqval : int
            Interpreted DQ value.

        Returns
        -------
        mask : ndarray
            Boolean mask where the affected pixels are marked `True`.

        Raises
        ------
        ValueError
            Missing interpreter result or invalid DQ value.

        """
        self._check_cache()

        if dqval not in self._dqs_by_flag:
            raise ValueError('DQ={0} not found in {1}'.format(
                dqval, sorted(self._dqs_by_flag)))

        mask = np.zeros_like(self.data, dtype=np.bool)
        mask[self._dqs_by_flag[dqval]] = True

        return mask

    def pixlist(self):
        """Convert cached results from Python indices to list of
        IRAF-style ``(X, Y)`` pixel coordinates.

        Returns
        -------
        pixlist_by_flag : dict
            Dictionary mapping each interpreted DQ value to a list
            of IRAF coordinates.

        """
        self._check_cache()
        pixlist_by_flag = {}

        for key, idx in iteritems(self._dqs_by_flag):
            pixlist_by_flag[key] = list(zip(idx[1] + 1, idx[0] + 1))

        return pixlist_by_flag
