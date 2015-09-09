"""The ``reftools.parsedq`` module contains functions for parsing HST DQ flags.

Examples
--------
>>> from reftools import parsedq

Print defined DQ flags for ACS:

>>> parsedq.print_dq_dict('ACS')
INSTRUMENT: ACS
    0: OK
    1: Lost during compression
    2: Replaced by fill value
    4: Bad detector pixel or beyond aperture or HRC upper-right defect
    8: Masked by aperture feature or HRC occulting finger
   16: Hot pixel
   32: CTE tail
   64: Warm pixel (since Oct 8, 2004) or permanent hot pixel (before Oct 8, 2004)
  128: Bad column
  256: Full-well or A-to-D saturated pixel
  512: Bad pixel in reference file (FLAT, polarizer, or dust mote)
 1024: Charge traps
 2048: A-to-D saturated pixels
 4096: Cosmic rays and detector artifacts (AstroDrizzle, CR-SPLIT)
 8192: Cosmic rays (ACSREJ)
16384: Manually flagged by user
32768: Not used

Parse WFC3 DQ flag obtained from a single pixel:

>>> dqs = parsedq.parse_one_dq(16658, instrument='WFC3', verbose=True)
INSTRUMENT: WFC3
DQ: 16658
    2: Replaced by fill value
   16: Hot pixel
  256: Full-well or A-to-D saturated pixel
16384: Pixel has more than max CRs, ghost, or crosstalk

Parse DQ flags from ACS/WFC1 detector in the given image and overwrite existing
output file:

>>> parsedq.parse_one_image(
...     'j12345678q_crj.fits', 'results.txt', ext=('DQ', 2), clobber=True)
Parsing DQ flags...
Done! results.txt written

The output file ("results.txt") from above would look something like this::

    IMAGE: j12345678q_crj.fits[('DQ', 2)]
    INSTRUMENT: ACS

    X=11, Y=1, DQ=16
           16: Hot pixel

    X=34, Y=613, DQ=8192
         8192: Cosmic rays (ACSREJ)

    X=368, Y=613, DQ=128
          128: Bad column

    X=1457, Y=170, DQ=1040
           16: Hot pixel
         1024: Charge traps

"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# STDLIB
import os
import warnings
from collections import OrderedDict

# THIRD-PARTY
import numpy as np
from astropy.io import fits
from astropy.utils.exceptions import AstropyUserWarning

__all__ = ['parse_one_dq', 'parse_one_image', 'print_dq_dict']

# All instrument-specific dictionary must have same keys as None
dqval = {
    None: {
        0: 'OK',
        1: 'Lost during compression',
        2: 'Replaced by fill value',
        4: 'Bad detector pixel or beyond aperture',
        8: 'Masked by aperture feature',
        16: 'Hot pixel',
        32: 'CTE tail',
        64: 'Warm pixel',
        128: 'Bad column',
        256: 'Full-well or A-to-D saturated pixel',
        512: 'Bad pixel in reference file (FLAT)',
        1024: 'Charge traps',
        2048: 'A-to-D saturated pixels',
        4096: 'Cosmic rays and detector artifacts (AstroDrizzle, CR-SPLIT)',
        8192: 'Cosmic rays (CRREJ)',
        16384: 'Manually flagged by user',
        32768: 'Not used'},
    'ACS': {
        0: 'OK',
        1: 'Lost during compression',
        2: 'Replaced by fill value',
        4: 'Bad detector pixel or beyond aperture or HRC upper-right defect',
        8: 'Masked by aperture feature or HRC occulting finger',
        16: 'Hot pixel',
        32: 'CTE tail',
        64: ('Warm pixel (since Oct 8, 2004) or permanent hot pixel '
             '(before Oct 8, 2004)'),
        128: 'Bad column',
        256: 'Full-well or A-to-D saturated pixel',
        512: 'Bad pixel in reference file (FLAT, polarizer, or dust mote)',
        1024: 'Charge traps',
        2048: 'A-to-D saturated pixels',
        4096: 'Cosmic rays and detector artifacts (AstroDrizzle, CR-SPLIT)',
        8192: 'Cosmic rays (ACSREJ)',
        16384: 'Manually flagged by user',
        32768: 'Not used'},
    'WFC3': {
        0: 'OK',
        1: 'Lost during compression',
        2: 'Replaced by fill value',
        4: 'Bad detector pixel or beyond aperture',
        8: 'Masked by occulting bar or deviant IR zero-read pixel',
        16: 'Hot pixel',
        32: 'UVIS CTE tail or IR unstable pixel',
        64: 'Warm pixel',
        128: 'Bad bias value',
        256: 'Full-well or A-to-D saturated pixel',
        512: 'Bad pixel in reference file (FLAT)',
        1024: ('UVIS charge trap, SINK pixel, or CR spike detected during '
               'cridcalc IR'),
        2048: 'A-to-D saturated pixels or IR zero-read signal correction',
        4096: 'Cosmic rays and detector artifacts (AstroDrizzle)',
        8192: 'Rejected during image combination UVIS or IR CR rejection',
        16384: 'Pixel has more than max CRs, ghost, or crosstalk',
        32768: 'Manually flagged by user'}
    }
revkeys = sorted(dqval[None], reverse=True)[:-1]  # Bad pixels only


def parse_one_dq(dval, instrument=None, raise_error=True, verbose=False):
    """Parse DQ values for a single pixel.

    .. warning::

        If given value contains unknown DQ flag(s), this function
        returns inaccurate results.

    Parameters
    ----------
    dval : uint16
        DQ value.

    instrument : {`None`, 'ACS', 'WFC3'}
        Use instrument-specific DQ definition, if given.

    raise_error : bool
        Raise ``KeyError`` instead of falling back to default
        DQ definition.

    verbose : bool
        Print results to screen.

    Returns
    -------
    dqs : dict
        Parsed DQ values and their meanings.

    Raises
    ------
    KeyError
        Instrument not found.

    """
    dqs = OrderedDict()

    if instrument not in dqval:
        if raise_error:
            raise KeyError(
                '{0} is not supported: {1}'.format(instrument, sorted(dqval)))
        else:
            warnings.warn(
                '{0} is not supported, using default'.format(instrument),
                AstropyUserWarning)
            instrument = None

    # Good pixel, nothing to do
    if dval == 0:
        dqs[dval] = dqval[instrument][dval]

    # Find all the possible DQ flags
    else:
        leftover = dval
        while leftover > 0:
            for dq in revkeys:
                delta = leftover - dq
                if delta < 0:
                    continue
                dqs[dq] = dqval[instrument][dq]
                leftover = delta

    if verbose:
        print('INSTRUMENT: {0}\nDQ: {1}'.format(instrument, dval))
        for key in sorted(dqs):
            print('{0:>5d}: {1}'.format(key, dqs[key]))

    return dqs


def parse_one_image(image, outfile, ext=('DQ', 1), clobber=False):
    """Parse DQ values for all flagged pixels in the given
    DQ extension in a HST FITS image.

    .. warning::

        If the image is large and has a lot of flagged pixels,
        this can be time intensive and produce a large text file.

    Parameters
    ----------
    image : str
        Image filename.

    outfile : str
        Text file to store parsed DQ values for all flagged pixels.

    ext : str or int or tuple
        Extension identifier, as accepted by ``astropy.io.fits``.

    clobber : bool
        Clobber existing output file.

    Raises
    ------
    OSError
        Output file exists.

    """
    if os.path.exists(outfile) and not clobber:
        raise OSError(
            '{0} exists, use clobber=True to overwrite'.format(outfile))

    with fits.open(image) as pf:
        try:
            instrument = pf['PRIMARY'].header['INSTRUME']
        except KeyError:
            warnings.warn('Failed to extract instrument name from '
                          'primary header, using default', AstropyUserWarning)
            instrument = None
        data = pf[ext].data

    # Get locations of all flagged pixels
    pixel = np.where(data > 0)
    pixlist = np.dstack(pixel)

    print('Parsing DQ flags...')

    with open(outfile, 'w') as fout:
        fout.write('IMAGE: {0}[{1}]\nINSTRUMENT: {2}\n\n'.format(
            image, ext, instrument))

        # Write out coordinates in IRAF format
        for y, x in pixlist[0]:
            pixdq = data[y, x]
            fout.write('X={0}, Y={1}, DQ={2}\n'.format(x + 1, y + 1, pixdq))
            dqs = parse_one_dq(
                pixdq, instrument=instrument, raise_error=False, verbose=False)
            for key in sorted(dqs):
                fout.write('    {0:>5d}: {1}\n'.format(key, dqs[key]))
            fout.write('\n')

    print('Done! {0} written'.format(outfile))


def print_dq_dict(instrument=None):
    """Print currently defined DQ flags for the given instrument.

    Parameters
    ----------
    instrument : {`None`, 'ACS', 'WFC3'}
        Use instrument-specific DQ definition, if given.

    """
    if instrument not in dqval:
        warnings.warn('{0} is not supported, using default'.format(instrument),
                      AstropyUserWarning)
        instrument = None

    dqdict = dqval[instrument]
    print('INSTRUMENT: {0}'.format(instrument))
    for key in sorted(dqdict):
        print('{0:>5d}: {1}'.format(key, dqdict[key]))


# ----- Interactive unit test -----

def test_one_dq(input_flags):
    """Test with an arbitrary DQ combo. All input flags must be defined."""
    pixdq = sum(input_flags)
    dqs = parse_one_dq(pixdq, verbose=True)
    assert sorted(dqs) == sorted(input_flags)
