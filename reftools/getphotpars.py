"""This module includes utilities for calculating the photometry keywords
PHOTZPT, PHOTFLAM, PHOTPLAM, and PHOTBW for a given observation mode (obsmode)
and IMPHTTAB.
The calculations are performed in the same way here as they are in HSTCAL
pipeline.

To calculate a single set of keywords use the function :func:`get_phot_pars`.

If you are calculating for several obsmodes from a single IMPHTTAB file it's
best to use the `GetPhotPars` class. For example::

    get_phot = GetPhotPars(imphttab)
    for obs in obsmodes:
        photzpt, photflam, photplam, photbw = get_phot.get_phot_pars(obs)
        ...
    get_phot.close()

"""
# 4/2014  MLS Made this more general to accept IMPHTTAB files which
# have any number of extensions, and where the NAME of the phot keyword
# is the name of the extension (as already is the practice). This is done
# to accommodate new WFC3 tables with 5 extensions
#
# Right now, the only values that stsynphot computes are
# PHOTFLAM, PHOTPLAM, and PHOTBW.

import numpy as np
from astropy.io import fits

try:
    from . import _computephotpars
except ImportError:
    _computephotpars = None

__version__ = '0.1.2'
__vdate__ = '15-Apr-2014'

__all__ = ['ImphttabError', 'get_phot_pars', 'GetPhotPars']


class ImphttabError(Exception):
    """Class for errors associated with the imphttab file."""
    pass


def get_phot_pars(obsmode, imphttab):
    """
    Return PHOTZPT, PHOTFLAM, PHOTPLAM, and PHOTBW and any other
    keywords outlined in the table for specified obsmode
    and imphttab.

    Parameters
    ----------
    obsmode : str
        Complete obsmode string including any parameterized values.
        For example obsmodes, see `GetPhotPars`.

    imphttab : str
        Path and filename of IMPHTTAB reference file.

    Returns
    -------
    results_dict : dict
        Dictionary containing the photometric keyword results.

    """
    get_phot = GetPhotPars(imphttab)  # Not a context manager
    try:
        results_dict = get_phot.get_phot_pars(obsmode)
    finally:
        get_phot.close()
    return results_dict


class GetPhotPars:
    """Class to be used to get photometry parameters from a given
    IMPHTTAB reference file. Initialize with the name of an IMPHTTAB
    reference file and then call this class or the :meth:`get_phot_pars`
    method with a complete obsmode to get the photometry parameters.

    Example obsmodes are:

    * ``'acs,wfc1,f625w,f660n'``
    * ``'acs,wfc1,f625w,f814w,MJD#55000.0'``
    * ``'acs,wfc1,f625w,fr505n#5000.0,MJD#55000.0'``

    Parameters
    ----------
    imphttab : str
        Filename and path of IMPHTTAB reference file.

    Attributes
    ----------
    imphttab_name : str
        Filename and path of IMPHTTAB reference file. Same as input
        ``imphttab``.

    imphttab_fits : `astropy.io.fits.HDUList`
        Open ``HDUList`` object from ``imphttab``.

    """
    def __init__(self, imphttab):
        self.imphttab_name = imphttab
        self.imphttab_fits = fits.open(imphttab, 'readonly')

        # Check the extensions in the file, each contains a phot key.
        self._nextend = self.imphttab_fits[0].header['NEXTEND']
        self._parkeys = [self.imphttab_fits[ext].header['EXTNAME']
                         for ext in range(1, self._nextend + 1)]

        # Everything not in compute keys will just return the row value.
        self._compute_keys = ['PHOTFLAM', 'PHOTPLAM', 'PHOTBW']

    def get_phot_pars(self, obsmode):
        """Return the required keywords for the specified obsmode.

        Parameters
        ----------
        obsmode : str
            Observation mode string.

        Returns
        -------
        results_dict : dict
            Dictionary of photometry keyword results.
            Each key corresponds to one keyword and its value.

        """
        npars, strp_obsmode, par_dict = self._parse_obsmode(obsmode)
        par_struct = self._make_par_struct(npars, par_dict)
        result_dict = {}

        for par in self._parkeys:
            row = self._get_row(strp_obsmode, par)
            row_struct = self._make_row_struct(row, npars)

            # compute_value returns a float
            if par in self._compute_keys:
                result_dict[par] = self._compute_value(row_struct, par_struct)
            else:
                result_dict[par] = row_struct['results'][0]

        result_dict["PHOTZPT"] = self.imphttab_fits[0].header['PHOTZPT']
        return result_dict

    def close(self):
        """Close ``imphttab_fits`` attribute."""
        self.imphttab_fits.close()

    def _parse_obsmode(self, obsmode):
        """Return number of parameterized variables in ``obsmode`` and
        obsmode string with the values of the parameterized variables removed.
        Also returns a dictionary in which the keys are the names of the
        parameterized variables and the dictionary values are are the values
        of the parameterized variables.

        Parameters
        ----------
        obsmode : str
            Observation mode string.

        Returns
        -------
        npars : int
            Number of parameterized variables in ``obsmode``.

        strp_obsmode : str
            Obsmode with parameterized values removed and converted to lower
            case. Retains order of input.

        pars : dict
            Keys are parameterized variable names and values are the
            parameterized variable values. Keys are all lower case.

        """
        npars = obsmode.count('#')
        strp_obsmode = ''
        pars = {}

        for mode in obsmode.split(','):
            if '#' not in mode:
                strp_obsmode += mode + ','
            else:
                hashind = mode.index('#')
                i = hashind + 1
                s = mode[:i]
                strp_obsmode += s + ','
                pars[s.lower()] = float(mode[i:])

        # Remove trailing comma and convert to lower case.
        strp_obsmode = strp_obsmode[:-1].lower()
        return npars, strp_obsmode, pars

    def _get_row(self, obsmode, ext):
        """Return the `astropy.io.fits.FITS_rec` object corresponding to the
        table row from extension ``ext`` that has matching ``obsmode``.

        Parameters
        ----------
        obsmode : str
            Observation mode string.

        ext : str or int
            Specifier of FITS table extension from which to return a row.

        Returns
        -------
        row : `astropy.io.fits.FITS_rec`
            Row matching input ``obsmode``.

        Raises
        ------
        ImphttabError
            If ``obsmode`` does not appear in ``imphttab_fits`` or
            it appears multiple times in ``imphttab_fits``.

        """
        o = np.char.strip(np.char.lower(
            self.imphttab_fits[ext].data['obsmode']))
        w = np.where(o == obsmode.lower())

        if len(w[0]) == 0:
            raise ImphttabError(
                f'Obsmode {obsmode} does not appear in {self.imphttab_name} '
                f'extension {ext}')
        elif len(w[0]) > 1:
            raise ImphttabError(
                f'Obsmode {obsmode} appears multiple times in '
                f'{self.imphttab_name} extension {ext}')

        return self.imphttab_fits[ext].data[w]

    def _make_row_struct(self, row, npars):
        """Construct a dictionary corresponding to the ``PhtRow C`` structure
        used in the ``_computephotpars`` C-extension.

        Parameters
        ----------
        row : `astropy.io.fits.FITS_rec`
            A single row from an IMPHTTAB extension.

        npars : int
            Number of parameterized variables for this row.

        Returns
        -------
        row_struct : dict
            Dictionary corresponding to the ``PhtRow C`` structure used in the
            ``_computephotpars`` C-extension. Has the same keys as ``PhtRow``
            structure members. All keys are lower case.

        """
        row_struct = {}
        row_struct['obsmode'] = row['obsmode'][0]
        row_struct['datacol'] = row['datacol'][0]
        row_struct['parnames'] = []
        row_struct['parnum'] = npars
        row_struct['nelem'] = []
        row_struct['parvals'] = []

        for i in range(1, npars + 1):
            row_struct['parnames'].append(row[f'par{i}names'][0])
            row_struct['nelem'].append(row[f'nelem{i}'][0])
            row_struct['parvals'].append(row[f'par{i}values'][0].tolist())

        if npars == 0:
            row_struct['results'] = row[row['datacol'][0]]
            row_struct['telem'] = 1
        else:
            row_struct['results'] = row[row['datacol'][0]][0].tolist()
            row_struct['telem'] = len(row_struct['results'])

        return row_struct

    def _make_par_struct(self, npars, par_dict):
        """Construct a dictionary corresponding to (part of) the ``PhotPar``
        structure used in the ``_computephotpars`` C-extension. Not all
        members of the structure are used in the extension, so we supply the
        necessary ones here.

        Parameters
        ----------
        npars : int
            Number of parameterized variables for this obsmode.

        par_dict : dict
            Keys are parameterized variable names and values are the
            parameterized variable values. Keys are all lower case.
            As returned by `_parse_obsmode`.

        Returns
        -------
        par_struct : dict
            Dictionary with keys corresponding to some of the structure members
            of the ``PhotPar`` structure used in the ``_computephotpars``
            C-extension.

        """
        par_struct = {}
        par_struct['npar'] = npars
        par_struct['parnames'] = list(par_dict.keys())
        par_struct['parvals'] = [par_dict[k] for k in par_struct['parnames']]
        return par_struct

    def _compute_value(self, row_struct, par_struct):
        """Compute a photometry parameter based on an obsmode and a given row
        from the IMPHTTAB reference file.

        Parameters
        ----------
        row_struct : dict
            Dictionary corresponding to the ``PhtRow C`` structure used in the
            ``_computephotpar`` C-extension. Has the same keys as ``PhtRow``
            structure members. All keys are lower case.
            Should be the same as returned by ``_make_row_struct``.

        par_struct : dict
            Dictionary with keys corresponding to some of the structure members
            of the ``PhotPar`` structure used in the ``_computephotpars``
            C-extension. Should be the same as returned by
            ``_make_par_struct``.

        Returns
        -------
        result : float
            Result returned by ``_computephotpars.compute_value``.

        """
        if _computephotpars is None:
            raise ImportError('_computephotpars C-extension not built')
        return _computephotpars.compute_value(row_struct, par_struct)
