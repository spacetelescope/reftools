#! /usr/bin/env python
"""Convert a COS or STIS TDS file to ``stsynphot`` throughput files.

This program reads a time-dependent sensitivity (TDS) file for either
COS or STIS and writes one throughput table in ``stsynphot`` format for each
input row.

.. note::

    A row number may be specified in order to create just a single throughput
    table for that row in the input TDS table.

To run this task from within Python::

    >>> from reftools import tdspysyn
    >>> tdspysyn.tdsToPysynphot(
    ...     'u7d20377l_tds.fits', 'xyz',
    ...     1100, 1700, 100, 55927, None, None)  # doctest: +SKIP

To run this task from the command line::

    tdspysyn xyz_tds.fits thru 1100 1700 100 55927

"""
# Phil Hodge, STScI, February 2011

import math
import os
import warnings

import numpy as np
from astropy.io import fits

__taskname__ = 'tdspysyn'
__version__ = '0.1'
__vdate__ = '2011-02-18'
__all__ = ['tdsToPysynphot']

DAYS_PER_YEAR = 365.25


def main():
    """Command line driver."""
    import argparse

    parser = argparse.ArgumentParser(
        prog=__taskname__,
        description='Convert a COS or STIS TDS file to stsynphot throughput files.')  # noqa
    parser.add_argument('input', type=str, help='Input TDS file name')
    parser.add_argument('outroot', type=str, help='Prefix for output names')
    parser.add_argument('min_wl', type=float,
                        help='Minimum wavelength for stsynphot files')
    parser.add_argument('max_wl', type=float,
                        help='Maximum wavelength for stsynphot files')
    parser.add_argument('dwl', type=float,
                        help='Wavelength increment for stsynphot files')
    parser.add_argument('thru_mjd', type=float,
                        help='Time (MJD) for THROUGHPUT column')
    parser.add_argument('last_mjd', type=float, default=None, help='Last MJD')
    parser.add_argument('row', type=int, default=None,
                        help='One-indexed row number in TDS table')
    parser.add_argument(
        '--version', action='version',
        version=f'{__taskname__} {__version__} ({__vdate__})')
    args = parser.parse_args()

    tdsToPysynphot(args.input, args.outroot, args.min_wl, args.max_wl,
                   args.dwl, args.thru_mjd, args.last_mjd, args.row)


def tdsToPysynphot(input, outroot, min_wl, max_wl, dwl, thru_mjd,
                   last_mjd=None, row=None):
    """Convert a TDS table to ``stsynphot`` format.

    Parameters
    ----------
    input : str
        The name of the input TDS table.

    outroot : str
        A root name for constructing the name(s) of the output ``stsynphot``
        file(s). For COS data, the full file name will be ``outroot`` plus
        the grating name, the segment or stripe name, the aperture name,
        and ``.fits``.  The grating, segment, and aperture names will be in
        lower case, and each will be preceded by ``_``.  For STIS data,
        the full file name will be ``outroot`` plus ``_``, the grating name
        (lower case), and ``.fits``.

    min_wl : float
        The minimum wavelength (Angstroms) for the WAVELENGTH column in
        the output ``stsynphot`` table.

    max_wl : float
        The maximum wavelength (Angstroms) for the WAVELENGTH column in
        the output ``stsynphot`` table.

    dwl : float
        The increment (Angstroms) from row to row in the WAVELENGTH column
        in the output ``stsynphot`` table.

    thru_mjd : float
        The time (MJD) to use for the THROUGHPUT column.  This column
        contains the default throughput values that ``stsynphot`` will use
        if no date was specified.
        This is the one time value that will not be rounded to an int.

    last_mjd : float or `None`
        A time (MJD) in the future. If `None`, no ``MJD#`` column
        will be created for this value; this would be appropriate, for
        example, if the input TDS table already has a time value that is
        far enough in the future.

    row : int or `None`
        If a row number (1-indexed) was specified, a ``stsynphot`` table
        will be written for just that one row of the input TDS table.
        If `None` (the default), a ``stsynphot`` table will be written
        for each row of the TDS table.

    """
    # convert row to zero indexing
    row -= 1

    tds = initTds(input, outroot, min_wl, max_wl, dwl, thru_mjd, last_mjd, row)

    for row in tds.rowlist:
        tds.makeFilename(row)
        time = tds.getTimes(row)
        wavelength = tds.getWavelengths(min_wl, max_wl, dwl)
        hdu = tds.createTable(time, wavelength)
        tds.assignThroughputs(hdu, time, row)
        tds.createFile(hdu, row)


def initTds(input, outroot, min_wl, max_wl, dwl, thru_mjd, last_mjd, row):
    """Create a ``ConvertTds`` object for either COS or STIS data.

    Parameters
    ----------
    input : str
        The name of the input TDS table.

    outroot : str
        A root name for constructing the name(s) of the output ``stsynphot``
        file(s).

    min_wl : float
        The minimum wavelength (Angstroms) for the WAVELENGTH column in
        the output ``stsynphot`` table.

    max_wl : float
        The maximum wavelength (Angstroms) for the WAVELENGTH column in
        the output ``stsynphot`` table.

    dwl : float
        The increment (Angstroms) from row to row in the WAVELENGTH column
        in the output ``stsynphot`` table.

    thru_mjd : float
        The time (MJD) to use for the THROUGHPUT column. This column
        contains the default throughput values that ``stsynphot`` will use
        if no date was specified.

    last_mjd : float or `None`
        A time (MJD) in the future. If `None`, no ``MJD#`` column
        will be created for this value.

    row : int or `None`
        Row number (0-indexed) to read from the input TDS table, or `None`
        if a ``stsynphot`` table for each TDS table row should be written.

    """
    with fits.open(expandFilename(input)) as ifd:
        filetype = ifd[0].header.get('filetype', 'missing')
        instrument = ifd[0].header.get('instrume', 'missing')
        if row is not None:
            nrows = len(ifd[1].data)
            if row < 0 or row >= nrows:
                raise RuntimeError(f'row is out of range; there are {nrows} '
                                   'rows in the TDS table')

    if (filetype != "missing" and
            filetype != "TIME DEPENDENT SENSITIVITY TABLE"):
        warnings.warn(f'FILETYPE = {filetype}')

    if instrument == "COS":
        tds = CosTds(input, outroot, min_wl, max_wl, dwl,
                     thru_mjd, last_mjd, row)
    elif instrument == "STIS":
        tds = StisTds(input, outroot, min_wl, max_wl, dwl,
                      thru_mjd, last_mjd, row)
    else:
        raise RuntimeError(f'INSTRUME is {instrument}; must be COS or STIS')

    return tds


def expandFilename(filename):
    """Get the real file name.

    Parameters
    ----------
    filename : str
        A file name.

    Returns
    -------
    real_file_name : str
        The real file name.

    """
    fname = filename
    done = False
    count = 0
    MAX_COUNT = 100
    while not done:
        temp = os.path.expandvars(fname)  # $stuff/file
        count += 1
        if temp == fname:
            done = True
        fname = temp
        if count >= MAX_COUNT:
            break
    if not done:
        raise RuntimeError(f'{MAX_COUNT} iterations exceeded while expanding '
                           f'variables in file name {filename}')
    fname = os.path.abspath(fname)            # ../file
    fname = os.path.expanduser(fname)         # ~/file
    real_file_name = os.path.normpath(fname)  # remove redundant strings

    return real_file_name


def toMjd(year, month, day):
    """Convert calendar date (UTC) to Modified Julian Date (MJD).

    This is based on expression 12.92-1 in the Explanatory Supplement
    to the Astronomical Almanac, 1992, p 604.

    Parameters
    ----------
    year : int
        Four-digit year.

    month : int
        One-indexed month of the year.

    day : int or float
        One-indexed day of the month, optionally including a fractional
        part. The time scale is assumed to be UTC.

    Returns
    -------
    mjd : float
        Modified Julian Date.

    """
    # In the ES, this term was given as (month - 14) / 12, with integer
    # division.  The quotient (before truncation) is always negative.
    # In C or Fortran, the result would be truncated toward zero; however,
    # in Python the result is the floor.  To work around this difference,
    # the expression was rewritten as follows.
    a = (month + 9) // 12 - 1
    imjd = ((1461 * (year + 4800 + a)) // 4 +
            (367 * (month - 2 - 12 * a)) // 12 -
            (3 * ((year + 4900 + a) // 100)) // 4)
    return float(imjd) - 2432076. + day


class ConvertTds:
    """Class to handle TDS conversion."""

    def __init__(self, input, outroot, min_wl, max_wl, dwl, thru_mjd,
                 last_mjd, row):
        self.input = expandFilename(input)
        self.data = None                # will be the input data block
        self.phdr_tds = None            # assigned below
        self.outroot = expandFilename(outroot)
        self.output = self.outroot      # constructed later from outroot, etc.
        self.orig_input = input         # save for a history record
        self.orig_outroot = outroot     # save for a history record
        self.min_wl = min_wl
        self.max_wl = max_wl
        self.dwl = dwl
        self.thru_mjd = thru_mjd
        self.last_mjd = last_mjd        # may be None
        self.row = row                  # save for a history record

        with fits.open(self.input) as ifd:
            self.phdr_tds = ifd[0].header.copy()
            self.useafter = ifd[0].header.get('useafter', 'missing')
            self.useafter_mjd = self.useafterMjd()
            self.ref_time = ifd[1].header.get('ref_time', 'missing')
            self.data = ifd[1].data.copy()
            if row is None:
                self.rowlist = list(range(len(self.data)))
            else:
                self.rowlist = [row]

    def makeFilename(self, row):
        raise NotImplementedError('To be defined in a subclass')

    def getTimes(self, row):
        """Create the array of times.

        Parameters
        ----------
        row : int
            Current row number (0-indexed) in input TDS table.

        Returns
        -------
        time : array_like
            MJD at useafter date, times from TDS table, a final (late) MJD
            (if last_mjd was specified).
            The values will be rounded to the nearest integer, sorted by
            date, and duplicates will be removed.

        """
        time_tds = self.data.field('time')[row]
        nt_tds = self.data.field('nt')[row]            # from TDS table
        # add one for the useafter date, and possibly add one for last_mjd
        if self.last_mjd is None:
            nt_pysyn = nt_tds + 1
        else:
            nt_pysyn = nt_tds + 2
        time = np.zeros(nt_pysyn, dtype=np.float64)    # for stsynphot table

        # NOTE: round the MJD values to the nearest integer
        time[0] = round(self.useafter_mjd)
        for i in range(nt_tds):
            time[i+1] = round(time_tds[i])             # from TDS table
        if self.last_mjd is not None:
            time[-1] = round(self.last_mjd)

        sorted_times = np.sort(time)

        # copy sorted, unique times back to time
        k = 0                                           # initial value
        time[k] = sorted_times[0]
        last_time = sorted_times[0]
        for i in range(1, nt_pysyn):
            if sorted_times[i] > last_time:
                k += 1
                time[k] = sorted_times[i]
                last_time = sorted_times[i]
        nelem = k + 1

        return time[0:nelem]

    def getWavelengths(self, min_wl, max_wl, dwl):
        """Create the array of wavelengths for the ``stsynphot`` table.

        Parameters
        ----------
        min_wl : float
            Minimum wavelength (Angstroms) for ``stsynphot`` table.

        max_wl : float
            Maximum wavelength (Angstroms) for ``stsynphot`` table.

        dwl : float
            Wavelength increment (Angstroms) for ``stsynphot`` table.

        Returns
        -------
        wavelength : array_like
            Array of wavelengths (float64).

        """
        nwl = int(math.ceil((max_wl - min_wl) / dwl)) + 1
        return np.arange(nwl, dtype=np.float64) * dwl + min_wl

    def createTable(self, time, wavelength):
        """Create the output bintable header/data unit.

        All the columns will be defined, and the values in ``wavelength``
        will be copied to the WAVELENGTH column.

        Parameters
        ----------
        time : array_like
            The array of times returned by :meth:`getTimes`.

        wavelength : array_like
            The array of wavelengths returned by :meth:`getWavelengths`.

        Returns
        -------
        hdu : `astropy.io.fits.BinTableHDU`
            The header/data unit containing the table.

        """
        nrows = len(wavelength)
        col = []
        col.append(fits.Column(name='WAVELENGTH', format='1D',
                               unit='ANGSTROM', disp='G10.7'))
        col.append(fits.Column(name='THROUGHPUT', format='1E',
                               unit='TRANSMISSION', disp='G15.7'))
        col.append(fits.Column(name='ERROR', format='1E',
                               unit='TRANSMISSION', disp='G15.7'))
        for i in range(len(time)):
            colname = f'MJD#{time[i]:.1f}'
            col.append(fits.Column(name=colname, format='1E',
                                   unit='TRANSMISSION', disp='G15.7'))
        cd = fits.ColDefs(col)
        hdu = fits.BinTableHDU.from_columns(cd, nrows=nrows)
        hdu.header.update('thru_mjd', self.thru_mjd,
                          comment='Time (MJD) for THROUGHPUT column')
        hdu.data.field('wavelength')[:] = wavelength
        return hdu

    def assignThroughputs(self, hdu, time, row):
        """Assign values to the throughput columns.

        The values in the throughput column for the useafter date will be
        assigned the value 1. Other throughput columns will be assigned
        values computed from the input TDS table.

        Parameters
        ----------
        hdu : `astropy.io.fits.BinTableHDU`
            The header/data unit containing the table.

        time : array_like
            The array of times: MJD at useafter date, times from TDS
            table, and possibly a final (late) MJD.

        row : int
            Current row number (0-indexed) in input TDS table.

        """
        wavelength = hdu.data.field('wavelength')      # stsynphot table

        colname = 'THROUGHPUT'
        thru = hdu.data.field(colname)
        thru[:] = self.interpTdsFactors(row, self.thru_mjd, wavelength)

        useafter_mjd = round(self.useafter_mjd)
        for i in range(len(time)):
            colname = f'MJD#{time[i]:.1f}'
            thru = hdu.data.field(colname)
            if abs(time[i] - useafter_mjd) < 0.5:
                # Throughput values at USEAFTER date.
                thru[:] = 1.
            else:
                thru[:] = self.interpTdsFactors(row, time[i], wavelength)

    def createFile(self, hdu, row, overwrite=False):
        """Create the ``stsynphot`` file.

        This creates a pyfits HDUList object using the primary header/data
        unit from the input TDS table, appends ``hdu``, then writes the
        HDUList to the output file.

        Parameters
        ----------
        hdu : `astropy.io.fits.BinTableHDU`
            Header/data unit for the ``stsynphot`` table.

        row : int
            Current row number (0-indexed) in input TDS table. This is
            used only for writing a history record in the output header.
        """
        primary_hdu = fits.PrimaryHDU(header=self.phdr_tds)
        ofd = fits.HDUList(primary_hdu)
        phdr = ofd[0].header
        # these keywords may have been copied from the input TDS file
        del phdr['origin']
        del phdr['iraf-tlm']
        del phdr['date']
        del phdr['obstype']
        del phdr['aperture']
        del phdr['cenwave']
        del phdr['filetype']
        del phdr['vcalcos']
        del phdr['coscoord']
        del phdr['detector']
        del phdr['comment']
        del phdr['history']
        phdr.update('filename', os.path.basename(self.output))
        if len(self.orig_input) + len(self.orig_outroot) < 48:
            phdr.add_history(f'tdspysyn:  input = {self.orig_input}, '
                             f'outroot = {self.orig_outroot}')
        else:
            phdr.add_history(f'tdspysyn:  input = {self.orig_input}')
            phdr.add_history(f'tdspysyn:  outroot = {self.orig_outroot}')
        phdr.add_history(f'tdspysyn:  min_wl = {self.min_wl:.8g}, '
                         f'max_wl = {self.max_wl:.8g}, dwl = {self.dwl:.5g}')
        if self.last_mjd is None:
            last_mjd = 'None'
        else:
            last_mjd = str(self.last_mjd)
        # convert the row number to one-indexed
        phdr.add_history(f'tdspysyn:  thru_mjd = {self.thru_mjd:.10g}, '
                         f'last_mjd = {last_mjd}, row = {row + 1}')

        ofd.append(hdu)
        ofd.writeto(self.output, overwrite=overwrite)

    def useafterMjd(self):
        """Get the USEAFTER date and convert it to MJD.

        Returns
        -------
        mjd : float
            Modified Julian Date.

        """
        if self.useafter == 'missing':
            raise RuntimeError('Keyword USEAFTER is missing '
                               'from the primary header.')

        # The format is expected to be 'month day year hh:mm:ss',
        # for example, "May 11 2009 00:00:00".
        words = self.useafter.split()
        len_words = len(words)
        if len_words < 3 or len_words > 6:
            raise RuntimeError(
                f"Can't interpret USEAFTER date '{self.useafter}'")
        if words[1].endswith(','):     # allow "May 11, 2009"
            words[1] = words[1][:-1]
        if words[2].endswith(','):     # allow "May 11 2009, 00:00:00"
            words[2] = words[2][:-1]
        month_str = words[0][0:3].lower()
        months = ['jan', 'feb', 'mar', 'apr', 'may', 'jun',
                  'jul', 'aug', 'sep', 'oct', 'nov', 'dec']
        if month_str not in months:
            raise RuntimeError(
                "Don't understand the month in USEAFTER date "
                f"'{self.useafter}'")
        month = months.index(month_str) + 1  # 1-indexed month
        day = int(words[1])
        year = int(words[2])
        if len_words == 4:
            hms = words[3].split(':')
            len_hms = len(hms)
            if len_hms < 1 or len_hms > 3:
                raise RuntimeError(
                    f"Can't interpret USEAFTER date '{self.useafter}'")
            if len_hms == 1:
                fraction = float(hms[0]) / 24.
            elif len_hms == 2:
                fraction = (float(hms[0]) + float(hms[1]) / 60.) / 24.
            elif len_hms == 3:
                fraction = ((float(hms[0]) +
                            (float(hms[1]) +
                             float(hms[2]) / 60.) / 60.) / 24.)
        elif len_words == 5:
            hours = float(words[3])
            minutes = float(words[4])
            fraction = (hours + minutes / 60.) / 24.
        elif len_words == 6:
            hours = float(words[3])
            minutes = float(words[4])
            seconds = float(words[5])
            fraction = (hours + (minutes + seconds / 60.) / 60.) / 24.
        else:
            fraction = 0.

        day = float(day) + fraction
        useafter_mjd = toMjd(year, month, day)

        return useafter_mjd


class CosTds(ConvertTds):
    """Class to handle COS TDS."""

    def __init__(self, *args):
        super().__init__(*args)

        if self.ref_time == 'missing':
            raise RuntimeError('Keyword REF_TIME is missing '
                               'from the table header.')

    def makeFilename(self, row):
        """Create a name for the output COS ``stsynphot`` file.

        The file name will be the root that was specified by the user,
        then the values of OPT_ELEM, SEGMENT, and APERTURE from ``row`` in
        the input TDS table (each converted to lower case and preceded
        by ``_``), and ending with ``.fits``.

        Parameters
        ----------
        row : int
            Zero-indexed row number in the input TDS table.

        """
        opt_elem = self.data.field('opt_elem')[row].lower()
        segment = self.data.field('segment')[row].lower()
        aperture = self.data.field('aperture')[row].lower()

        self.output = f'{self.outroot}_{opt_elem}_{segment}_{aperture}.fits'

    def interpTdsFactors(self, row, mjd, wavelength):
        """Read and interpret the specified row of the input COS TDS table.

        Parameters
        ----------
        row : int
            Current row number (0-indexed) in input TDS table.

        mjd : float
            Modified Julian Date.

        wavelength : array_like
            Array of wavelengths for the ``stsynphot`` table (float64).

        Returns
        -------
        factor : array_like
            Array of time-dependent sensitivity factors for the current
            row, interpreted from the input TDS table and interpolated
            to correspond to the values in ``wavelength``.
        """
        from calcos import ccos

        nwl = self.data.field('nwl')[row]
        nt = self.data.field('nt')[row]
        wl_tds = self.data.field('wavelength')[row]    # 1-D array
        time_tds = self.data.field('time')[row]        # 1-D array
        # slope_pct is percent per year
        slope_pct = self.data.field('slope')[row]      # 2-D array
        intercept = self.data.field('intercept')[row]  # 2-D array

        # This section is needed because pyfits currently ignores TDIMi.
        maxt = len(time_tds)
        maxwl = len(wl_tds)
        slope_pct = np.reshape(slope_pct, (maxt, maxwl))
        intercept = np.reshape(intercept, (maxt, maxwl))

        # Find the time interval that includes the time of observation.
        if nt == 1 or mjd >= time_tds[nt-1]:
            i = nt - 1
        else:
            for i in range(nt-1):
                if mjd < time_tds[i+1]:
                    break

        # The slope in the tdstab is in percent per year.  Convert the time
        # interval to years, and convert the slope to fraction per year.
        # If the time of observation is before the first time in the table
        # or after the last time, the correction factor is to be the factor
        # at the first time or the last time respectively.  This is done by
        # setting delta_t to be the difference from the reference time to
        # the first or last time.
        if mjd < time_tds[0]:
            delta_t = (time_tds[0] - self.ref_time) / DAYS_PER_YEAR
        elif mjd > time_tds[nt-1]:
            delta_t = (time_tds[nt-1] - self.ref_time) / DAYS_PER_YEAR
        else:
            delta_t = (mjd - self.ref_time) / DAYS_PER_YEAR
        slope = slope_pct / 100.        # slope as a fraction per year

        # Take the slice [0:nwl] to avoid using elements that may not be valid,
        # and because the array of factors should be the same length as the
        # set of wavelengths that have been specified.
        wl_tds = wl_tds[0:nwl]
        factor_tds = delta_t * slope[i][0:nwl] + intercept[i][0:nwl]
        factor = np.zeros(len(wavelength), dtype=np.float64)

        # Interpolate factor_tds at each wavelength.
        ccos.interp1d(wl_tds, factor_tds, wavelength, factor)

        return factor


class StisTds(ConvertTds):
    """Class to handle STIS TDS."""

    def makeFilename(self, row):
        """Create a name for the output STIS ``stsynphot`` file.

        The file name will be the root that was specified by the user, then
        ``_`` and the value of OPT_ELEM (converted to lower case) from ``row``
        in the input TDS table, and ending with ``.fits``.

        Parameters
        ----------
        row : int
            Zero-indexed row number in the input TDS table.

        """
        opt_elem = self.data.field('opt_elem')[row].lower()

        self.output = f'{self.outroot}_{opt_elem}.fits'

    def interpTdsFactors(self, row, mjd, wavelength):
        """Read and interpret the specified row of the input STIS TDS table.

        Parameters
        ----------
        row : int
            Current row number (0-indexed) in input TDS table.

        mjd : float
            Modified Julian Date.

        wavelength : array_like
            Array of wavelengths for the ``stsynphot`` table (float64).

        Returns
        -------
        factor : array_like
            Array of time-dependent sensitivity factors for the current
            row, interpreted from the input TDS table and interpolated
            to correspond to the values in ``wavelength``.

        """
        from calcos import ccos

        nwl = self.data.field('nwl')[row]
        nt = self.data.field('nt')[row]
        wl_tds = self.data.field('wavelength')[row]    # 1-D array
        time_tds = self.data.field('time')[row]        # 1-D array
        # slope_pct is percent per year
        slope_pct = self.data.field('slope')[row]      # 2-D array

        # This section is needed because pyfits currently ignores TDIMi.
        maxt = len(time_tds)
        maxwl = len(wl_tds)
        slope_pct = np.reshape(slope_pct, (maxt, maxwl))

        # Convert the slope from percent to a fraction,
        # and convert from per year to per day.
        slope = slope_pct / (DAYS_PER_YEAR * 100.)

        factor_tds = np.ones(nwl, dtype=np.float64)
        for i in range(nwl):
            for j in range(nt):
                if j == nt-1 or mjd <= time_tds[j+1]:
                    factor_tds[i] += (mjd - time_tds[j]) * slope[j, i]
                    break
                else:
                    factor_tds[i] += ((time_tds[j + 1] - time_tds[j]) *
                                      slope[j, i])

        # Take the slice [0:nwl] to avoid using elements that may not be valid,
        # and because the array of factors should be the same length as the
        # set of wavelengths that have been specified.
        wl_tds = wl_tds[0:nwl]
        factor = np.zeros(len(wavelength), dtype=np.float64)

        # Interpolate factor_tds at each wavelength.
        ccos.interp1d(wl_tds, factor_tds, wavelength, factor)

        return factor


if __name__ == '__main__':
    main()
