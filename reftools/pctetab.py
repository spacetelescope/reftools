"""Functions for ACS PCTETAB reference file.

Examples
--------
>>> from reftools import pctetab
>>> pctetab.MakePCTETab(
...     'pctetab_pcte.fits', 'pctetab_dtdel.txt',
...     ['pctetab_chgleak-1.txt', 'pctetab_chgleak-2.txt'],
...     'pctetab_levels.txt', 'pctetab_scaling.txt',
...     'pctetab_column_scaling.txt', history_file='pctetab_history.txt')

"""
# :Authors: Pey Lian Lim, Matt Davis
# :Organization: Space Telescope Science Institute
# :History:
#    * 2010-08-31 PLL created this module.
#    * 2010-11-09 PLL added RN2_NIT keyword and updated documentation.
#    * 2011-04-25 MRD updated for new CTE algorithm parameters
#    * 2011-07-18 MRD updated to handle time dependence
#    * 2011-11-29 MRD updated with column-by-column CTE scaling
#    * 2013-08-13 PLL removed depreciated PyFITS calls and cleaned up codes.
from __future__ import absolute_import, division, print_function

# STDLIB
import glob
import sys
import os

# THIRD-PARTY
import numpy as np
from astropy.io import fits

__version__ = '1.2.0'
__vdata__ = '13-Aug-2013'


class PCTEFileError(Exception):
    """Generic exception for errors in this module."""
    pass


class _Text2Fits(object):
    """Helper class for making the CTE parameters file (PCTETAB) from a
    collection of data saved in text files. The resulting fits file will
    have information in the primary header and in at least four table extensions.

    """
    def __init__(self):
        self.header = None
        self.dtde = None
        self.charge_leak = []
        self.levels = None
        self.scale = None
        self.col_scale = None
        self.out_name = None

    def make_header(self, out_name, sim_nit, shft_nit, read_noise, noise_model,
                    oversub_thresh, nchg_leak, useafter, pedigree, creatorName,
                    history_file, detector):
        """Make the primary extension header for the pctetab.

        Parameters
        ----------
        out_name : str
            Name of pcte fits file being created. May include path.

        sim_nit : int
            Value for ``SIM_NIT`` keyword in PCTEFILE header.
            Number of iterations of readout simulation per column.

        shft_nit : int
            Value for ``SHFT_NIT`` keyword in PCTEFILE header.
            Number of shifts each readout simulation is broken up into.
            A large number means pixels are shifted a smaller number of rows
            before the CTE is evaluated again.

        read_noise : float
            Value for ``RN_CLIP`` keyword in PCTEFILE
            EXT 0. This is the maximum amplitude of read noise
            used in the read noise mitigation. Unit is in electrons.

        noise_model : {0, 1, 2}
            Select the method to be used for readnoise removal::

                0: no read noise smoothing
                1: standard smoothing
                2: strong smoothing

        oversub_thresh : float
            Value for ``SUBTHRSH`` keyword in PCTEFILE header. CTE corrected
            pixels taken below this value are re-corrected. Unit is in electrons.

        useafter : str
            Value for ``USEAFTER`` keyword.

        pedigree : str
            Value for ``PEDIGREE`` keyword.

        creatorName : str
            Name of the person generating this `fitsFile`.

        historyFile : str
            ASCII file containing ``HISTORY`` lines for
            EXT 0. Include path. Each row will produce one
            ``HISTORY`` line.

        detector : str
            Supported detector.

        """
        self.out_name = out_name

        self.header = fits.PrimaryHDU()

        self.header.header['ORIGIN'] = 'STScI-STSDAS/TABLES'
        self.header.header['FILENAME'] = os.path.basename(out_name)
        self.header.header['FILETYPE'] = 'PIXCTE'
        self.header.header['TELESCOP'] = 'HST'
        self.header.header['USEAFTER'] = useafter
        self.header.header['PEDIGREE'] = pedigree
        self.header.header['DESCRIP'] = 'Parameters needed for pixel-based CTE correction ------------------'  # Must have 67 char
        self.header.header.add_comment('= \'Created or updated by {0:s}\''.format(creatorName), before='ORIGIN')
        self.header.header['NCHGLEAK'] = (nchg_leak, 'number of chg_leak extensions')

        # This is detector-specific
        if detector == 'WFC':
            self.header.header['INSTRUME'] = 'ACS'
            self.header.header['DETECTOR'] = 'WFC'
        else:
            raise PCTEFileError('Detector not supported: {0:s}'.format(
                str(detector)))

        # Optional HISTORY
        if os.path.isfile(history_file):
            with open(history_file) as fin:
                for line in fin:
                    self.header.header.add_history(line[:-1])

        # the number of readout simulations done per column
        self.header.header['SIM_NIT'] = (int(sim_nit), 'number of readout simulations done per column')

        # the number of shifts each column readout simulation is broken up into
        self.header.header['SHFT_NIT'] = (int(shft_nit), 'num shifts col readout sim is broken up into')

        # read noise level
        self.header.header['RN_CLIP'] = (float(read_noise), 'Read noise level in electrons.')

        # read noise smoothing algorithm
        self.header.header['NSEMODEL'] = (int(noise_model), 'Read noise smoothing algorithm.')

        # over-subtraction correction threshold
        self.header.header['SUBTHRSH'] = (float(oversub_thresh), 'Over-subtraction correction threshold.')

    def make_dtde(self,dtde_file):
        """Make fits extension containing the dtde data that describes the
        marginal loss to CTE at the charge levels given in the file.

        The input file should have two columns with headers DTDE and Q,
        in that order.
        The first column is the dtde data and the second is the corresponding
        charge level for that dtde value.

        The file should have format::

            DTDE  Q
            float int
            ...   ...

        Lines beginning with # are ignored.

        Parameters
        ----------
        dtde_file : str
            Path to text file containing dtde data.

        """
        if not os.path.isfile(dtde_file):
            raise IOError('Invalid dtde file: {0:s}'.format(str(dtde_file)))

        lRange, colName, colData, colForm, colUnit = 0, {}, {}, {}, {}

        # read in dtde data from text file
        with open(dtde_file) as fin:
            for line in fin:
                # skip comments
                if line.startswith('#'):
                    continue

                row = line.split()

                # column names
                if row[0] == 'DTDE':
                    colRange = range(len(row))
                    for i in colRange:
                        colName[i] = row[i]
                        colData[i] = []
                # data
                else:
                    for i in colRange:
                        colData[i].append(row[i])

        # convert data to numpy arrays
        colData[0] = np.array(colData[0], dtype=np.float32)
        colForm[0] = 'E'
        colUnit[0] = ''

        colData[1] = np.array(colData[1], dtype=np.int32)
        colForm[1] = 'J'
        colUnit[1] = 'DN/S'

        c0 = fits.Column(name=colName[0], format=colForm[0], array=colData[0])
        c1 = fits.Column(name=colName[1], format=colForm[1], unit=colUnit[1], array=colData[1])

        self.dtde = fits.BinTableHDU.from_columns(fits.ColDefs([c0, c1]))
        self.dtde.header['EXTNAME'] = 'DTDE'
        self.dtde.header['DATAFILE'] = (os.path.basename(dtde_file), 'data source file')

    def make_charge_leak(self, chg_leak_file, num):
        """Make fits extension containing parameterization of CTE losses along
        the CTE tail and across different charge levels.

        The input file should contain 5 columns with following format::

            NODE LOG_Q_1 LOG_Q_2 LOG_Q_3 LOG_Q_4
            int  float   float   float   float
            ...  ...     ...     ...     ...

        Lines beginning with # are ignored.

        Parameters
        ----------
        chg_leak_file : str
            Path to text file containing charge leak data.

        num : int
            Number to append to extension name since there may be more than
            one charge_leak extension.

        """
        if not os.path.isfile(chg_leak_file):
            raise IOError('Invalid charge leak file: {0:s}'.format(
                str(chg_leak_file)))

        colRange, colName, colData, colForm, colUnit = 0, {}, {}, {}, {}

        mjd1 = None
        mjd2 = None

        # read in dtde data from text file
        with open(chg_leak_file) as fin:
            for line in fin:
                # skip comments
                if line.startswith('#'):
                    continue

                row = line.split()

                # MJD parameters
                if row[0] == 'MJD1':
                    mjd1 = float(row[1])
                elif row[0] == 'MJD2':
                    mjd2 = float(row[1])

                # column names
                elif row[0] == 'NODE':
                    colRange = range(len(row))
                    for i in colRange:
                        colName[i] = row[i]
                        colData[i] = []
                # data
                else:
                    for i in colRange:
                        colData[i].append(row[i])

        # make sure we got our MJD values
        if not mjd1:
            raise PCTEFileError('MJD1 parameter not correctly specified in {0:s}'.format(chg_leak_file))
        elif not mjd2:
            raise PCTEFileError('MJD2 parameter not correctly specified in {0:s}'.format(chg_leak_file))

        # Convert data to Numpy arrays
        colData[0] = np.array(colData[0], dtype=np.int16)
        colForm[0] = 'I'
        colUnit[0] = 'PIXEL'

        for i in colRange[1:]:
            colData[i] = np.array(colData[i], dtype=np.float32)
            colForm[i] = 'E'
            colUnit[i] = 'FRACTION'

        # Write to FITS table extension
        tabData = [fits.Column(name=colName[i], format=colForm[i], unit=colUnit[i], array=colData[i]) for i in colRange]

        self.charge_leak.append(fits.BinTableHDU.from_columns(fits.ColDefs(tabData)))
        self.charge_leak[-1].header['EXTNAME'] = 'CHG_LEAK{0:d}'.format(num)
        self.charge_leak[-1].header['MJD1'] = (mjd1, 'start valid time range for data')
        self.charge_leak[-1].header['MJD2'] = (mjd2, 'end valid time range for data')
        self.charge_leak[-1].header['DATAFILE'] = (os.path.basename(chg_leak_file), 'data source file')

    def make_levels(self, levels_file):
        """Make fits extension containing charge levels at which to evaluate CTE
        losses (as opposed to every level from 0 - 99999).

        The input file should have a single column with the following format::

            LEVEL
            int
            ...

        Columns beginning with # are ignored.

        Parameters
        ----------
        levels_file : str
            Text file containing charge levels at which to do CTE evaluation.

        """
        if not os.path.isfile(levels_file):
            raise IOError('Invalid levels file: {0:s}'.format(
                str(levels_file)))

        colData = []

        # read in data from text file
        with open(levels_file) as fin:
            for line in fin:
                # skip comments
                if line.startswith('#'):
                    continue

                row = line.split()

                # column heading
                if row[0] == 'LEVEL':
                    colName = row[0]
                else:
                    colData.append(row[0])

        colData = np.array(colData, dtype=np.int32)
        colForm = 'J'

        c1 = fits.Column(name=colName, format=colForm, array=colData)

        self.levels = fits.BinTableHDU.from_columns(fits.ColDefs([c1]))
        self.levels.header['EXTNAME'] = 'LEVELS'
        self.levels.header['DATAFILE'] = (os.path.basename(levels_file), 'data source file')

    def make_scale(self, scale_file):
        """Make fits extension containing time dependent CTE scaling.

        The input file should have two columns with the following format::

            MJD     SCALE
            float   float
            ...     ...

        Columns beginning with # are ignored.

        Parameters
        ----------
        scale_file : str
            Text file containing time dependent CTE scaling parameters.

        """
        if not os.path.isfile(scale_file):
            raise IOError('Invalid scale file: {0:s}'.format(
                str(scale_file)))

        lRange, colName, colData, colForm, colUnit = 0, {}, {}, {}, {}

        # read in dtde data from text file
        with open(scale_file) as fin:
            for line in fin:
                # skip comments
                if line.startswith('#'):
                    continue

                row = line.split()

                # column names
                if row[0] == 'MJD':
                    colRange = range(len(row))
                    for i in colRange:
                        colName[i] = row[i]
                        colData[i] = []
                # data
                else:
                    for i in colRange:
                        colData[i].append(row[i])

        # convert data to numpy arrays
        colData[0] = np.array(colData[0], dtype=np.float32)
        colForm[0] = 'E'
        colUnit[0] = 'DAYS'

        colData[1] = np.array(colData[1], dtype=np.float32)
        colForm[1] = 'E'
        colUnit[1] = 'FRACTION'

        c0 = fits.Column(name=colName[0], format=colForm[0], unit=colUnit[0], array=colData[0])
        c1 = fits.Column(name=colName[1], format=colForm[1], unit=colUnit[1], array=colData[1])

        self.scale = fits.BinTableHDU.from_columns(fits.ColDefs([c0,c1]))
        self.scale.header['EXTNAME'] = 'CTE_SCALE'
        self.scale.header['DATAFILE'] = (os.path.basename(scale_file), 'data source file')

    def make_column_scale(self, column_file):
        """Make fits extension containing column by column CTE scaling.

        The input file should have 5 columns with the following format::

            COLUMN  AMPA    AMPB    AMPC    AMPD
            int     float   float   float   float
            ...     ...     ...     ...     ...

        Lines beginning with # are ignored.

        Parameters
        ----------
        column_file : str
            Text file containing CTE column-by-column scaling.

        """
        if not os.path.isfile(column_file):
            raise IOError('Invalid column scale file: {0:s}'.format(
                str(column_file)))

        lRange, colName, colData, colForm, colUnit = 0, {}, {}, {}, {}

        # read in dtde data from text file
        with open(column_file) as fin:
            for line in fin:
                # skip comments
                if line.startswith('#'):
                    continue

                row = line.split()

                # column names
                if row[0] == 'COLUMN':
                    colRange = range(len(row))
                    for i in colRange:
                        colName[i] = row[i]
                        colData[i] = []
                # data
                else:
                    for i in colRange:
                        colData[i].append(row[i])

        # convert data to numpy arrays
        colData[0] = np.array(colData[0], dtype=np.int32)
        colForm[0] = 'J'
        colUnit[0] = 'COLUMN NUMBER'

        colData[1] = np.array(colData[1], dtype=np.float32)
        colForm[1] = 'E'
        colUnit[1] = 'FRACTION'

        colData[2] = np.array(colData[2], dtype=np.float32)
        colForm[2] = 'E'
        colUnit[2] = 'FRACTION'

        colData[3] = np.array(colData[3], dtype=np.float32)
        colForm[3] = 'E'
        colUnit[3] = 'FRACTION'

        colData[4] = np.array(colData[4], dtype=np.float32)
        colForm[4] = 'E'
        colUnit[4] = 'FRACTION'

        c0 = fits.Column(name=colName[0], format=colForm[0], unit=colUnit[0], array=colData[0])
        c1 = fits.Column(name=colName[1], format=colForm[1], unit=colUnit[1], array=colData[1])
        c2 = fits.Column(name=colName[2], format=colForm[2], unit=colUnit[2], array=colData[2])
        c3 = fits.Column(name=colName[3], format=colForm[3], unit=colUnit[3], array=colData[3])
        c4 = fits.Column(name=colName[4], format=colForm[4], unit=colUnit[4], array=colData[4])

        self.col_scale = fits.BinTableHDU.from_columns(fits.ColDefs([c0,c1,c2,c3,c4]))
        self.col_scale.header['EXTNAME'] = 'COL_SCALE'
        self.col_scale.header['DATAFILE'] = (os.path.basename(column_file), 'data source file')

    def make_fits(self):
        """Combine primary and table extensions into an HDU List and
        save to fits file.

        The methods make_header, make_dtde, make_charge, and make_levels
        must have been succesfully run before calling this method.

        Raises
        ------
        PCTEFileError
            If any of the necessary extensions have not been made.

        """
        if not self.header:
            raise PCTEFileError('Fits header has not been prepared: '
                                'call make_header method first.')
        if not self.dtde:
            raise PCTEFileError('DTDE extension has not been prepared: '
                                'call make_dtde method first.')
        if not self.charge_leak:
            raise PCTEFileError('Charge leak extension has not been prepared: '
                                'call make_charge_leak method first.')
        if not self.levels:
            raise PCTEFileError('Levels extension has not been prepared: '
                                'call make_levels method first.')
        if not self.scale:
            raise PCTEFileError('Scale extension has not been prepared: '
                                'call make_scale method first.')
        if not self.col_scale:
            raise PCTEFileError('Column scaline extension has not been '
                                'prepared: call make_column_scale method first')

        hduList = fits.HDUList([self.header, self.dtde, self.levels, self.scale,
                                self.col_scale] + self.charge_leak)

        hduList.writeto(self.out_name, clobber=True)


def MakePCTETab(out_name, dtde_file, chg_leak_file, levels_file, scale_file,
                column_file, sim_nit=7, shft_nit=7, read_noise=5.0,
                noise_model=1, oversub_thresh=-10,
                useafter='Mar 01 2002 00:00:00',
                pedigree='INFLIGHT 01/03/2002 22/07/2010',
                creatorName='ACS Team', history_file='', detector='WFC'):
    """Make the CTE parameters reference file.

    Parameters
    ----------
    out_name : str
        Name of pcte fits file being created. May include path.

    dtde_file : str
        Path to text file containing dtde data.

        The file should have 2 columns with the following format::

            DTDE  Q
            float int
            ...   ..

        Lines beginning with # are ignored.

    chg_leak_file : str or list of str
        Path to text file(s) containing charge leak data. If passed as a string
        the string may contain wild cards so that multiple files are specified.

        The input file should contain 5 columns with following format::

            NODE LOG_Q_1 LOG_Q_2 LOG_Q_3 LOG_Q_4
            int  float   float   float   float
            ...  ...     ...     ...     ...

        Lines beginning with # are ignored.

    levels_file : str
        Text file containing charge levels at which to do CTE evaluation.

        The input file should have a single column with the following format::

            LEVELS
            int
            ...

        Lines beginning with # are ignored.

    scale_file : str
        Text file containing CTE scaling parameters

        The input file should have two columns with the following format::

            MJD     SCALE
            float   float
            ...     ...

        Lines beginning with # are ignored.

    column_file : str
        Text file containing CTE column-by-column scaling.

        The input file should have 5 columns with the following format::

            COLUMN  AMPA    AMPB    AMPC    AMPD
            int     float   float   float   float
            ...     ...     ...     ...     ...

        Lines beginning with # are ignored.

    sim_nit : int, optional
        Number of iterations of readout simulation per column.

    shft_nit : int, optional
        Number of shifts each readout simulation is broken up into.
        A large number means pixels are shifted a smaller number of rows
        before the CTE is evaluated again.

    read_noise : float
        Value for ``RN_CLIP`` keyword in PCTEFILE
        EXT 0. This is the maximum amplitude of read noise
        used in the read noise mitigation. Unit is in electrons.

    noise_model : {0, 1, 2}
        Select the method to be used for readnoise removal::

            0: no read noise smoothing
            1: standard smoothing
            2: strong smoothing

    oversub_thresh : float
        Value for ``SUBTHRSH`` keyword in PCTEFILE header. CTE corrected
        pixels taken below this value are re-corrected. Unit is in electrons.

    useafter : str, optional
        Value for ``USEAFTER`` keyword.
        Defaults to 'Mar 01 2002 00:00:00'

    pedigree : str, optional
        Value for ``PEDIGREE`` keyword.
        Defaults to 'INFLIGHT 01/03/2002 22/07/2010'

    creatorName : str, optional
        Name of the person generating this `fitsFile`.
        Defaults to 'ACS Team'

    historyFile : str, optional
        ASCII file containing ``HISTORY`` lines for
        EXT 0. Include path. Each row will produce one
        ``HISTORY`` line.
        Defaults to ''

    detector : str, optional
        Supported detector. Defaults to 'WFC'

    Examples
    --------
    Saving file pctetab_pcte.fits with the command:

    >>> MakePCTETab(
    ...     'pctetab_pcte.fits', 'pctetab_dtdel.txt',
    ...     ['pctetab_chgleak-1.txt', 'pctetab_chgleak-2.txt'],
    ...     'pctetab_levels.txt', 'pctetab_scaling.txt',
    ...     'pctetab_column_scaling.txt', history_file='pctetab_history.txt')

    """
    # give the output file it's official suffix
    if out_name.find('_pcte.fits') == -1:
        out_name = out_name + '_pcte.fits'

    # test for the presence of the input files
    if not os.path.isfile(dtde_file):
        raise IOError('Invalid dtde file: {0:s}'.format(str(dtde_file)))

    if isinstance(chg_leak_file, str):
        chg_leak_file = glob.glob(chg_leak_file)

    for f in chg_leak_file:
        if not os.path.isfile(f):
            raise IOError('Invalid charge leak file: {0:s}'.format(
                str(chg_leak_file)))

    nchg_leak = len(chg_leak_file)

    if not os.path.isfile(levels_file):
        raise IOError('Invalid levels file: {0:s}'.format(str(levels_file)))

    if not os.path.isfile(scale_file):
        raise IOError('Invalid scale file: {0:s}'.format(str(scale_file)))

    if not os.path.isfile(column_file):
        raise IOError('Invalid column scaling file: {0:s}'.format(
            str(column_file)))

    # make Text2Fits object and run it's methods to construct fits extensions
    t2f = _Text2Fits()
    t2f.make_header(out_name, sim_nit, shft_nit, read_noise, noise_model,
                    oversub_thresh, nchg_leak, useafter, pedigree,
                    creatorName, history_file, detector)
    t2f.make_dtde(dtde_file)

    for i, f in enumerate(chg_leak_file, 1):
        t2f.make_charge_leak(f, i)

    t2f.make_levels(levels_file)

    t2f.make_scale(scale_file)

    t2f.make_column_scale(column_file)

    # have t2f save the fits file
    print('Saving file {0:s}'.format(str(out_name)))
    t2f.make_fits()
