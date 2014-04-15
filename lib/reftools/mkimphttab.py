"""Use this module to create an IMPHTTAB for an instrument."""
from __future__ import print_function, division

# STDLIB
import os
import re
import sys
import tempfile
import time as _time
from collections import OrderedDict

# THIRD-PARTY
import astropy
import numpy as np
from astropy import log
from astropy.io import fits

# STSCI
import pysynphot as S
from pysynphot import refs

# LOCAL
from . import graphfile as sgf

__all__ = ['create_table', 'create_nicmos_table', 'create_table_from_table']
__version__ = '0.5'
__vdate__ = '05-Nov-2013'


def compute_values(obsmode, component_dict):
    """Compute the 3 basic photometric values needed for a given
    observation mode string using `pysynphot`.

    Values calculated:

        * PHOTFLAM - Unit response in FLAM.
        * PHOTPLAM - Pivot wavelength.
        * PHOTBW - Bandwidth.

    Parameters
    ----------
    obsmode : str
        Observation mode for which to calculate values.

    component_dict : dict
        A dictionary in which to cache opened component objects.
        May be empty.

    Returns
    -------
    valdict : dict
        Dictionary with photometry keywords as keys.

    """
    # Define the bandpass for this obsmode
    bp = S.ObsBandpass(obsmode, component_dict=component_dict)

    # compute the photometric values
    uresp = bp.unit_response()
    if not np.isfinite(uresp):
        uresp = 0.0

    return {'PHOTFLAM': uresp, 'PHOTPLAM': bp.pivot(), 'PHOTBW': bp.photbw()}


def compute_synphot_values(obsmode):
    """Calculate the same values as :func:`compute_values` but
    using IRAF SYNPHOT.

    .. note::

        This is replaced by :func:`compute_values` but kept for debugging.

    Parameters
    ----------
    obsmode : str
        Observation mode for which to calculate values.

    Returns
    -------
    valdict : dict
        Dictionary with photometry keywords as keys.

    """
    from pyraf import iraf
    from iraf import stsdas, hst_calib, synphot

    synphot.bandpar(obsmode, output=tmpfits, Stdout=1)
    tmpfits = os.path.join(tempfile.gettempdir(), 'temp.fits')

    with fits.open(tmpfits) as f:
        d = f[1].data
        photflam = d['uresp'][0]
        pivot = d['pivwv'][0]
        rmswidth = d['bandw'][0]

    os.remove(tmpfits)

    return {'PHOTFLAM': photflam, 'PHOTPLAM': pivot, 'PHOTBW': rmswidth}


def expand_obsmodes(basemode, pardict):
    """Generate a set of observation mode strings spanning all the
    combinations of the ``basemode`` as specified in ``pardict``.

    Parameters
    ----------
    basemode : str
        Root of the observation mode string including all
        non-parameterized components. Example::

            acs,wfc1,f850lp

    pardict : dict
        Dictionary with the name of each parameterized variable
        as the key and a list of values for that key. For example::

            {'mjd': [52334, 53919.99, 53920, 55516],
             'fr853n': [8158.0, 8531.5, 8905.0]}

    Returns
    -------
    olist : list of str
        Expanded observation modes.

    """
    basemode = basemode.lower()
    obsmode_str = '%s,%s#%0.4f'
    olist = list()

    if len(pardict) > 0:
        for k in pardict.keys():
            basemode = basemode.replace('{0},'.format(k), '')

    # we don't have any parameterized variables, so just return the basemode
    if len(pardict) == 0:
        olist.append(basemode.rstrip(','))

    # Build up list of OBSMODE covering all combinations of parameterized
    # variable values
    elif len(pardict) == 1:
        key = pardict.keys()[0]
        for val in pardict[key]:
            ostr = basemode.replace(key.lower(), key.lower() + str(val))
            olist.append(ostr)

    else:
        nkeys = len(pardict)
        for nkey in range(nkeys - 1):
            key = pardict.keys()[nkey]
            for val in pardict[key]:
                pdict = {}
                for k in pardict.keys()[nkey+1:]:
                    pdict[k] = pardict[k]
                ostr = basemode.replace(key.lower(), key.lower()+str(val))
                olist.extend(expand_obsmodes(ostr, pdict))

    return olist


def interpret_obsmode(obsmode):
    """Convert a full observation mode string with parameterized values
    into a string which only lists the parameterized variable names
    without values for comparison with the observation mode strings
    in the IMPHTTAB table.

    Parameters
    ----------
    obsmode : str
        Full observation mode string.

    Returns
    -------
    omode : str
         The OBSMODE value for each row in IMPHTTAB.

    Examples
    --------
    >>> interpret_obsmode('acs,wfc1,mjd#52334.0000,fr853n#8158.0000')
    'acs,wfc1,mjd#,fr853n#'

    """
    ospl = obsmode.split(',')
    omode = ''
    for o in ospl:
        if '#' in o:
            o = o.split('#')[0] + '#'
        omode += o + ','
    omode = omode.rstrip(',')
    return omode


def parse_filters(filters):
    """Parse the filters specification.

    Parameters
    ----------
    filters : str
        Filter specification.

    Returns
    -------
    fval : str
        Non-parameterized filter names, separated by commas.

    fpars : list of str
        A list of parameterized filter names.

    """
    fspl = filters.split(',')
    fpars = list()
    fval = ''
    for f in fspl:
        if f.find('#') < 0:
            fval += f + ','
        else:
            fpars.append(f)
    fval = fval[:-1]
    return fval, fpars


def get_date():
    """Returns a formatted string with the current date and time.

    .. note:: Copied from :func:`stsci.tools.fileutil.get_date`.

    Returns
    -------
    date_str : str
        Date and time in the format of YYYY-MM-DDTHH:MM:SS.

    """
    return _time.strftime('%Y-%m-%dT%H:%M:%S', _time.localtime(_time.time()))


def make_pri_hdu(filename, numpars, instrument, detector, pedigree, useafter):
    """Create a primary header for the multi-extension FITS reference table.

    Parameters
    ----------
    filename : str
        Output filename that the header belongs to.

    numpars : int
        Number of parameterized variables.

    instrument : str
        Instrument name.

    detector : str
        Detector name.

    pedigree : str
        Data pedigree.

    useafter : str
        Useafter date in the format of ``MMM DD YYYY HH:MM:SS``.

    """
    d = refs.getref()
    phdu = fits.PrimaryHDU()

    phdu.header['DATE'] = (get_date(), 'Date FITS file was generated')
    phdu.header['ORIGIN'] = ('astropy-fits',
                             'astropy version {0}'.format(astropy.__version__))
    phdu.header['FILENAME'] = (filename, 'Name of file')
    phdu.header['FILETYPE'] = ('Image photometry table', 'File type')
    phdu.header['NEXTEND'] = (3, 'Number of extensions in file')
    phdu.header['PHOTZPT'] = (-21.1, 'Photometric zero-point for STMAG system')
    phdu.header['PARNUM'] = (numpars, 'Number of parameterized variables')
    phdu.header['DBTABLE'] = ('IMPHTTAB', 'Database table')
    phdu.header['INSTRUME'] = (instrument, 'Instrument name')
    phdu.header['DETECTOR'] = (detector, 'Detector name')
    phdu.header['SYNSWVER'] = (S.__version__,
                               'Version of synthetic photometry software')
    phdu.header['MKTABVER'] = (__version__, 'Version of reftools.mkimphttab')
    phdu.header['GRAPHTAB'] = (os.path.basename(d['graphtable']),
                               'HST Graph Table')
    phdu.header['COMPTAB'] = (os.path.basename(d['comptable']),
                              'HST Components Table')
    phdu.header['USEAFTER'] = (useafter, 'Useafter date')
    phdu.header['PEDIGREE'] = (pedigree, 'Data pedigree')

    # Must be 67 char long
    phdu.header['DESCRIP'] = 'photometry keywords reference file' + 33 * '-'

    return phdu


def save_skipped_obsmodes(output, obsmodes):
    """Save skipped observation modes to a text file.
    Existing file is overwritten.

    Parameters
    ----------
    output : str
        IMPHTTAB filename. Output text file will
        replace ``_imp.fits`` with ``_skipped.txt``.

    obsmodes : list of str
        Skipped observation modes.

    """
    ind = output.find('_imp.fits')
    if ind != -1:
        output = output[:ind] + '_skipped.txt'
    else:
        output = output + '_skipped.txt'
    with open(output,'w') as f:
        for skipped in obsmodes:
            f.write(skipped + '\n')


def create_table(output, basemode, detector, useafter, tmgtab=None,
                 tmctab=None, tmttab=None, mode_list=[], nmodes=None,
                 clobber=False, verbose=True):
    """Create an IMPHTTAB reference file for a specified base
    configuration, ``basemode``.

    Parameters
    ----------
    output : str
        Output IMPHTTAB filename.
        (``_imp.fits`` will be appended if only prefix is given.)

    basemode : str
        Base observation mode for which to generate IMPHTTAB
        (e.g., ``acs,hrc``). This is ignored if ``mode_list``
        is given.

    detector : str
        Detector name.

    useafter : str
        Useafter date in the format of ``MMM DD YYYY HH:MM:SS``.

    tmgtab, tmctab, tmttab : str, optional
        Graph (TMG), component (TMC), and thermal component (TMT)
        tables to use. If `None`, the most recent version in CDBS
        is used.

    mode_list : list of str, optional
        A list of observation modes which should be used to make an
        IMPHTTAB reference file. If given, ``basemode`` is ignored.

    nmodes : int, optional
        Set to limit the number of modes to calculate.
        This is for testing only, otherwise set to `None`.

    clobber : bool, optional
        Overwrite existing IMPHTTAB?

    verbose : bool, optional
        Display extra information.


    """
    if not output.endswith('_imp.fits'):
        output = output + '_imp.fits'

    # If cannot overwrite, raise exception here rather than at the end.
    if os.path.exists(output) and not clobber:
        raise IOError('Output file already exists. Please delete/rename '
                      'before restarting.')

    # Define graph and component tables.
    if tmgtab is None:
        tmgtab = refs.GRAPHTABLE
    if tmctab is None:
        tmctab = refs.COMPTABLE
    if tmttab is None:
        tmttab = refs.THERMTABLE

    refs.setref(graphtable=tmgtab, comptable=tmctab, thermtable=tmttab)
    x = sgf.read_graphtable(tmgtab, tmctab, tmttab)

    if len(mode_list) == 0:
        # start by getting the full list of obsmodes before
        # expanding the parameterized elements
        x.get_obsmodes(basemode, prefix=True)
        obsmodes = x.obsmodes
    else:
        obsmodes = mode_list

    # start building obsmodes for each row
    if nmodes is not None:
        nrows = nmodes
        obsmodes = obsmodes[:nmodes]
    else:
        nrows = len(obsmodes)

    fpars_vals = list() # list of parameterized variable names for each obsmode/row
    npar_vals = list() # number of parameterized variables for each obsmode/row
    flam_datacol_vals = list()
    plam_datacol_vals = list()
    bw_datacol_vals = list()
    fpars_sz = 1

    # Compute 'globally' required values: max number of parameterized variables,...
    for filt in obsmodes:
        # For each filter combination (row in the table)...
        basename,fpars = parse_filters(filt)
        # keep track of how many parameterized variables are used in this obsmode
        npars = len(fpars)
        npar_vals.append(npars)
        fpars_vals.append(fpars)
        fpars_len = [len(f) for f in fpars]

        if len(fpars_len) == 0:
            fpars_max = 0
        else:
            fpars_max = max(fpars_len)

        if fpars_max > fpars_sz:
            fpars_sz = fpars_max

        if npars == 0:
            nstr = ''
        else:
            nstr = str(npars)

        flam_datacol_vals.append('PHOTFLAM' + nstr)
        plam_datacol_vals.append('PHOTPLAM' + nstr)
        bw_datacol_vals.append('PHOTBW' + nstr)
    #
    # At this point, all the interpretation for the following columns has been done:
    # OBSMODE, DATACOL (for all 3 tables), PEDIGREE, DESCRIP
    #
    # Start by determining the maximum number of parameters in any given obsmode
    max_npars = np.array(npar_vals,np.int32).max()
    log.info('MAX_NPARS: {0}   NROWS: {1}'.format(max_npars, nrows))

    #
    # Now, define empty lists for NELEM* and PAR*VALUES columns
    #
    nelem_rows = np.zeros([nrows,max_npars], np.int16) # nelem_rows[i] for each column i
    parvals_rows = list()
    filtdata_set = dict()
    parnames_rows = np.chararray([nrows,max_npars], itemsize=fpars_sz) # create columns for PAR*NAMES
    parnames_rows[:] = ''*fpars_sz # initialize with blanks, just to be safe

    for nr in xrange(nrows):
        # create path through graphtab for this obsmode, reading in values for
        # all parameterized variables as well
        obspath = x.traverse(obsmodes[nr], verbose=False)
        filtdata = obspath._params

        # Create a master set of parameterized variables and their ranges of values
        for p in filtdata:
            if (p.upper(),obsmodes[nr]) not in filtdata_set.keys():
                filtdata_set[(p.upper(),obsmodes[nr])] = filtdata[p]

        fpars = fpars_vals[nr]
        npars = npar_vals[nr]
        pvals = list()

        #extract entries from 'filtdata' for only the values given in 'fpars'
        for i in range(max_npars):
            if len(fpars) == 0:
                pvals.append(np.array([0]))
            else:
                if i < len(fpars):
                    f = fpars[i].upper()
                    nelem = len(filtdata[f])
                    nelem_rows[nr,i] = nelem
                    if filtdata[f] == '':
                        pvals.append(np.array([0]))
                    else:
                        pvals.append(np.array(filtdata[f]))
                        parnames_rows[nr,i] = f
                else:
                    pvals.append(np.array([0]))

        parvals_rows.append(pvals)

    #
    # All NELEM* and PAR*VALUES columns are correctly populated up to this point
    # in the code.
    #
    # Now, define the values for the actual results columns: PHOT*LAM, PHOTBW
    #
    flam_rows = list()
    plam_rows = list()
    bw_rows = list()

    nmode_vals = list()

    # dictionary to hold optical components
    # (pysynphot.observationmode._Component objects)
    component_dict = {}

    # list to hold skipped obsmodes
    skipped_obs = []

    if verbose:
        log.info('Computing photmetry values for each row\'s obsmode...')
        sys.stdout.flush()

    for nr in xrange(nrows):
        if verbose:
            log.info('Row: {0}'.format(nr + 1))  # Provide some indication of which row is being worked
            sys.stdout.flush()

        obsmode = obsmodes[nr]
        fpars = fpars_vals[nr]
        npars = npar_vals[nr]
        filtdict = OrderedDict()
        lenpars = list()
        for f in fpars:
            f = f.upper()
            filtdict[f] = filtdata_set[(f,obsmode)]
            lenpars.append(len(filtdict[f]))

        # Now build up list of all obsmodes with all combinations of
        # parameterized variables values
        olist = expand_obsmodes(obsmode,filtdict)

        nmodes = len(olist)

        # If there are warnings from here, try updating rules_dict
        # in graphtab.py
        if nmodes == 0 and verbose:
            log.warn('No info for {0}'.format(obsmode))

        pflam = np.zeros(nmodes, np.float64)
        pplam = np.zeros(nmodes, np.float64)
        pbw = np.zeros(nmodes, np.float64)

        skip = False

        for n,fullmode in enumerate(olist):
            try:
                value = compute_values(fullmode, component_dict)
            except ValueError as e:
                if e.message == 'Integrated flux is <= 0':
                    # integrated flux is zero, skip this obsmode
                    skip = True
                    skipped_obs.append(obsmode)

                    flam_datacol_vals.pop(nr)
                    plam_datacol_vals.pop(nr)
                    bw_datacol_vals.pop(nr)
                    parvals_rows.pop(nr)
                    nelem_rows = np.delete(nelem_rows, nr, 0)
                    parnames_rows = np.delete(parnames_rows, nr, 0)

                    if verbose:
                        log.info('\tSkipping {0}'.format(obsmode))

                    break
                elif e.message == 'math domain error':
                    skip = True
                    skipped_obs.append(obsmode)

                    if verbose:
                        log.info('\tSkipping {0}'.format(obsmode))

                    flam_datacol_vals.pop(nr)
                    plam_datacol_vals.pop(nr)
                    bw_datacol_vals.pop(nr)
                    parvals_rows.pop(nr)
                    nelem_rows = np.delete(nelem_rows, nr, 0)
                    parnames_rows = np.delete(parnames_rows, nr, 0)

                    break
                else:
                    raise

            if verbose:
                log.info('\tPHOTFLAM({0}) = {1}'.format(
                        fullmode, value['PHOTFLAM']))

            pflam[n] = value['PHOTFLAM']
            pplam[n] = value['PHOTPLAM']
            pbw[n] = value['PHOTBW']

        if skip is True:
            continue

        nmode_vals.append(nmodes)

        # Re-order results so that fastest varying variable is the last index
        # when accessed as a numpy array later by the C code
        photflam = ((pflam.reshape(lenpars)).transpose()).ravel()
        photplam = ((pplam.reshape(lenpars)).transpose()).ravel()
        photbw = ((pbw.reshape(lenpars)).transpose()).ravel()
        fvals = list()
        pvals = list()
        bvals = list()
        if npars == 0:
            fvals.append(photflam[0])
            pvals.append(photplam[0])
            bvals.append(photbw[0])
        else:
            fvals.append(0)
            pvals.append(0)
            bvals.append(0)
        for col in range(1,max_npars+1):
            if col == npars:
                fvals.append(np.array(photflam,np.float64))
                pvals.append(np.array(photplam,np.float64))
                bvals.append(np.array(photbw,np.float64))
            else:
                fvals.append(np.array([0]))
                pvals.append(np.array([0]))
                bvals.append(np.array([0]))
        flam_rows.append(fvals)
        plam_rows.append(pvals)
        bw_rows.append(bvals)

        del photflam,photplam,photbw,filtdict,lenpars

    del component_dict

    # remove any skipped obsmodes from the obsmodes list
    for sk in skipped_obs:
        obsmodes.remove(sk)

    # save skipped obs to a file
    if len(skipped_obs) > 0:
        save_skipped_obsmodes(output, skipped_obs)

    del skipped_obs

    log.info('Creating table columns from photometry values...')

    # Convert nelem information from row-oriented to column oriented
    nelem_cols = nelem_rows.transpose()
    parnames_cols = parnames_rows.transpose()

    parvals_cols = list()
    flam_cols = list()
    plam_cols = list()
    bw_cols = list()
    for col in range(max_npars):
        pvals = list()
        for row in range(len(parvals_rows)):
            pvals.append(parvals_rows[row][col])
        parvals_cols.append(pvals)

    for col in range(max_npars+1):
        fvals = list()
        plvals = list()
        bvals = list()
        for row in range(len(flam_rows)):
            fvals.append(flam_rows[row][col])
            plvals.append(plam_rows[row][col])
            bvals.append(bw_rows[row][col])
        if col == 0:
            fvals = np.array(fvals)
            plvals = np.array(plvals)
            bvals = np.array(bvals)
        flam_cols.append(fvals)
        plam_cols.append(plvals)
        bw_cols.append(bvals)

    ped_vals = [fits.getval(tmctab, 'pedigree', 0)] * len(nmode_vals)
    descrip_str = 'Generated {0} from {1}, mkimphttab version {2}, ' \
        'pysynphot version {3}'.format(
        get_date(), os.path.basename(tmgtab), __version__, S.__version__)
    descrip_vals = [descrip_str] * len(nmode_vals)

    # Finally, create the structures needed to define this row in the FITS table

    # Define each column in the table based on max_npars which are not different
    # from one extension to the other
    obsmode_col = fits.Column(name='OBSMODE', format='40A', array=obsmodes)
    pedigree_col = fits.Column(name='PEDIGREE', format='30A', array=ped_vals)
    descrip_col = fits.Column(name='DESCRIP', format='110A', array=descrip_vals)
    datacol_col = {}
    datacol_col['PHOTFLAM'] = fits.Column(
        name='DATACOL', format='12A', array=flam_datacol_vals)
    datacol_col['PHOTPLAM'] = fits.Column(
        name='DATACOL', format='12A', array=plam_datacol_vals)
    datacol_col['PHOTBW'] = fits.Column(
        name='DATACOL', format='12A', array=bw_datacol_vals)

    parvals_tabcols = list()
    nelem_tabcols = list()
    parnames_tabcols = list()
    parnames_format = str(fpars_sz) + "A"
    # for each parameterized element, create a set of columns specifying the
    # range of values for that parameter and the number of elements covering
    # that range namely, the PAR<n>VALUES and NELEM<n> columns
    for p in range(max_npars):
        nelem_tabcols.append(fits.Column(
            name="NELEM"+str(p+1), format="I", array=nelem_cols[p]))
        parvals_tabcols.append(fits.Column(
            name="PAR"+str(p+1)+"VALUES", format="PD()", array=parvals_cols[p]))
        parnames_tabcols.append(fits.Column(
            name="PAR"+str(p+1)+"NAMES", format=parnames_format,
            array=parnames_cols[p]))

    # create the set of results columns
    flam_tabcols = list()
    plam_tabcols = list()
    bw_tabcols = list()
    for p in range(max_npars + 1):
        if p == 0:
            format_str = 'D'
            pstr = ''
            fcols = flam_cols[p]
            pcols = plam_cols[p]
            bcols = bw_cols[p]
        else:
            format_str = 'PD()'
            pstr = str(p)
            fcols = flam_cols[p]
            pcols = plam_cols[p]
            bcols = bw_cols[p]

        flam_tabcols.append(fits.Column(
            name='PHOTFLAM'+pstr, format=format_str, array=fcols))
        plam_tabcols.append(fits.Column(
            name='PHOTPLAM'+pstr, format=format_str, array=pcols))
        bw_tabcols.append(fits.Column(
            name='PHOTBW'+pstr, format=format_str, array=bcols))

    # Now create the FITS file with the table in each extension
    log.info('Creating full table: {0}'.format(output))

    phdu = make_pri_hdu(output, max_npars, basemode.split(',')[0], detector,
                        ped_vals[0], useafter)

    flam_tab = fits.new_table(
        [obsmode_col, datacol_col['PHOTFLAM']] + flam_tabcols +
        parnames_tabcols + parvals_tabcols + nelem_tabcols +
        [pedigree_col, descrip_col])
    flam_tab.header['EXTNAME'] = ('PHOTFLAM', 'Extension name')
    flam_tab.header['EXTVER'] = (1, 'Extension number')

    plam_tab = fits.new_table(
        [obsmode_col, datacol_col['PHOTPLAM']] + plam_tabcols +
        parnames_tabcols + parvals_tabcols + nelem_tabcols +
        [pedigree_col, descrip_col])
    plam_tab.header['EXTNAME'] = ('PHOTPLAM', 'Extension name')
    plam_tab.header['EXTVER'] = (1, 'Extension number')

    bw_tab = fits.new_table(
        [obsmode_col, datacol_col['PHOTBW']] + bw_tabcols +
        parnames_tabcols + parvals_tabcols + nelem_tabcols +
        [pedigree_col, descrip_col])
    bw_tab.header['EXTNAME'] = ('PHOTBW', 'Extension name')
    bw_tab.header['EXTVER'] = (1, 'Extension number')

    ftab = fits.HDUList()
    ftab.append(phdu)
    ftab.append(flam_tab)
    ftab.append(plam_tab)
    ftab.append(bw_tab)
    ftab.writeto(output, clobber=clobber)


def create_nicmos_table(output, detector, useafter, pht_table, **kwargs):
    """Use a NICMOS ``_pht.fits`` table to generate an IMPHTTAB table
    for observation modes listed in the given table.

    Parameters
    ----------
    output, detector, useafter : str
        See :func:`create_table`.

    pht_table : str
        File name of ``_pht.fits`` table from which to take
        observation modes.

    kwargs : dict
        Keywords accepted by :func:`create_table`, except ``mode_list``.

    """
    with fits.open(pht_table) as pht:
        modes = np.char.strip(pht[1].data['photmode']).tolist()
    kwargs['mode_list'] = modes
    create_table(output, 'nicmos', **kwargs)


def create_table_from_table(output, useafter, imphttab, **kwargs):
    """Use a previously created IMPHTTAB reference file to generate a new
    IMPHTTAB reference file.

    Parameters
    ----------
    output, useafter : str
        See :func:`create_table`.

    imphttab : str
        File name of ``_imp.fits`` IMPHTTAB table from which to
        take observation modes.

    kwargs : dict
        Keywords accepted by :func:`create_table`, except ``mode_list``.


    """
    nextend=3 #default number of computed extensions
    extra_exten=list()
    
    with fits.open(imphttab) as imp:
        detector=imp[0].header['DETECTOR']
        basemode=''
        modes = np.char.strip(imp[1].data['OBSMODE']).tolist()
        #check if there are more than the 3 computed extensions
        nextend=imp[0].header['NEXTEND']
        print(nextend)
        if nextend > 3:
            for ext in range(4,nextend+1,1):
                extra_exten.append(imp[ext].header['extname'])            
    
    kwargs['mode_list'] = modes
    create_table(output,basemode,detector,useafter, **kwargs)
        
    #if there are more than 3 extensions assume they are static and append them to the output file
    if extra_exten:
        junkfile=tempfile.NamedTemporaryFile(dir='./',delete=False)
        print("Adding static extensions to output",extra_exten)
        with fits.open(imphttab) as imp:
            with fits.open(output) as outfile:
                for ext in range(4,nextend+1,1):
                    outfile.append(imp[ext])
                outfile.writeto(junkfile.name)
        os.remove(output)
        os.rename(junkfile.name,output)

