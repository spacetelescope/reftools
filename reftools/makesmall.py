import os

import numpy as np
from astropy.io import fits

__version__ = '0.3.1'
__vdate__ = '2010-04-26'
__all__ = ['dxy']


def dxy(dgeofile, filter=None, colcorr=None, corrext='DX', minsize=32,
        debug=False):
    """Build subsampled CDBS ``_dxy`` files from full-frame (Jay's)
    DXM and DYM images.

    If there is a column/row correction provided, it will be removed
    before resampling the full DXM,DYM images. In that case, the
    ``corrext`` parameter specifies the extension from the full-sized
    image that needs to have this column/row correction removed.

    """
    from stsci.tools import fileutil

    wheel1 = ['F475W', 'F555W', 'F606W', 'F850LP', 'F625W',
              'F658N', 'F775W', 'F502N', 'F550M']
    wheel2 = ['F435W', 'F814W', 'F660N']

    odgeofile = dgeofile
    # Interpret input filename, expanding any variables used for path to file
    dgeofile = fileutil.osfn(dgeofile)
    dgeoroot = os.path.split(dgeofile)[1]
    if not os.path.exists(dgeofile):
        raise ValueError(f'{dgeofile} is not a valid input')

    # If column/row correction image is specified, open it for removal from
    # DGEOFILE.
    if colcorr is not None:
        if not os.path.exists(colcorr):
            raise ValueError('No valid column/row correction image found!')

        # User specified a column/row correction image to be subtracted from
        # full DGEOFILE before resampling
        corrimg = fits.open(colcorr, 'readonly')
    else:
        print('======')
        print('No column/row correction removed from DGEOFILE before '
              'resampling.')
        print('======')

    # Open file for processing
    dxyfile = fits.open(dgeofile, 'readonly')  # Model _DXY file
    detector = dxyfile[0].header['detector']
    print(f'DETECTOR {detector}')
    print('')

    # Read in filter name from DGEOFILE directly, if not specified by user
    filter1 = None
    filter2 = None
    if filter is None:
        if 'FILTER1' in dxyfile[0].header:
            filter1 = dxyfile[0].header['FILTER1']
            filter2 = dxyfile[0].header['FILTER2']
            if filter2 == 'N/A':
                filter2 = 'CLEAR2S'
            if filter2.find('CLEAR') > -1:
                filter = filter1
            else:
                filter = filter2
        else:
            # Error case, filter name not provided in file or by user
            dxyfile.close()
            raise ValueError('Filter name not found in DGEOFILE')

    filter = filter.upper()  # On input, case can be upper or lower
    print(f'Filter {filter}')

    # Get the shape of each chip from the first DGEOFILE array
    dxy_shape = dxyfile[1].data.shape

    # compute step size needed to create DXY file no smaller than 32x32
    stepsize = min(dxy_shape[1] / minsize, dxy_shape[0] / minsize)
    # create the grid in such a way that it contains the idices of the
    # first and last element of the full size array (0-4096) inclusive
    grid = [dxy_shape[1], dxy_shape[0], stepsize, stepsize]
    xpts = np.array(list(range(0, grid[0] + 1, grid[2])), np.float32)
    ypts = np.array(list(range(0, grid[1] + 1, grid[3])), np.float32)
    xygrid = np.meshgrid(xpts, ypts)
    # this padding is necessary so that the end points in the small npol file
    # match the end point of the large dgeo file.
    if xygrid[0][:, -1][-1] >= dxy_shape[1]:
        xygrid[0][:, -1] = [dxy_shape[1] - 1] * len(xygrid[0][:, -1])
    if xygrid[1][-1][-1] >= dxy_shape[0]:
        xygrid[1][-1] = [dxy_shape[0] - 1] * len(xygrid[1][-1])

    # count the number of chips in DGEOFILE
    numchips = 0
    extname = dxyfile[1].header['EXTNAME']

    for extn in dxyfile:
        if 'extname' in extn.header and extn.header['extname'] == extname:
            numchips += 1

    # process each chip
    for chip in range(1, numchips + 1):
        # Process DX and DY for each chip
        for xy in ['DX', 'DY']:
            print(f'Processing chip from extension: {xy},{chip}')
            onaxis1 = dxyfile[xy, chip].header['NAXIS1']
            onaxis2 = dxyfile[xy, chip].header['NAXIS2']
            if 'CCDCHIP' not in dxyfile[xy, chip].header:
                ccdchip = 1
            else:
                ccdchip = dxyfile[xy, chip].header['CCDCHIP']
            cdelt1 = grid[2]
            cdelt2 = grid[3]

            dxy = dxyfile[xy, chip].data.copy()
            if colcorr is not None and xy == corrext:
                # Remove column/row correction from this extension
                xycorr = corrimg[0].data
            else:
                xycorr = None

            # define default chipname for debug results; None indicates no
            # debugging
            cname = None
            if debug:
                cname = f'full_68corr_{xy}{chip}_dxy.fits'

            # CDELT1/2 are the stepsize used to create the npl files.
            # It's best to write cdelt in the headers of the npl files instead
            # of try to calculate cdelt from naxis and onaxis.
            # Keep onaxis for reference.
            dxyfile[xy, chip].data = resample_chip(
                dxy, xygrid, corrxy=xycorr, chipname=cname)
            dxyfile[xy, chip].header.update(
                'ONAXIS1', onaxis1, 'NAXIS1 of full size dgeo file')
            dxyfile[xy, chip].header.update(
                'ONAXIS2', onaxis2, 'NAXIS2 of full size dgeo file')
            dxyfile[xy, chip].header.update(
                'CCDCHIP', ccdchip, 'CCDCHIP from full size dgeo file')
            dxyfile[xy, chip].header.update(
                'CDELT1', cdelt1, 'Coordinate increment along axis')
            dxyfile[xy, chip].header.update(
                'CDELT2', cdelt2, 'Coordinate increment along axis')

    # Get filter info ready for use in updating output image header
    if filter1 is None:
        if filter in wheel1:
            filter1 = filter
            filter2 = 'CLEAR2L'

        if filter in wheel2:
            filter2 = filter
            filter1 = 'CLEAR1L'

    print(filter1, filter2)

    # Update keywords
    newname = dgeoroot[:dgeoroot.find('_dxy.fits')] + '_npl.fits'
    if os.path.exists(newname):
        os.remove(newname)
    dxyfile[0].header['filename'] = newname
    dxyfile[0].header['filter1'] = filter1
    dxyfile[0].header['filter2'] = filter2
    dxyfile[0].header.update('pedigree', 'INFLIGHT 01/03/2002 01/10/2005')
    dxyfile[0].header.update('date', fileutil.getDate())
    dxyfile[0].header.add_history(f'File generated from DGEOFILE: {odgeofile}',
                                  after='pedigree')
    dxyfile.writeto(newname)

    # close open file handles
    if colcorr is not None:
        corrimg.close()
    dxyfile.close()

    # finish by cleaning up output image header by removing unnecessary
    # keywords
    print('Cleaning up header')
    fixheader(filter, newname, odgeofile)
    print(f'Finished creating new file: {newname}')


def resample_chip(chip, mesh, corrxy=None, chipname=None, overwrite=True):
    """Apply resampling operations to a single chip's array.
    If provided, also remove column/row correction from array prior to
    resampling.

    """
    if corrxy is not None:
        print(f'Removing column correction {corrxy.shape} from full chip '
              'dxy array...')
        chip -= corrxy
    if chipname is not None:
        outimg = fits.PrimaryHDU(data=chip)
        outimg.writeto(chipname, overwrite=overwrite)

    return chip[mesh[1].astype(np.int), mesh[0].astype(np.int)]


def fixheader(filter, filename, oldname):
    """Clean up dxy headers."""
    required_keys = ['INSTRUME', 'DETECTOR', 'FILTER1', 'FILTER2',
                     'COMMENT', 'HISTORY', 'PROPOSID', 'CAL_VER']
    with fits.open(filename, mode='update') as dxyfile:
        h0 = dxyfile[0].header
        # remove extraneous keywords from original DGEOFILE
        for h0key in list(h0.keys())[15:]:
            if h0key not in required_keys:
                del h0[h0key]
        del h0['']
        h0.update('FILETYPE', 'DXY GRID')
        h0.update('COMMENT',
                  'Accuracy to 0.01 pixels when dxy corrections included')
        h0.add_history('Improved solution as reported in 2005 Cal Workshop')
        h0.add_history(
            f'Average of 64x64 blocks from full DXY image {oldname}')
