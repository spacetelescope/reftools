"""This module contains a function that reads a graph table,
or set of segmented graph tables, and produces a Graph instance
that can be traversed.

"""
from __future__ import absolute_import, print_function, division

# STDLIB
import os
import re

# ASTROPY
from astropy.io import fits
from astropy import log

# LOCAL
from .graphtab import Graph, Node, Edge

USE_PYSYNPHOT = os.environ.get('MKIMPHTTAB_USE_PYSYNPHOT', '0') == '1'

__all__ = ['read_graphtable']

# search for string embedded in parentheses
__re_compvar_match = re.compile(r'\[(?P<varname>\w+#)\]$')


def read_graphtable(fname, tmcfile=None, tmtfile=None, offset=0, verbose=False):
    """Reads a file containing a graphtable and returns a Graph instance.
    Follows links to sub-graphs if present.

    offset is an optional parameter specifying the innode number
    at which you want the resulting Graph to start.

    """
    if USE_PYSYNPHOT:
        from pysynphot.locations import irafconvert
    else:
        from stsynphot.stio import irafconvert

    fname = irafconvert(fname)
    if tmcfile is not None:
        tmcfile = irafconvert(tmcfile)
    else:
        # try to get definition from header of graphtab itself
        tmcfile = irafconvert(fits.getval(fname, 'tmcfile'))
        if tmcfile is None:
            raise ValueError('No valid TMCFILE provided as input...')
    if tmtfile is not None:
        tmtfile = irafconvert(tmtfile)
    else:
        tmtfile = irafconvert(fits.getval(fname, 'tmtfile'))
        if tmtfile is None:
            raise ValueError('No valid TMTFILE provided as input...')

    # TODO: FITS vs ascii
    with fits.open(fname) as f:
        incol = f[1].data.field('innode')  # get the innode column here
        d = f[1].data[incol >= offset]  # only grab desired rows based on `offset`

        log.info('\nTMCFILE: {0}   TMTFILE: {1}\n'.format(tmcfile, tmtfile))

        # Preload empty nodes into the table, to simplify
        # what follows.
        # Preserve numerical order because of worries about "top"
        # in the recursive calls.
        numlist = list(set(d.field('innode')))  # make a list of node numbers for selected rows
        numlist.sort()
        nodelist = list()
        for num in numlist:
            nodelist.append(Node(num))
        del numlist

        ans = Graph(nodelist)
        del nodelist

        # Add some metadata
        ans.name = fname
        ans.tmcfiles.append(tmcfile)
        ans.tmtfiles.append(tmtfile)
        ans.area = f[0].header.get('effarea')    # in sq. cm.
        # fill the lookup
        ans.lookup = make_lookup(ans.tmcfiles + ans.tmtfiles)

        # Process each row of the table
        for row in d:
            # Get the innode
            innum = row.field('innode')
            node = ans.get(innum)

            if row.field('compname').endswith('_graph'):
                # Then we have to jump to another graph table.
                try:
                    subgraph = read_graphtable(row.field('compname'),
                                               offset=row.field('outnode'))
                    outnode = subgraph.top
                    ans.tmcfiles.extend(subgraph.tmcfiles)
                    ans.tmtfiles.extend(subgraph.tmtfiles)
                    check_area(ans, subgraph)
                except IOError as e:
                    log.warn("{0}: skipping".format(str(e)))
                    outnode = None

            else:
                # Get the outnode
                outnum = row.field('outnode')
                if outnum in ans:
                    outnode = ans[outnum]
                else:
                    outnode = None

            # If tmcfile has been defined, extract the throughput filename
            # for the component as well
            compkey = make_key(tmcfile, row.field('compname'))
            cname = ''
            parvar = ''
            if compkey in ans.lookup and compkey is not None:
                try:
                    cname = irafconvert(ans.lookup[compkey])
                except KeyError:
                    # No environment variable defined for this component, so do not
                    # try to read that component file
                    cname = ans.lookup[compkey]
                    if verbose:
                        log.warn('environment variable not defined for '
                                 '{0}'.format(cname))
                # extract the parameterized variable name from the compname
                varname = __re_compvar_match.findall(cname)
                if len(varname) > 0:
                    parvar = varname[0].upper()
                    cname = cname[:cname.find('[')]

            # Make the edge that corresponds to this row.
            # Compnames are qualified by their lookup table names.
            edge = Edge(row.field('keyword'),
                        [compkey,
                         make_key(tmtfile, row.field('thcompname')),
                         cname, parvar],
                        outnode)

            # and hook it up to the innode
            node.add_edge(edge)

    return ans


# The below are helper functions to handle the complexities of
# area handling and filename lookup with segmented files.
#
def check_area(main, sub):
    """ main = a Graph containing a sub-Graph
    sub = the contained sub-Graph
    Graph tables may specify the area used to convert to counts.
    This function checks that if both main and sub specify
    the area, they have the same value. If not, an exception
    is raised.
    """
    if main.area == sub.area:
        #all is well
        return
    elif main.area is not None and sub.area is None:
        #all is well: use main value
        return
    elif main.area is None and sub.area is not None:
        #then trickle the subarea up?
        #TODO: decide if this is the right thing to do
        main.area = sub.area
    else:
        #area mismatch, raise exception
        raise ValueError('Graph tables contain mismatched values for '
                         'effective area: {0}.effarea = {1}, {2}.effarea = '
                         '{3}'.format(main.name, main.area, sub.name, sub.area))


def make_key(fname, compname):
    """ fname = path to lookup file (tmc or tmt)
    compname = key in lookup file
    result = <basename of lookup>.<key>

    """
    if compname == 'clear':
        result = None
    else:
        result = '.'.join([os.path.basename(fname), compname])
    return result


def make_lookup(inlist):
    """ Each item in tmcfiles or tmtfiles is the
    name of a file containing a list. This method
    fills the lookup dictionary with the contents of
    all these files. A flattened dictionary with decorated
    keys is used rather than a hierarchical dict of dicts
    to simplify lookups.
    """
    lookup = dict()
    # "None" is a special key that should always return None
    lookup[None] = None

    # backwards compatibility: not all tables have these
    # header keywords set
    flist = [k for k in inlist if os.path.isfile(k)]
    # Process each file in the list
    for fname in flist:
        # Extract relevant data from file
        base = os.path.basename(fname)
        with fits.open(fname) as f:
            compname = f[1].data.field('compname')
            filename = f[1].data.field('filename')

        # Add items to the dictionary
        for (key, val) in zip(compname, filename):
            decorated_key = ".".join([base, key])
            lookup[decorated_key] = val

    if len(lookup) == 1:  # we put None in, remember
        log.warn('Warning, no lookup tables found. They must be filled by '
                 'hand.')

    return lookup
