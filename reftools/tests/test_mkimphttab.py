import pytest
pytest.importorskip('stsynphot')

from astropy import log  # noqa
from astropy.io.fits.verify import VerifyWarning  # noqa

from .. import mkimphttab  # noqa


def setup_module():
    # Do not spam pytest log
    log.setLevel('ERROR')


def teardown_module():
    log.setLevel('INFO')


@pytest.mark.remote_data
def test_mkimphttab(tmpdir):
    """Make sure that when making IMPHTTAB, there are no errors."""
    output = str(tmpdir.join('test_out_imp.fits'))
    with pytest.warns(VerifyWarning, match='comment will be truncated'):
        mkimphttab.create_table(output, 'acs,sbc', 'sbc', 'DUMMY')
