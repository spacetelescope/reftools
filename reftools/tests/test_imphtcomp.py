from astropy.utils.data import get_pkg_data_filename
from numpy.testing import assert_allclose

from .. import imphtcomp


def test_same_file_diffs():
    """Test comparing a file to itself."""
    input_file = get_pkg_data_filename('data/acs_sbc_dev_imp.fits')
    comp = imphtcomp.ImphttabComp(input_file, input_file)
    assert_allclose(comp.flamdiff, 0)
    assert_allclose(comp.plamdiff, 0)
    assert_allclose(comp.bwdiff, 0)
