"""
.. note::

  Some tests are marked as expected failures but they should be fixed
  if ``getphotpars`` is still used.

"""
import pytest
from astropy.utils.data import get_pkg_data_filename
from numpy.testing import assert_allclose

from .. import getphotpars


class TestGetPhotPars:
    def setup_class(self):
        self.get_pars = getphotpars.GetPhotPars(
            get_pkg_data_filename('data/test_wfc1_dev_imp.fits'))

    def teardown_class(self):
        self.get_pars.close()

    def test_get_row_len_obs(self):
        obsmode = 'acs,wfc1,f625w,f814w,MJD#'
        row = self.get_pars._get_row(obsmode, 'photflam')
        assert len(row[0]) == 13
        assert row['obsmode'] == obsmode

    def test_parse_obsmode_and_make_row_struct_and_compute_value_0(self):
        obsmode = 'acs,wfc1,f625w,f660n'
        ext = 'photflam'
        npars, strp_obsmode, par_dict = self.get_pars._parse_obsmode(obsmode)
        row = self.get_pars._get_row(strp_obsmode, ext)
        rd = self.get_pars._make_row_struct(row, npars)
        ps = self.get_pars._make_par_struct(npars, par_dict)

        assert npars == 0
        assert strp_obsmode == obsmode
        assert len(par_dict) == 0

        assert rd['obsmode'].lower() == obsmode
        assert rd['datacol'].lower() == ext
        assert len(rd['parnames']) == 0
        assert rd['parnum'] == 0
        assert_allclose(rd['results'], 5.8962401031019617e-18)
        assert rd['telem'] == 1
        assert len(rd['nelem']) == 0
        assert len(rd['parvals']) == 0

        assert ps['npar'] == 0
        assert len(ps['parnames']) == 0
        assert len(ps['parvals']) == 0

        assert_allclose(self.get_pars._compute_value(rd, ps),
                        5.8962401031019617e-18)

        rd_2 = self.get_pars._make_row_struct(self.get_pars._get_row(
            strp_obsmode, 'photplam'), npars)
        assert_allclose(self.get_pars._compute_value(rd_2, ps),
                        6599.6045327828697)

        rd_3 = self.get_pars._make_row_struct(self.get_pars._get_row(
            strp_obsmode, 'photbw'), npars)
        assert_allclose(self.get_pars._compute_value(rd_3, ps),
                        13.622138313347964)

    def test_parse_obsmode_and_make_par_struct_1(self):
        npars, strp_obsmode, par_dict = self.get_pars._parse_obsmode(
            'acs,wfc1,f625w,f814w,MJD#55000.0')
        assert npars == 1
        assert strp_obsmode == 'acs,wfc1,f625w,f814w,mjd#'
        assert list(par_dict)[0] == 'mjd#'
        assert par_dict['mjd#'] == 55000

        ps = self.get_pars._make_par_struct(npars, par_dict)
        assert ps['npar'] == 1
        assert ps['parnames'] == ['mjd#']
        assert_allclose(ps['parvals'], [55000])

    def test_parse_obsmode_and_make_par_struct_2(self):
        npars, strp_obsmode, par_dict = self.get_pars._parse_obsmode(
            'acs,wfc1,f625w,fr505n#5000.0,MJD#55000.0')
        keys = list(par_dict)
        assert npars == 2
        assert strp_obsmode == 'acs,wfc1,f625w,fr505n#,mjd#'
        assert 'mjd#' in keys
        assert 'fr505n#' in keys
        assert par_dict['mjd#'] == 55000
        assert par_dict['fr505n#'] == 5000

        ps = self.get_pars._make_par_struct(npars, par_dict)
        assert ps['npar'] == 2
        assert sorted(ps['parnames']) == ['fr505n#', 'mjd#']
        assert_allclose(sorted(ps['parvals']), [5000, 55000])

    @pytest.mark.xfail(reason='new type not compatible with array')
    def test_make_row_struct_1(self):
        obsmode = 'acs,wfc1,f625w,f814w,MJD#55000.0'
        npars, strp_obsmode, par_dict = self.get_pars._parse_obsmode(obsmode)
        row = self.get_pars._get_row(strp_obsmode, 'photflam')
        rd = self.get_pars._make_row_struct(row, npars)
        assert rd['obsmode'].lower() == strp_obsmode.lower()
        assert rd['datacol'].lower() == 'photflam1'
        assert rd['parnames'] == ['MJD#']
        assert rd['parnum'] == 1
        assert_allclose(rd['results'],
                        [8.353948620387228e-18, 8.353948620387228e-18,
                         8.439405629259882e-18, 8.439405629259882e-18])
        assert rd['telem'] == 4
        assert rd['nelem'] == [4]
        assert_allclose(rd['parvals'], [[52334.0, 53919.0, 53920.0, 55516.0]])

    @pytest.mark.xfail(reason='new type not compatible with array')
    def test_make_row_struct_2(self):
        obsmode = 'acs,wfc1,f625w,fr505n#5000.0,MJD#55000.0'
        npars, strp_obsmode, par_dict = self.get_pars._parse_obsmode(obsmode)
        row = self.get_pars._get_row(strp_obsmode, 'photflam')
        rd = self.get_pars._make_row_struct(row, npars)

        assert rd['obsmode'].lower() == strp_obsmode.lower()
        assert rd['datacol'].lower() == 'photflam2'
        assert rd['parnames'] == ['FR505N#', 'MJD#']
        assert rd['parnum'] == 2
        assert rd['telem'] == 44
        assert rd['nelem'] == [11, 4]

        ans = [
            7.003710657512407e-14, 6.944751784992699e-14, 5.933229935875258e-14,
            6.903709791285467e-14, 6.679623972708115e-14, 5.70322329920528e-14,
            6.871738863125632e-14, 6.430043639034794e-14, 5.2999039030883145e-14,  # noqa
            6.984240267831951e-14, 6.158040091027122e-14, 6.903709791285467e-14,
            6.679623972708115e-14, 5.70322329920528e-14, 6.777606555875925e-14,
            6.344417923365138e-14, 5.230630810980452e-14, 6.984240267831951e-14,
            6.158040091027122e-14, 7.102142054088e-14, 7.038429771136416e-14,
            6.012566486354809e-14, 6.777606555875925e-14, 6.344417923365138e-14,
            5.230630810980452e-14, 6.889989612586935e-14, 6.076197306179641e-14,
            7.102142054088e-14, 7.038429771136416e-14, 6.012566486354809e-14,
            7.000118709415907e-14, 6.7696605688248e-14, 5.778989112976214e-14,
            6.889989612586935e-14, 6.076197306179641e-14, 7.003710657512407e-14,
            6.944751784992699e-14, 5.933229935875258e-14, 7.000118709415907e-14,
            6.7696605688248e-14, 5.778989112976214e-14, 6.871738863125632e-14,
            6.430043639034794e-14, 5.2999039030883145e-14]
        assert_allclose(rd['results'], ans)
        assert_allclose(
            rd['parvals'],
            [[4824, 4868.2, 4912.4, 4956.6, 5000.8, 5045, 5089.2, 5133.4,
              5177.6, 5221.8, 5266],
             [52334, 53919, 53920, 55516]])

    @pytest.mark.xfail(reason='new type not compatible with array')
    def test_compute_value_1(self):
        npars, strp_obsmode, par_dict = self.get_pars._parse_obsmode(
            'acs,wfc1,f625w,f814w,MJD#55000.0')
        ps = self.get_pars._make_par_struct(npars, par_dict)

        ext = 'photflam'
        row = self.get_pars._get_row(strp_obsmode, ext)
        rd = self.get_pars._make_row_struct(row, npars)
        result = self.get_pars._compute_value(rd, ps)
        assert_allclose(result, 8.43940563e-18)

        ext = 'photplam'
        row = self.get_pars._get_row(strp_obsmode, ext)
        rd = self.get_pars._make_row_struct(row, npars)
        result = self.get_pars._compute_value(rd, ps)
        assert_allclose(result, 6992.37762323)

        ext = 'photbw'
        row = self.get_pars._get_row(strp_obsmode, ext)
        rd = self.get_pars._make_row_struct(row, npars)
        result = self.get_pars._compute_value(rd, ps)
        assert_allclose(result, 58.85223114)

    @pytest.mark.xfail(reason='new type not compatible with array')
    def test_compute_value_2(self):
        npars, strp_obsmode, par_dict = self.get_pars._parse_obsmode(
            'acs,wfc1,f625w,fr505n#5000.0,MJD#55000.0')
        ps = self.get_pars._make_par_struct(npars, par_dict)

        ext = 'photflam'
        row = self.get_pars._get_row(strp_obsmode, ext)
        rd = self.get_pars._make_row_struct(row, npars)
        result = self.get_pars._compute_value(rd, ps)
        assert_allclose(result, 5.99660350e-14)

        ext = 'photplam'
        row = self.get_pars._get_row(strp_obsmode, ext)
        rd = self.get_pars._make_row_struct(row, npars)
        result = self.get_pars._compute_value(rd, ps)
        assert_allclose(result, 5737.95131007)

        ext = 'photbw'
        row = self.get_pars._get_row(strp_obsmode, ext)
        rd = self.get_pars._make_row_struct(row, npars)
        result = self.get_pars._compute_value(rd, ps)
        assert_allclose(result, 643.42528984)


@pytest.mark.parametrize(
    ('obsmode', 'imphttab', 'photflam', 'photplam', 'photbw'),
    [('acs,wfc1,f625w,f660n', 'data/test_wfc1_dev_imp.fits',
      5.8962401031019617e-18, 6599.6045327828697, 13.622138313347964),
     ('acs,wfc1,f625w,f660n', 'data/test_acs_wfc1_dev_imp.fits',
      5.8962401031019617e-18, 6599.6045327828697, 13.622138313347964)])
def test_get_phot_pars_func(obsmode, imphttab, photflam, photplam, photbw):
    results = getphotpars.get_phot_pars(
        obsmode, get_pkg_data_filename(imphttab))
    assert results["PHOTZPT"] == -21.1
    assert_allclose(results["PHOTFLAM"], photflam)
    assert_allclose(results["PHOTPLAM"], photplam)
    assert_allclose(results["PHOTBW"], photbw)


# TODO: Merge this with test_get_phot_pars_func() when fixed.
@pytest.mark.xfail(reason='new type not compatible with array')
@pytest.mark.parametrize(
    ('obsmode', 'imphttab', 'photflam', 'photplam', 'photbw'),
    [('acs,wfc1,f625w,f814w,MJD#55000.0', 'data/test_wfc1_dev_imp.fits',
      8.43940563e-18, 6992.37762323, 58.85223114),
     ('acs,wfc1,f625w,f814w,MJD#55000.0', 'data/test_acs_wfc1_dev_imp.fits',
      8.43940563e-18, 6992.37762323, 58.85223114),
     ('acs,wfc1,f625w,fr505n#5000.0,MJD#55000.0',
      'data/test_wfc1_dev_imp.fits',
      5.99660350e-14, 5737.95131007, 643.42528984),
     ('acs,wfc1,f625w,fr505n#5000.0,MJD#55000.0',
      'data/test_acs_wfc1_dev_imp.fits',
      5.99660350e-14, 5737.95131007, 643.42528984),
     ('acs,wfc1,fr931n#8905', 'data/test_acs_wfc1_dev_imp.fits',
      1.4852052585262792e-18, 8903.8757202468823, 64.255181883310669)])
def test_get_phot_pars_func_2(obsmode, imphttab, photflam, photplam, photbw):
    results = getphotpars.get_phot_pars(
        obsmode, get_pkg_data_filename(imphttab))
    assert results["PHOTZPT"] == -21.1
    assert_allclose(results["PHOTFLAM"], photflam)
    assert_allclose(results["PHOTPLAM"], photplam)
    assert_allclose(results["PHOTBW"], photbw)
