"""Test interpretdq module."""

import pytest
from astropy.utils.exceptions import AstropyUserWarning

from ..interpretdq import DQParser


def test_one_flag():
    with pytest.warns(AstropyUserWarning, match='using default'):
        parsedq = DQParser.from_instrument(None)
    dqs = parsedq.interpret_dqval(16658)
    assert sorted(dqs['DQFLAG']) == [2, 16, 256, 16384]
