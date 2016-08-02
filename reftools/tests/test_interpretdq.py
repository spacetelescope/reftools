"""Test interpretdq module."""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ..interpretdq import DQParser

__all__ = ['test_one_flag']


def test_one_flag():
    parsedq = DQParser.from_instrument(None)  # Default
    dqs = parsedq.interpret_dqval(16658)
    assert sorted(dqs['DQFLAG']) == [2, 16, 256, 16384]
