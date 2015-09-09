"""Test parsedq module."""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from reftools.parsedq import parse_one_dq


def test_one_flag():
    dqs = parse_one_dq(16658)
    assert sorted(dqs) == [2, 16, 256, 16384]
