#!/usr/bin/env python
from glob import glob
from setuptools import setup, Extension

setup(
    use_scm_version={'write_to': 'reftools/version.py'},
    ext_modules=[
        Extension('reftools._computephotpars', glob('reftools/src/*.c'))
    ]
)
