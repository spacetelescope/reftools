#!/usr/bin/env python
from glob import glob
from setuptools import setup, Extension

MACROS = [
    ("Py_LIMITED_API", 0x03090000),  # PY_VERSION_HEX for 3.11
]

setup(
    use_scm_version={"write_to": "reftools/version.py"},
    ext_modules=[
        Extension(
            "reftools._computephotpars", glob("reftools/src/*.c"), define_macros=MACROS
        )
    ],
    options={'bdist_wheel': {'py_limited_api': 'cp39'}},
)
