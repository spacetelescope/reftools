#!/usr/bin/env python
import sysconfig
from glob import glob

from setuptools import Extension, setup

FREE_THREADED_PYTHON = sysconfig.get_config_var("Py_GIL_DISABLED") == 1

MACROS = []
if not FREE_THREADED_PYTHON:
    MACROS.append(("Py_LIMITED_API", 0x03090000))  # PY_VERSION_HEX for 3.9

SETUPTOOLS_OPTIONS = {}
if not FREE_THREADED_PYTHON:
    SETUPTOOLS_OPTIONS["bdist_wheel"] = {"py_limited_api": "cp39"}

setup(
    use_scm_version={"write_to": "reftools/version.py"},
    ext_modules=[
        Extension(
            "reftools._computephotpars",
            glob("reftools/src/*.c"),
            define_macros=MACROS,
            py_limited_api=not FREE_THREADED_PYTHON,
        )
    ],
    options=SETUPTOOLS_OPTIONS,
)
