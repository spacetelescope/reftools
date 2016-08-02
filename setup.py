#!/usr/bin/env python
import os
import subprocess
import sys
from glob import glob
from setuptools import setup, find_packages, Extension

# Submodule
sys.path.insert(1, 'relic')
import relic.release

version = relic.release.get_info()
relic.release.write_template(version, 'lib/reftools')

setup(
    name='reftools',
    version=version.pep386,
    author=('Warren Hack, Nadezhda Dencheva, Vicki Laidler, Matt Davis, '
            'Megan Sosey, Pey Lian Lim, Mihai Cara'),
    author_email='help@stsci.edu',
    description=('Set of tools used to support creation of calibration '
                 'reference files for Hubble Space Telescope'),
    url='https://github.com/spacetelescope/reftools',
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    install_requires=[
        'astropy',
        'calcos',
        'nose',
        'numpy',
        'scipy',
        'sphinx',
        'stsci.imagestats',
        'stsci.sphinxext',
        'stsci.tools',
        'stwcs'
    ],
    package_dir={
        '': 'lib'
    },
    packages=find_packages('lib'),
    package_data={
        'reftools': [
            'data/*',
            'pars/*',
            '*.help',
            '*.rules',
            'LICENSE.txt',
        ]
    },
    entry_points={
        'console_scripts': ['tdspysyn=reftools.tdspysyn:main']
    },
    ext_modules=[
        Extension('reftools._computephotpars', glob('src/*.c'))
    ]
)
