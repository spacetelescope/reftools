#!/usr/bin/env python
import sys
from glob import glob
from setuptools import setup, Extension

# Submodule
sys.path.insert(1, 'relic')
import relic.release

version = relic.release.get_info()
relic.release.write_template(version, 'reftools')

setup(
    name='reftools',
    version=version.pep386,
    author=('Warren Hack, Nadezhda Dencheva, Vicki Laidler, Matt Davis, '
            'Megan Sosey, Pey Lian Lim, Mihai Cara'),
    author_email='help@stsci.edu',
    description=('Set of tools used to support creation of calibration '
                 'reference files for Hubble Space Telescope'),
    url='https://github.com/spacetelescope/reftools',
    license='BSD',
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    install_requires=[
        'astropy>=1.1',
        'numpy'
    ],
    tests_require=['nose'],
    packages=['reftools', 'reftools.tests'],
    package_dir={'reftools': 'reftools'},
    package_data={
        'reftools': [
            'data/*.*',
            'data/pctetab/*.*',
            'pars/*.*',
            '*.help'
        ],
        'reftools.tests': ['data/*.*']
    },
    entry_points={
        'console_scripts': ['tdspysyn=reftools.tdspysyn:main']
    },
    ext_modules=[
        Extension('reftools._computephotpars', glob('reftools/src/*.c'))
    ]
)
