#!/usr/bin/env python
from glob import glob
from setuptools import setup, Extension

setup(
    name='reftools',
    use_scm_version=True,
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
    setup_requires=['setuptools_scm'],
    install_requires=[
        'astropy>=2',
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
