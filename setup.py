#!/usr/bin/env python

import glob

from setuptools import setup

scripts = glob.glob('bin/*')

description = "Simulated Absorption for Cosmology with Lyman-Alpha from the Yvette Mocks"

version="0.1"
setup(name="SaclayMocks",
    version=version,
    description=description,
    url="https://github.com/igmhub/SaclayMocks",
    author="Thomas Etourneau, Jean-Marc Le Goff et al",
    author_email="<your email here>",
    packages=['SaclayMocks'],
    package_dir = {'': 'py'},
    package_data = {'SaclayMocks': ['etc/']},
    install_requires=['numpy','scipy','iminuit','healpy','fitsio',
                      'numba','future','setuptools', 'pyfftw', 'pyigm'],
    test_suite='SaclayMocks.test.test_cor',
    scripts = scripts
    )
