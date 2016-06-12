#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from setuptools import setup, find_packages
import os
import io

import ebsel

here = os.path.abspath(os.path.dirname(__file__))


def read(*filenames, **kwargs):
    encoding = kwargs.get('encoding', 'utf-8')
    sep = kwargs.get('sep', '\n')
    buf = []
    for filename in filenames:
        with io.open(os.path.join(here, filename), encoding=encoding) as f:
            buf.append(f.read())
    return sep.join(buf)

long_description = read('README.md')

setup(
    name='ebsel',
    version=ebsel.version,
    url='https://github.com/jaimergp/ebsel',
    license='MIT',
    author='Jaime Rodriguez-Guerra Pedregal',
    author_email='jaime.rogue@gmail.com',
    description='Emsl basis set exchange local (Py3 fork)',
    long_description=long_description,
    packages=find_packages(),
    include_package_data=True,
    platforms='any',
    classifiers=[
        'Programming Language :: Python',
        'Development Status :: 2 - Pre-Alpha',
        'Natural Language :: English',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: Free For Educational Use',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Chemistry',
    ],
    scripts=['scripts/EMSL_api.py', 'scripts/basis_data_from_logs.py']
    ,
)
