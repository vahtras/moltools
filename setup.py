#!/usr/bin/env python

from distutils.core import setup


setup_kwargs = {}

setup_kwargs[ 'setup_requires' ] = [ 'numpy>=1.9.2', 'cython>=0.20.1' ]
setup_kwargs[ 'install_requires' ] = [ 'numpy>=1.9.2',
        'cython>=0.20.1',
        'matplotlib>=1.4.3',
        'loprop>=0.1',
        'nose',
        'applequist>=0.1',
        ]
try:
    import pkg_resources
    from setuptools import setup, Command
    _have_setuptools = True
except:
    from distutils.core import setup, Command
    _have_setuptools = False

setup(name="moltools",
    version="1.0",
    packages=[ "moltools", ],
    scripts=[ "moltools/pdbreader.py", ],
    author="Ignat Harczuk",
    author_email="harczuk@kth.se",
    license = 'MIT',
    description = 'Running point dipole calculations on chemical systems',
    platforms = 'Linux',
    maintainer_email = 'ignathe@gmail.com',
    url = 'https://github.org/fishstamp82/moltools.git',
    **setup_kwargs
    )
