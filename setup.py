#!/usr/bin/env python

from distutils.core import setup


setup_kwargs = {}

setup_kwargs[ 'setup_requires' ] = [ 'numpy==1.9.2', 'cython==0.20.1' ]
setup_kwargs[ 'install_requires' ] = [ 'numpy==1.9.2',
        'cython==0.20.1',
        'matplotlib>=1.4.3',
        'loprop>=0.1',
        'pd>=0.1', ]
setup_kwargs[ 'dependency_links' ] = [ 'git+https://github.com/fishstamp82/loprop.git@e0eb1b2796da77b44177138bc84c2f44329e9d44#egg=loprop-0.1',
        'git+https://github.com/fishstamp82/pd.git@af469199c15f7ca73f770ad9ae4527753a7b1928#egg=pd-0.1' ]

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
