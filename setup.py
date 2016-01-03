#!/usr/bin/env python
from distutils.core import setup

setup(name="moltools",
    version="1.0",
    packages=["src/test", "src/loprop", "src/loprop/daltools", "src/loprop/daltools/util", "src/pd" ],
    scripts=["src/molecules.py", "src/pdbreader.py", ],
    author="Ignat Harczuk",
    author_email="harczuk@kth.se",
    )



