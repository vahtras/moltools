#!/usr/bin/env python
from distutils.core import setup

setup(name="moltools",
    version="1.0",
    packages=["src", "src/test", "src/loprop", ],
    scripts=["src/molecules.py", "src/pdbreader.py", ],
    author="Ignat Harczuk",
    author_email="harczuk@kth.se",
    )



