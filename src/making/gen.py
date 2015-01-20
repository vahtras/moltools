#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from molecules import *
from use_generator import *


g = Generator()

m1 = g.get_mol( model = 'spc')
open('spc.mol','w').write( m1.get_mol_string())

