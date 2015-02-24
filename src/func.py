#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re, sys, os
from matplotlib import pyplot as plt
import numpy as np
from read_dal import o_filter
from gaussian import *

outs = [f for f in os.listdir( os.getcwd()) if f.endswith('.out') ]

def run_mpl(
        x_dim = "r",
        y_dims = ["rel","qm_abs"],
        comps = ["x"],
        props = ["d"],
        levels = [0],
        mol_model = "TIP3P",
        freqs = ["0.0"],
        models = ["gaussian"],
        max_l = [1],
        basis = "ANOPVDZ",
        Rqs = [0.0001],
        Rps = [0.0001],
        in_AA = False,
        out_AA = False,
        ):
    for i in enumerate( [comps, props] ):
        print i

if __name__ == '__main__':
    run_mpl()
