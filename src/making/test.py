#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from molecules import *

def main():
    ref = []
    mine = []
    for i in open('ref').readlines():
        isp = i.split()
        ref.append( [ isp[0], isp[1], isp[2], isp[3] ] )

    for i in open('mine').readlines():
        isp = i.split()
        mine.append( [isp[0], isp[1], isp[2], isp[3] ] )

    for i in ref:
        lis, rev = i, list(reversed(i))
        #raise SystemExit
    print "ref not in mine:"
    for i in ref:
        if (i not in mine) and (list(reversed(i)) not in mine):
            print i
    print "mine not in ref:"
    for i in mine:
        if (i not in ref) and (list(reversed(i)) not in ref):
            print i


if __name__ == '__main__':
    main()

