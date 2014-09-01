#!/usr/bin/env python
#-*- coding: utf-8 -*-

import numpy

def main():
    num = 2
    cnt = 0

    tmp = []
    while cnt < num:
        for i in range(52):
            tmp.append( i )
        cnt += 1
    print tmp

if __name__ == '__main__':
    main()
