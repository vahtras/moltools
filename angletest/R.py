#!/usr/bin/env python
#-*- coding: utf-8 -*-

import numpy as np

class R:
    def __init__(self):
        import rpy2.robjects as robj
        self.obj = robj
    def __call__(self):
        return self.obj.r
    def getRVectors(self, *args ):

        tmp = []
        if args is not None:
            for i in args:
                if type (i) == list:
                    if type( i[0] ) == float:
                        tmp.append( self.obj.FloatVector( i ) )
                    if type( i[0] ) == str:
                        tmp.append( self.obj.StrVector( i ) )
                    if type( i[0] ) == int:
                        tmp.append( self.obj.IntVector( i ) )
                else:
                    if type( i ) == float:
                        tmp.append( self.obj.FloatVector( [i] ) )
                    if type( i ) == str:
                        tmp.append( self.obj.StrVector( [i] ) )
                    if type( i ) == int:
                        tmp.append( self.obj.IntVector( [i] ) )
        return tmp

def main():

    r = R()

    a = r.getRVectors( range(10) )
    m = r().matrix( [0, 2] , nrow = 3, ncol=2, byrow = True,
            dimnames = r.getRVectors( ["r", "theta", "tau"], ["1","second"] )  )

    r().X11()
    r().layout( r().matrix( [1] ))

    y = r().rnorm(10)
    r().plot( r().runif(10), y, xlab="runif", ylab="foo/bar", col="red")



if __name__ == '__main__':
    main()
