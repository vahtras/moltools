#!/usr/bin/env python
#-*- coding: utf-8 -*-

import os, re, time

from multiprocessing import Queue, Process, Pool
from loprop import MolFrag, penalty_function, shift_function

mollist = [0, 0, 0, 0, 0, 0, ]

class Worker(Process):
    def __init__(self, queue):
        super(Worker, self).__init__()
        self.queue= queue
    def run(self):
        print 'Worker started'
        # do some initialization here
        print 'Computing things!'
        for data in iter( self.queue.get, None ):
            out = MolFrag( tmpdir = "/tmp/%d" %data, max_l = 0, pol = 0 , freqs = None, pf = penalty_function, sf = shift_function ).output_potential_file( 0, 0, None, False)
            print out
            mollist[data] = out

            # Use data

#pat = re.compile(r'(\d).*(\d)')

mollist = [ 0, 0, 0, 0, 0, 0 ]
q = Queue()

def f( num, ):
    global mollist
    mollist[num] = MolFrag( tmpdir = "/tmp/%d" %num, max_l = 0, pol = 0 , freqs = None, pf = penalty_function, sf = shift_function ).output_potential_file( 0, 0, None, False)

def fq(q, num):
    q.put( MolFrag( tmpdir = "/tmp/%d" %num, max_l = 1, pol = 2 , freqs = None, pf = penalty_function, sf = shift_function ).output_potential_file( 1, 0, None, False) )

#for i in [f for f in os.listdir(os.getcwd()) if f.endswith('gz') ]:
#    num = pat.search( i ).group(2)
#    if not os.path.isdir( "/tmp" + "/%s" %str(num)):
#        os.mkdir( os.path.join ( "/tmp", str(num)))

t1 = time.time()

for i in range(4):
    Worker( q ).start()

for data in range(1, 6 ):
    q.put( data )

for i in range(4):
    q.put( None ) 

t2 = time.time()

print t2 - t1

