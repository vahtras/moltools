#!/usr/bin/env python
#-*- coding: utf-8 -*-

class Node:
    def __init__(self):
        self.n = None
        self.p = None
        self.v = None

class Circularlinkedlist:
    def __init__(self):
        self.first = None
        self.first = None
        self.last = None

    def exists(self, v ):
        if self.first is None:
            return False
        else:
            return self.existsTrue( self.first, v )

    def existsTrue(self, node, v ):
        if node is self.last:
            if node.v == v:
                return True
            else:
                return False
        else:
            return self.existsTrue( node.p, v )

    def remove(self, v):
        if self.first is None:
            return
        else:
            if self.first.v == v:
                self.first.n.p = self.first.p
                self.first.p.n = self.first.n
                self.first = self.first.p
            else:
                self.removeTrue(self.first, v)

    def removeTrue(self, n, v ):

        if n is self.last:
            if n.v == v:
                n.n.p = self.first
                self.last = n.n
            return
        else:
            if n.v == v:
                n.n.p = n.p
                n.p.n = n.n
                return
            else:
                self.removeTrue( n.p, v )
    
    def put(self, v):
        if self.first is None:
            n = Node()
            n.v = v
            self.first = n
            self.last = n
        else:
            n = Node()
            n.v = v

            self.first.n = n
            self.last.p = n

            n.p = self.first
            n.n = self.last
            self.first = n
            self.first = n

    def write(self):
        if self.first is None:
            return
        else:
            self.writeTrue( self.first )
    def writeTrue( self, node ):
        if node == self.last:
            print node.v
            return
        print node.v
        self.writeTrue( node.p )
            
if __name__ == '__main__':
    C = Circularlinkedlist()

