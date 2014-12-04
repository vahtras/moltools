#!/usr/bin/env python
#-*- coding: utf-8 -*-


class Node:
    def __init__(self):
        self.p = None
        self.l = None
        self.r = None
        self.v = None

class MaxHeap:
    def __init__(self):
        self.root = None

    def put(self, v):

        if self.root is None:
            n = Node()
            n.v = v
            self.root  = n
        else:
            if v <= self.root.v:
                if self.root.l is None:
                    self.root.l = self.insert( self.root.l, v )
                    self.root.l.p = self.root
                else:
                    if self.root.r is None:
                        self.root.r = self.insert( self.root.r, v )
                        self.root.r.p = self.root
                    else:
                        self.root.l = self.insert( self.root.l, v )
                        self.root.l.p = self.root
            else:
                n = Node()
                n.v = v
                n.l = self.root
                self.root.p = n
                self.root = n

    def insert(self, node, v ):

        if node is None:
            n = Node()
            n.v = v
            return n

        if v <= node.v:
            if node.l is None:
                node.l = self.insert( node.l, v )
                node.l.p = node
            else:
                if node.r is None:
                    node.r = self.insert( node.r, v )
                    node.r.p = node
                else:
                    node.l = self.insert( node.l, v )
                    node.l.p = node
        else:
            n = Node()
            n.v = v
            n.l = node
            node.p = n
            node = n
        return node

    def print_heap(self):
        string = ""
        self.printVals( self.root, string )

    def printVals(self, node, string ):
        """Prints vues of nodes recursively"""
        if node is None:
            return
        print string + str(node.v )
        string += " - " 
        self.printVals( node.l, string )
        self.printVals( node.r, string )

if __name__ == '__main__':
    h = MaxHeap()

