#!/usr/bin/env python

from rootpy.tree import Tree, TreeModel
from rootpy.io import open
from rootpy.types import *
from random import gauss

f = open("test.root", "recreate")


# define the model
class Event(TreeModel):

    x = Float()
    y = Float()
    z = Float()
    i = Int()

tree = Tree("test", model=Event)

# fill the tree
for i in xrange(10000):
    tree.x = gauss(.5, 1.)
    tree.y = gauss(.3, 2.)
    tree.z = gauss(13., 42.)
    tree.i = i
    tree.fill()
tree.write()

f.close()
