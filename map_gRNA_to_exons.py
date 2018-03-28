#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys

from collections import defaultdict

class Intvl:
   def __init__(self, x, y, data=None):
      # Swap x and y if x > y.
      if x > y: (x,y) = (y,x)
      self.x = x
      self.y = y
      self.data = data

   def mid(self):
      return (self.x + self.y)/2

   def __contains__(self, z):
      return self.x < z < self.y

   def __lt__(self, z):
      return self.y < z

   def __gt__(self, z):
      return self.x > z


class IntvlNode:
   def __init__(self, intvl):
      self.mid = intvl.mid()     # The threshold for left and right.
      self.cross = set([intvl])  # A set of intervals crossing 'mid'.
      self.left = None           # Left child IntvlNode.
      self.right = None          # Right child IntvlNode.

   def addIntvl(self, intvl):
      '''Recursively find the place to add interval.'''
      if intvl < self.mid:
         # Interval is left of midpoint.
         if self.left is None: self.left = IntvlNode(intvl)
         else: self.left.addIntvl(intvl)
      elif intvl > self.mid:
         # Interval is right or midpoint.
         if self.right is None: self.right = IntvlNode(intvl)
         else: self.right.addIntvl(intvl)
      else:
         # Interval intersects midpoint.
         self.cross.add(intvl)

   def query(self, x):
      '''Recursively query the node and its descent.'''
      # Yield all intervals in cross.    
      for intvl in self.cross:
         if x in intvl: yield intvl
      # Query left.
      if x < self.mid and self.left:
         for found in self.left.query(x): yield found
      # Query right.
      elif x > self.mid and self.right:
         for found in self.right.query(x): yield found

   @staticmethod
   def create_tree(intvl_set):
      root = IntvlNode(intvl_set.pop())
      for intvl in intvl_set: root.addIntvl(intvl)
      return root


def exons_to_interval_tree(f):
   set_dict = defaultdict(set)
   # Lines to parse as the following.
   # uc001aaa.3,chr1,+,11873,12227
   for line in f:
      (gene,chr,strand,start,end) = line.rstrip().split(',')
      set_dict[chr].add(Intvl(int(start), int(end), gene))
   tree_dict = dict()
   for chr,intvl_set in set_dict.items():
      tree_dict[chr] = IntvlNode.create_tree(intvl_set)
   return tree_dict


if __name__ == '__main__':
   with open(sys.argv[1]) as f:
      tree_dict = exons_to_interval_tree(f)
   for hit in tree_dict['chr1'].query(12000): print hit.data
