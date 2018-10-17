#!/usr/bin/env python
# -*- coding:utf-8 -*-

import random
import re
import sys

from collections import defaultdict
from math import log

from gzopen import gzopen

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


def get_pair(stringlist):
   tx = gname = None
   for string in stringlist:
      if string.startswith('transcript_id'):
         tx = re.search(r'ENST\d+', string).group()
         continue
      if string.startswith('gene_name'):
         (gname,) = re.search(r'"([^"]+)"', string).groups()
         continue
   return tx, gname


def read_dict(f):
   cast = dict()
   for line in f:
      if line[0] == '#':
         continue
      items = line.split('\t')
      annot = re.split(r';\s*', items[8])
      try:
         tx, gname = get_pair(annot)
      except AttributeError:
         continue
      if tx is not None and gname is not None:
         cast[tx] = gname
   return cast


def exons_to_interval_tree(f, cast):
   set_dict = defaultdict(set)
   # Lines to parse as the following.
   # uc001aaa.3,chr1,+,11873,12227
   for line in f:
      (tx,chrom,strand,start,end) = line.rstrip().split(',')
      tx = re.sub(r'\.\d+$', '', tx)
      if tx not in cast: continue
      gene = cast[tx]
      set_dict[chrom].add(Intvl(int(start), int(end), gene))
   dTree = dict()
   for chrom,intvl_set in set_dict.items():
      dTree[chrom] = IntvlNode.create_tree(intvl_set)
   return dTree


def map_reads_to_exons(f, dTree):
   sys.stdout.write('sgRNA\tgene\tcount\n')
   for line in f:
      # File is assumed to be in sam format.
      items = line.split()
      sgRNA = items[9]
      nreads = items[0].split('_')[-1]
      chrom = items[2]
      pos = int(items[3])
      if chrom not in dTree:
         # Do not consider alternate chromosomes
         # do not add them to 'dNread'.
         continue
      seen = set()
      for hit in dTree[chrom].query(pos):
         if hit.data in seen: continue
         # Every hit is a gene intersected by the position
         # of the gRNA (more specifically one exon of the
         # gene is intersected by the position of the gRNA).
         sys.stdout.write('%s\t%s\t%s\n' % (sgRNA, hit.data, nreads))
         seen.add(hit.data)


if __name__ == '__main__':
   sys.setrecursionlimit(10000)
   # Prepare transcript-gene lookup
   with gzopen(sys.argv[1]) as f:
      cast = read_dict(f)
   # Convert exon locations to interval tress for fast search.
   with gzopen(sys.argv[2]) as f:
      dTree = exons_to_interval_tree(f, cast)
   # Intersect read positions with exons.
   with gzopen(sys.argv[3]) as f:
      map_reads_to_exons(f, dTree)
