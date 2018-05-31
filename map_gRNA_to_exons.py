#!/usr/bin/env python
# -*- coding:utf-8 -*-

import random
import sys

from collections import defaultdict
from math import log

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
      (gene,chrom,strand,start,end) = line.rstrip().split(',')
      set_dict[chrom].add(Intvl(int(start), int(end), gene))
   dTree = dict()
   for chrom,intvl_set in set_dict.items():
      dTree[chrom] = IntvlNode.create_tree(intvl_set)
   return dTree


def map_reads_to_exons(f, dTree):
   dScore = defaultdict(float)
   dNread = dict()
   for line in f:
      # File is assumed to be in sam format.
      items = line.split()
      nreads = int(items[0].split('_')[-1])
      chrom = items[2]
      pos = int(items[3])
      dNread[(chrom,pos)] = nreads
      for hit in dTree[chrom].query(pos):
         dScore[hit.data] += log(nreads)
   # Perform 100 bootstraps.
   dBootstrap = defaultdict(lambda: [0.0] * 100)
   for (chrom,pos) in dNread.keys():
      genes = [hit.data for hit in dTree[chrom].query(pos)]
      if not genes: continue
      for bootn in range(100):
         # Sample a random number of reads.
         nreads = random.choice(dNread.values())
         for gene in genes:
            dBootstrap[gene][bootn] += log(nreads)
   # Compute bootstrap average and extremes.
   dBootav = dict()
   dBoothi = dict()
   for gene in dBootstrap:
      dBootav[gene] = sum(dBootstrap[gene]) / 100
      dBoothi[gene] = max(dBootstrap[gene])
   genes_in_score_order = sorted(dScore, key=dScore.get, reverse=True)
   genes_in_boot_order = sorted(dBootav, key=dBootav.get, reverse=True)
   # Print scores (and rank-matched bootstrap scores).
   for rank in range(len(genes_in_score_order)):
      gene  = genes_in_score_order[rank]
      score = dScore[gene]
      boot  = dBootav[genes_in_boot_order[rank]]
      hi    = dBoothi[genes_in_boot_order[rank]]
      print gene, score, boot, hi


if __name__ == '__main__':
   # Convert exon locations to interval tress for fast search.
   with open(sys.argv[1]) as f:
      dTree = exons_to_interval_tree(f)
   # Intersect read positions with exons.
   with open(sys.argv[2]) as f:
      map_reads_to_exons(f, dTree)
