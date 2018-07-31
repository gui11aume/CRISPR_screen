#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys

import seeq

from gzopen import gzopen

def main(f):
   # This is the sequence right in front of the gRNA
   # It is used as oligo for the CPR and is present in all reads.
   oligo = seeq.compile("TTGTGGAAAGGACGAAACACCG", 3)
   plasmid = seeq.compile("TTTTAGAGCTAGAAAT", 2)
   for lineno,line in enumerate(f):
      if lineno % 4 != 1: continue
      try:
         # Use Seeq to extract the gRNA as suffix of the read.
         ignore, suffix = oligo.match(line.rstrip()).split()
         gRNA, ignore = plasmid.match(suffix).split()
      except (AttributeError, ValueError):
         continue
      # Remove short sequences.
      print gRNA


if __name__ == '__main__':
   with gzopen(sys.argv[1]) as f:
      main(f)
