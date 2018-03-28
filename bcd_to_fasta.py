#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys

def main(f):
   for lineno,line in enumerate(f):
      bcd, cnt = line.split()
      # Do not use gRNAs with less than 5 reads.
      if int(cnt) < 5: break
      print '>%d_%s\n%s' % (lineno+1, cnt, bcd)

if __name__ == '__main__':
   with open(sys.argv[1]) as f:
      main(f)

