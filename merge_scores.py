#!/usr/bin/env python
# -*- coding:utf-8 -*-

import re
import sys

from collections import defaultdict

# From https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#Python
def levenshtein_dist_less_than(s1, s2, threshold):
   if len(s1) < len(s2):
      return levenshtein_dist_less_than(s2, s1, threshold)

   # len(s1) >= len(s2)
   if len(s2) == 0:
      return len(s1)

   previous_row = range(len(s2) + 1)
   for i, c1 in enumerate(s1):
      current_row = [i + 1]
      for j, c2 in enumerate(s2):
         insertions = previous_row[j + 1] + 1
         deletions = current_row[j] + 1
         substitutions = previous_row[j] + (c1 != c2)
         current_row.append(min(insertions, deletions, substitutions))
      previous_row = current_row
      if previous_row[min(i, len(s2))] > threshold: return False

   return True



def normalize(seq, group_of_seq):
   for s in group_of_seq:
      if levenshtein_dist_less_than(seq, s, 5): return s
   return seq


def read_file(fname, gene_dict):
   with open(fname) as f:
      # Remove header.
      discard = next(f)
      for line in f:
         sgRNA, gene, count = line.split()
         if gene in gene_dict:
            sgRNA = normalize(sgRNA, gene_dict[gene].keys())
         gene_dict[gene][sgRNA][fname] = int(count)


def main():
   fname_list = sys.argv[1:]
   # Read in the files.
   gene_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
   for fname in fname_list:
      read_file(fname, gene_dict)

   # Output.
   genes = sorted(gene_dict.keys())
   header = 'sgRNA\tgene\t' + \
         '\t'.join(re.sub(r'\..*', '', fname) for fname in fname_list)
   sys.stdout.write(header + '\n')

   for gene in genes:
      for sgRNA in gene_dict[gene].keys():
         local_dict = gene_dict[gene][sgRNA]
         line = ('%s\t%s\t' % (sgRNA, gene)) + \
            '\t'.join(['%d' % local_dict[x] for x in fname_list])
         sys.stdout.write(line + '\n')

if __name__ == '__main__':
   main()
