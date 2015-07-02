#!/usr/bin/python

# This script takes the megarray of solution obtained from genscale scaffolding tools
# and converts it to a megascaffold sequence thanks to the contig fasta file

import argparse
import ast

parser = argparse.ArgumentParser(description='Construct of the mega-scaffold from the genscaf \
                                            output (ordered list of contigs).')
parser.add_argument('-c', help='Original contig file, fasta format')
parser.add_argument('-s', help='Genscaf array_of_multiple_solutions file')

args = parser.parse_args()

# FILES
contig_file = args.c
megarray_file = args.s

# FUNCTIONS
def extract_cnr(c_name):
  cnr = c_name.split('_')[0]
  return cnr
def contruct_megascaf_name(element1):
  name=''
  for c in element1:
    name += c.split('_')[0]+c.split('_')[1]+'_'
  return name[:-1]
def construct_megascaf_seq(ordered_contigs):
  megascaffold_seq = ''
  for cnr in ordered_contigs:
    megascaffold_seq += dict_nrc_2_cseq[cnr].rstrip()
  return megascaffold_seq

# CONSTRUCTING THE SEQ
dict_nrc_2_cseq = {}

c = open(contig_file, 'r')
while True:
    c_name = c.readline()
    cnr = (c_name.split('__')[0]).replace('>','')
    cseq = c.readline()
    dict_nrc_2_cseq[cnr] = cseq
    if not cnr: 
      break
c.close()

s = open(megarray_file, 'r')
str_array_of_multiple_solutions = s.readline()
s.close()
array_of_multiple_solutions = ast.literal_eval(str_array_of_multiple_solutions)

for element in array_of_multiple_solutions:
  ordered_contigs = element[1]
  ordered_contigs = [extract_cnr(c) for c in ordered_contigs]
  megascaffold_name = contruct_megascaf_name(element[1])
  megascaffold_seq = construct_megascaf_seq(ordered_contigs)
  print '>Scaffold_megagenscaf_'+megascaffold_name
  print megascaffold_seq

# Multifasta split it