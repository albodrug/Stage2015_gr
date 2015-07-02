#!/usr/bin/python

import os
import os.path
import argparse
import re
from collections import Counter


parser = argparse.ArgumentParser(description='Parses the pacanier file and transforms it into a nice gurobi input file ')
parser.add_argument('-f', help='The file')

args = parser.parse_args()

# FILE
input_file = args.f


# EXTRACTING
f = open(input_file, 'r')
relevant_lines = []

for line in f:
  if 'len' in line:
    relevant_lines.append(line)

f.close()

contigs_and_coverage = []
links_between_contigs = []

for element in relevant_lines:
  if 'R' not in element and 'F' not in element:
    arrline = element.split()
    cnr, t, clen = arrline[0].split('__')[0], arrline[0].split('__')[1], arrline[0].split('__')[2]
    ccovmin, ccovmax = arrline[1], arrline[2]

    contigs_and_coverage.append([arrline[0], clen, arrline[1], arrline[2]])
  if ('R' in element or 'F' in element) and len(element.split()) == 3:
    c1 = element.split()[0]
    c2 = element.split()[1]
    ovlp = element.split()[2]
    cnr1, t, clen1, cord1 = c1.split('__')[0], c1.split('__')[1], c1.split('__')[2][:-2], c1.split('_')[-1]
    cnr2, t, clen2, cord2 = c2.split('__')[0], c2.split('__')[1], c2.split('__')[2][:-2], c2.split('_')[-1]
    
    links_between_contigs.append([cnr1+'__len__'+clen1+'__'+cord1, cnr2+'__len__'+clen2+'__'+cord2, '-'+ovlp])
  if ('R' in element or 'F' in element) and len(element.split()) == 4:
    c1 = element.split()[0]
    c2 = element.split()[1]
    dmin = element.split()[2]
    dmax = element.split()[3]
    cnr1, t, clen1, cord1 = c1.split('__')[0], c1.split('__')[1], c1.split('__')[2][:-2], c1.split('_')[-1]
    cnr2, t, clen2, cord2 = c2.split('__')[0], c2.split('__')[1], c2.split('__')[2][:-2], c2.split('_')[-1]

    links_between_contigs.append([cnr1+'__len__'+clen1+'__'+cord1, cnr2+'__len__'+clen2+'__'+cord2, str((int(dmin)+int(dmax))/2)])


# OUTPUTTING

print '\n'
print 'Contigs and coverages \n'
for element in contigs_and_coverage:
  print element[0], element[1], element[2], element[3]
print '\n'
print 'Links between contigs \n'
for element in links_between_contigs:
  print element[0], element[1], element[2]
print '\n'











