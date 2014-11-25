#!/usr/bin/python
# Given an input BED file, will subsample the lines down.
# arguments are infile N1 N2 N3...
# where each Ni is a number (in millions) of lines that we want in a subsampled file...
# or a fraction of original output size!

import random
import gzip
import sys
import re

infile_path = sys.argv[1]
if not '.bed' in infile_path:
    print 'ERROR: expecting a bed file.'
    sys.exit()

# first, check what sort of input
num_list = map(float,sys.argv[2:])
subsamps = []
# assume if all are <=1, it's a list of fractions
# else, assume it's a list of read-depths in millions
for num in num_list:
    if num > 1:
        millions = True
        break
    millions = False

num_lines = sum(1 for line in open(infile_path))
print 'There are', num_lines,'lines in', infile_path+'!'

if millions == True:
    print 'Interpreting input as millions.'
    for num in num_list:
        prob = num*(1000000.0)/num_lines
        path = re.sub('.bed', '_'+str(int(num))+'m.bed', infile_path)
        if not '.gz' in infile_path:
            path = path+'.gz'
        outfile = gzip.open(path, 'w')
        subsamps.append((prob, outfile))
else:
    print 'Interpreting input as fractions.'
    for num in num_list:
        prob = num
        path = re.sub('.bed', '_'+str(int(num*100))+'pc.bed', infile_path)
        if not '.gz' in infile_path:
            path = path+'.gz'
        outfile = gzip.open(path, 'w')
        subsamps.append((prob, outfile))

# prepare to run through!
i=0
print 'Opening file', infile_path
if '.gz' in infile_path: f = gzip.open(infile_path)
else: f = open(infile_path)
for line in f:
    if i%5000000==0:
        print str(100*float(i)/num_lines)+'%'
    i+=1
    if 'rRNA' in line:
        # in theory we can exclude anything else here, too
        continue
    for option in subsamps:
        r = random.random()
        prob = option[0]
        outfile = option[1]
        if r <= prob:
            outfile.write(line)
