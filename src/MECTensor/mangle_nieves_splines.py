#!/usr/bin/env python

import sys

infile=sys.argv[1]
outfile=sys.argv[2]

# Assumes the input file is one line with space-separated numbers
xs=[ float(x) for x in open(infile).readline().split() ]

outfile=open(outfile, "w")

for x in xs:
    if x==0:
        print >>outfile, 0,
    else:
        print >>outfile, "%.5g" % x,
