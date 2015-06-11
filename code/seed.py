#!/usr/bin/python

import sys

f = open(sys.argv[1],'r')
line = f.readline()
f.close()
p = [float(st) for st in line.split()]
x = [(p[i],i) for i in xrange(len(p))]
x.sort()
x.reverse()

seed = [t[1] for t in x]

out=sys.argv[1].split('.')[0] + '.seed'
g = open(out, 'w')
for s in seed:
	g.write(str(s) + '\t')
g.close()
