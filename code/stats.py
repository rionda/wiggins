#!/usr/bin/python

import sys

f = open(sys.argv[1],'r')
line = f.readline()
V = set([])
E = set([])
while line:
	u, v = line.split()[0], line.split()[1]
	V.add(int(u))
	V.add(int(v))
	E.add((int(u),int(v)))
	line = f.readline();
f.close()

print "number of nodes: ", len(V)
print "number of edges: ", len(E)
deg = {v:0 for v in V}
for u,v in E:
	deg[u] += 1

V1000 = 0
V500 = 0
V100 = 0

for u in V:
	if deg[u] >= 1000:
		V1000 += 1
	elif deg[u] >= 500:
		V500 +=1
	elif deg[u] >= 100:
		V100 += 1

print "V-1000: ", V1000
print "V-500: ", V500
print "V-100: ", V100
