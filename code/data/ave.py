#!/usr/bin/python
import sys

input_file = sys.argv[1]
f = open(input_file)
line = f.readline()
s = 0.0
c = 0
while line:
	c += 1
	s += int(line.split()[0])
	line = f.readline()
f.close()
print s/c
