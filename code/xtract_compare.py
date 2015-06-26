#!/usr/bin/python
import sys, os

Index = range(10);
Samples = [50, 100, 200, 500, 1000, 2000,3000,4000,5000,6000,7000,8000,9000]


dataset = sys.argv[1]
g = open(dataset+"COMP", "w")
for samp in Samples:
	g.write(str(samp)+"\t")
	for ind in Index:
		f = open(dataset + "indexed/S-"+str(samp)+"-"+str(ind)+"_probes-1.compare")
		cost = f.readline().split()[1]
		f.close()
		g.write(cost + "\t")
	g.write("\n")

g.close()


g = open(dataset+"COMP-unif", "w")
for samp in Samples:
	g.write(str(samp)+"\t")
	for ind in Index:
		f = open(dataset + "indexed/S-"+str(samp)+"-"+str(ind)+"_probes-1.compare")
		tmp = f.readline();
		cost = f.readline().split()[1]
		f.close()
		g.write(cost + "\t")
	g.write("\n")

g.close()

g = open(dataset+"COMP-out", "w")
for samp in Samples:
	g.write(str(samp)+"\t")
	for ind in Index:
		f = open(dataset + "indexed/S-"+str(samp)+"-"+str(ind)+"_probes-1.compare")
		tmp = f.readline();
		tmp = f.readline();
		cost = f.readline().split()[1]
		f.close()
		g.write(cost + "\t")
	g.write("\n")

g.close()



g = open(dataset+"COMP-in", "w")
for samp in Samples:
	g.write(str(samp)+"\t")
	for ind in Index:
		f = open(dataset + "indexed/S-"+str(samp)+"-"+str(ind)+"_probes-1.compare")
		tmp = f.readline();
		tmp = f.readline();
		tmp = f.readline();
		cost = f.readline().split()[1]
		f.close()
		g.write(cost + "\t")
	g.write("\n")

g.close()
