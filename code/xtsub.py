#!/usr/bin/python
import sys, os

def clean(string):
	return ''.join(e for e in string if e.isalnum())
	

Index = range(10);
lens  = [5,10,20,50,100,200,500,1000]
KS = [1,2]




for k in KS:
	g = open('RES-'+str(k), 'w')
	for l in lens:
		g.write(str(l)+'\t')
		for i in Index:
						
			f = open('D100/S-' + str(k) + '-' + str(l) + '-' + str(i) + '_probes-1.compare' , 'r')
			c = f.readline().split()[1]
			g.write(c+'\t')
			u = f.readline().split()[1]
			f.close()


		g.write('\n')

	g.write('\n\nuniform_cost\t'+u)
	g.close()
