#!/usr/bin/python
import sys, os

Index = range(10);
# Samples = [50, 100, 200, 500, 1000, 2000,10000]

f = open(sys.argv[1], 'r')

line = f.readline()
while (line.strip())[0] == '#':
	line = f.readline()

def clean(string):
	return ''.join(e for e in string if e.isalnum())

while line:

	if line.strip()[0] == '#':
		continue

	# for samp in Samples:
	for ind in Index:
		new_line = line.strip() + " -IND " + str(ind)
		
		x_name = ''.join([clean(s)+'-' for s in new_line.split()])
		job_name = 'GRID/jobs/j' + x_name + '.sh'
		g = open(job_name, 'w')
		g.write('#!/bin/bash' + '\n')
		g.write(new_line)
		g.close()
		os.system("chmod +x " + job_name)

		time	= 'day'
		memory	= '12G'
		host	= 'mblade12'
		

		
		command = "qsub -cwd -l %s -l vf=%s -o GRID/outs/%s -e GRID/errors/%s -q '*@@%s' %s" % (time, memory, x_name, x_name, host, job_name)
		# print command
		os.system(command)

	line = f.readline()
