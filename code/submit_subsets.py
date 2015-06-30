#!/usr/bin/python
import sys, os

def clean(string):
	return ''.join(e for e in string if e.isalnum())
	

Index = range(10);
lens  = [1,2,5,10,20,50,100,200,500,1000,2000]
KS = [1,2,3]





for l in lens:
	for i in Index:
		for k in KS:
			new_line = './phsp -TASK solver-with-file -PROBES 1 -SAMPLE_FILE D100/S-' + str(k) + '-' + str(l) + '-' + str(i)

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



