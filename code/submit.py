#!/usr/bin/python
import sys, os

f = open(sys.argv[1], 'r')

line = f.readline()
while (line.strip())[0] == '#':
	line = f.readline()

def clean(string):
	return ''.join(e for e in string if e.isalnum())

while line:

	if line.strip()[0] == '#':
		continue

	x_name = ''.join([clean(s)+'-' for s in line.split()])
	job_name = 'GRID/jobs/j' + x_name + '.sh'
	g = open(job_name, 'w')
	g.write('#!/bin/bash' + '\n')
	g.write(line)
	g.close()
	os.system("chmod +x " + job_name)

	time	= 'inf'
	memory	= '12G'
	host	= 'mblade12'
	

# qsub -cwd -l inf -l vf=32G -o outs/out%i -e errors/err%i -m abes -q '*@@mblade12' jobs/j%i.sh \n""" % (i, i, i))
	command = "qsub -cwd -l %s -l vf=%s -o GRID/outs/%s -e GRID/errors/%s -m abes -q '*@@%s' %s" % (time, memory, x_name, x_name, host, job_name)
	#print command
	os.system(command)

	line = f.readline()
