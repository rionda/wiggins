import sys

if __name__ == "__main__":

	input_graph = sys.argv[1]
	f = open(input_graph,'r')
	g = open(input_graph+'.switched','w')
	line = f.readline()
	while line:
		str_u, str_v = line.split()[0], line.split()[1]
		tmp = '\t'.join(line.split()[2:])
		g.write(str_v + '\t' + str_u + '\t' + tmp + '\n')
		line = f.readline();

	f.close()
	g.close()
