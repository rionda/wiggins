#include <cmath>
#include "graph.h"

string my_get_dir_name(const string filename) {
	string dir;
	const size_t last_slash_idx = filename.rfind('/');
	if (string::npos != last_slash_idx)
	{
    	dir = filename.substr(0, last_slash_idx) + "/";
	}
	return dir;
}

int main(int argc, char *argv[]) {
	
	string input_file = argv[1];
	string idx_st = argv[2];
	// double eps = 0.1;
	double eps = 1;
	// int num_iter = 50;
	int num_iter = 15;
	// double theta = stod(argv[3]);
	double theta = 0.75;

	Graph g(input_file, theta);
	vdoub_t cost = g.simul_process(eps, num_iter,""+idx_st);
	ofstream fout;
	fout.open(my_get_dir_name(input_file) + "simul-"+idx_st);
	for (double c : cost)
		fout << c << "\t";
	fout << endl;
	fout.close();
	
	

	



	// cout << "number of nodes: " << g.number_of_nodes << " " << g.number_of_edges << endl;
	// g.gen_sample(5, string(argv[1])+".samp");
	return 0;
}
