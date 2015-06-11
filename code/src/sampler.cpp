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
	ofstream fout;
	string input_file = argv[1];
	string dir = my_get_dir_name(input_file);
	// double theta = stod(argv[2]);
	double theta = 0.75;

	Graph g(input_file, theta);
	double eps = 0.1;
	// int r0 = ceil(log(g.number_of_nodes) * pow(theta,2) /(2*pow(eps,2)*pow(1-theta,2)));
	// int r0 = ceil(3*(log(g.number_of_nodes)+log(2))/ ((1-theta)*pow(eps,2)));
	int r0 = 2000;

	vector<int> tries = {1,5};
	for (int i : tries) {
		cout << "workin on i: " << i << endl; cout.flush();
		// g.gen_sample(i*r0, input_file+"_theta-"+to_string(theta)+"_S-"+to_string(i));
		g.gen_sample(i*r0, dir+"S-"+to_string(i));

		// fout.open(input_file+"_theta-"+to_string(theta)+"_S-"+to_string(i)+".stats");
		fout.open(dir+"S-"+to_string(i)+".stats");

		fout << g.number_of_nodes << endl;
		fout << i*r0;
		fout.close();
	}

	// cout << "workin on i: 20" << endl; cout.flush();
	// g.gen_sample(20*r0, input_file+"eps-0.1_S-20");

	
	// fout.open(input_file+"eps-0.1_S-20.stats");
	// fout <<  g.number_of_nodes << endl;
	// fout <<  20*r0;
	// fout.close();


	// cout << "number of nodes: " << g.number_of_nodes << " " << g.number_of_edges << endl;
	// g.gen_sample(5, string(argv[1])+".samp");
	return 0;
}
