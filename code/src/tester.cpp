#include <iostream>
#include <vector>
#include <deque>
#include <numeric>
#include <cstring>
#include <fstream>
#include <cmath>

using namespace std;

string my_get_dir_name(const string filename) {
	string dir;
	const size_t last_slash_idx = filename.rfind('/');
	if (string::npos != last_slash_idx)
	{
    	dir = filename.substr(0, last_slash_idx) + "/";
	}
	return dir;
}


vector<double> cost(vector<vector<double>>& pp, const string test_file, const double theta, const int time_length, const int cc) {
	ifstream fin;
	fin.open(test_file);
	if (! fin.is_open()) {
		cout << "cost: file not found: " << test_file << endl;
		// return;
	}
	
	vector<double> cost_val(pp.size(),0);

	int k, u;
	double pS;
	
	
	while (fin >> k) {
		vector<double> pS(pp.size(), 0);
		for (int i=0; i<k; ++i) {			
			fin >> u;
			for (int r=0; r<pp.size(); ++r) {
				pS[r] += pp[r][u];
			}			
		}
		for (int r=0; r<pp.size(); ++r) {
			cost_val[r] += 1.0/(1 - theta*pow(1-pS[r],cc));
		}
	}

	for (int r=0; r < pp.size(); ++r) {
		//cost_val[r] /= theta;
		cost_val[r] /= time_length;	
	}

	return cost_val;
}

vector<vector<double>> load_pp(const string pp_file, const int number_of_nodes, const int num_iter) {
	vector<vector<double>> pp(num_iter);
	ifstream fin;
	fin.open(pp_file);
	if (! fin.is_open()) {
		cout << "load_pp: file not found: " << pp_file << endl;
		// return;
	}
	double p_val;
	for (int r=0; r < num_iter; ++r) {
		pp[r] = vector<double>(number_of_nodes);
		for (int i=0; i<number_of_nodes; ++i) {
			fin >> p_val;
			pp[r][i] = p_val;
		}
	}
	return pp;
}






int main(int argc, char *argv[]) {
	
	vector<clock_t> iter_time;
	int num_iter = 51;
	int number_of_nodes;
	int time_length;
	double theta = 0.75;

	string input_graph = argv[1];
	int cc = stoi(argv[2]);
	// string test_file = input_graph+"eps-0.1_S-20";
	// string test_file = input_graph+"_theta-"+to_string(theta)+"_S-5";//+to_string(i);	


	string test_file = my_get_dir_name(input_graph)+"S-5";

	ifstream fin;
	fin.open(test_file + ".stats");
	fin >> number_of_nodes >> time_length;
	fin.close();

	vector<int> tries = {1};

	ofstream fout_test_cost;

	fout_test_cost.open(my_get_dir_name(input_graph) + to_string(cc) +".test");
	for (int i : tries) {
		cout << "working on i: " << i << endl;
		// string pp_file = input_graph+"_theta-"+to_string(theta)+"_S-"+to_string(i) +".p";
		string pp_file = my_get_dir_name(input_graph)+"S-"+to_string(i) +".p";

		vector<vector<double>> pp = load_pp(pp_file, number_of_nodes, num_iter);
		cout << "loaded" << endl; cout.flush();
		vector<double> cost_val = cost(pp, test_file, theta, time_length, cc);
		cout << "should be done" << endl;cout.flush();
		for (double c : cost_val) {
			cout << c << "\t";
			fout_test_cost << c << "\t";
		}
		fout_test_cost << endl;

	}
	fout_test_cost.close();



	return 0;
}
