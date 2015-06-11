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


// void update(vector<double>& p, vector<double>& W) {
// 	for (int i=0; i < p.size(); ++i) {
// 		p[i] *= sqrt(W[i]);
// 	}
// 	double sum = accumulate(p.begin(), p.end(), 0.0);
// 	for (int i=0; i < p.size(); ++i) {
// 		p[i] /= sum;
// 	}
// }

void update(vector<double>& p, vector<double>& W) {
	for (int i=0; i < p.size(); ++i) {
		p[i] *= W[i];
	}
	double sum = accumulate(p.begin(), p.end(), 0.0);
	for (int i=0; i < p.size(); ++i) {
		p[i] /= sum;
	}
}

void solver(string sample_file, int number_of_nodes, double theta, int time_length, 
	int num_iter, vector<double>& cost, vector<clock_t>& iter_time, const int cc) {
	
	vector<double> p(number_of_nodes, 1.0/number_of_nodes);	

	cost = vector<double>(num_iter, 0);
	iter_time = vector<clock_t>(num_iter);
	
	ifstream fin;	
	int k, u;
	double pS, tmp;
	clock_t t;

	ofstream fout_p;
	fout_p.open(sample_file+".p");

	for (int r=0; r < num_iter; ++r) {		
		cout << r << " "; cout.flush();
		// saving the schedule:
		for (double d : p) {
			fout_p << d << "\t";
		}
		fout_p << endl;


		t = clock();
		fin.open(sample_file);
		vector<double> W(number_of_nodes, 0);		
		while (fin >> k) {
			pS = 0;
			deque<int> S; // TODO: check if it gets new everytime
			for (int i=0; i<k; ++i) {
				fin >> u;
				pS += p[u];
				S.push_back(u);
			}
			tmp = 1.0/(1 - theta*pow(1-pS,cc));
			cost[r] += tmp;

			// update W:
			tmp = pow(tmp,2)*pow(1-pS,cc-1);
			for (int v : S) {
				W[v] += tmp;
			}
		}
		fin.close();
		// cost[r] /= theta;
		cost[r] /= time_length;
		update(p, W);
		iter_time[r] = clock()-t;
	}
	fout_p.close();
}


int main(int argc, char *argv[]) {

	vector<double> cost;
	vector<clock_t> iter_time;
	int num_iter = 51;
	int number_of_nodes;
	int time_length;
	

	double theta = 0.75;
	

	string input_graph = argv[1];
	int cc = stoi(argv[2]);
	string sample_file;

	
	vector<int> tries = {1};

	
	ofstream fout_conv;
	
	fout_conv.open(my_get_dir_name(input_graph) + to_string(cc) + ".conv");

	for (int i : tries) {
		cout << "\n\nworking on i: " << i << endl; cout.flush();
		cout << "\nrounds: \n";

		sample_file = my_get_dir_name(input_graph)+"S-"+to_string(i);

		ifstream fin;
		fin.open(sample_file + ".stats");
		fin >> number_of_nodes >> time_length;
		fin.close();		

		solver(sample_file, number_of_nodes, theta, time_length, num_iter, cost, iter_time, cc);

		for (double c : cost)
			fout_conv << c << "\t";
		fout_conv << endl;

		for (clock_t t : iter_time)
			fout_conv << t << "\t";
		fout_conv << endl;
	}
	fout_conv.close();
	cout << endl;
	return 0;
}
