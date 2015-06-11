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

void solver(string sample_file, int number_of_nodes, double theta, int time_length, int num_iter, vector<double>& cost, vector<clock_t>& iter_time) {
	
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
			tmp = 1.0/(1 - theta*(1-pS));
			cost[r] += tmp;

			// update W:
			tmp = pow(tmp,2);
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
	int num_iter = 15;
	int number_of_nodes;
	int time_length;
	clock_t t1, t2;
	// double theta = stod(argv[2]);
	// double theta = 0.75;
	double theta = 0.0001;

	string input_graph = argv[1];
	string sample_file;

	// vector<int> tries = {1,2,3,4,5,20};
	vector<int> tries = {1};

	
	ofstream fout_conv;
	// fout_cost.open(input_graph + "_theta-" + to_string(theta)+ "_results.cost");
	// fout_time.open(input_graph + "_theta-" + to_string(theta)+ "_results.time");
	
	fout_conv.open(my_get_dir_name(input_graph) +  "conv");
	// fout_time.open(my_get_dir_name(input_graph) +  "conv.time");

	ofstream fout_time_inf;
	fout_time_inf.open(my_get_dir_name(input_graph) +  "inf.time");
	for (int i : tries) {
		// sample_file = input_graph+"eps-0.1_S-"+to_string(i);
		// sample_file = input_graph+"_theta-"+to_string(theta)+"_S-"+to_string(i);
		sample_file = my_get_dir_name(input_graph)+"S-"+to_string(i);

		ifstream fin;
		fin.open(sample_file + ".stats");
		fin >> number_of_nodes >> time_length >> t1;
		fin.close();		

		t2 = clock();
		solver(sample_file, number_of_nodes, theta, time_length, num_iter, cost, iter_time);
		

		fout_time_inf << clock() - t2 + t1 ;
		// fout_cost << i << "\t";
		for (double c : cost)
			fout_conv << c << "\t";
		fout_conv << endl;

		// fout_time << i << "\t";
		for (clock_t t : iter_time)
			fout_conv << t << "\t";
		fout_conv << endl;
	}
	fout_conv.close();
	fout_time_inf.close();
	

	cout << endl;




	// cout << endl;
	// cout << "costs: ";
	// for (double c : cost)
	// 	cout << c << " " ;
	// cout << endl << endl;
	// cout << "time: ";
	// for (clock_t t : iter_time)
	// 	cout << t << " ";
	// cout << endl;

	
	
	// cout << "costs:\n";
	// for (double c : cost)
	// 	cout << c << "\t" ;
	// cout << endl << endl;
	// cout << "time:\n";
	// for (clock_t t : iter_time)
	// 	cout << t << "\t";
	// // cout << endl;
	return 0;
}