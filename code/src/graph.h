#include <iostream>
#include <vector>
#include <deque>
#include <cstring>
#include <unordered_set>
#include <fstream>
#include <algorithm>
#include <tuple>
#include <set>

using namespace std;
typedef vector<vector<deque<int>>> prev_t;
typedef vector<vector<int>>	vvint_t;
typedef vector<vector<double>>	vvdoub_t;
typedef vector<deque<int>>	vdint_t;
typedef vector<double>		vdoub_t;
typedef deque<int>			deint_t;
typedef unordered_set<int> 	unset_t;


class Graph
{
private:
	vvint_t nodes; // here we assume nodes are 0,1,...,N-1		
	string graph_name;
	uniform_int_distribution<> distr;
	vector<vector<bernoulli_distribution>> bern;
	set<int> random_spot(int source);
		
public:
	int number_of_nodes;
	int number_of_edges;	
	double theta;
	Graph(const string input_file, const double theta);
	// Graph(const string input_file);
	/*
		input_file: each row is an edge with probability.
		graph is "always directed"
		each edge (a,b): a --> b : a influences b.
		influence_type: for edge (a,b) is 0 if b influences a, otherwise it is 1.
	*/
	void gen_sample(const int time_interval, const string output_file);
	void gen_sample_with_time(const int time_interval, const string output_file, const int time_idx);
	void impose(vvint_t& schedule, set<int>& spot, int current_time, vector<double>& cost);
	void simulate(vvint_t & schedule, vector<double>& cost);
	vector<double> simul_with_file(const int time_interval, const string sample_file, const int num_iter, const int time_idx);
	vector<double> solver(const string sample_file, const int time_length, const int num_iter);
	vdoub_t simul_process(double eps, int num_iter, string idxstr);
	void go_online(vvint_t & schedule, vector<double>& cost, const int time_idx, const int time_interval);
	void go_offline(vvint_t & schedule, vector<double>& cost, const string sample_file);
	void inf_gen_sample(const int num_samples, const string output_file);
};


