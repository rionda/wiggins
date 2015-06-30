#include <cmath>
#include "graph.h"


string my_get_dir_name(const string filename) {
	string dir;
	const size_t LSI = filename.rfind('/');
	if (string::npos != LSI)
	{
    	dir = filename.substr(0, LSI) + "/";
	}
	return dir;
}

// ************************************************
//	Functins for Sampler task
// ************************************************

void sampler(Graph g, const string dir, const uint64_t len_samp) {
	/*
		len_samp: the lenght of the time interval of samples
	*/
	string samp_file = dir+"samples/S-"+to_string(len_samp);
	g.gen_sample(len_samp, samp_file);
	
	ofstream fout(samp_file + ".stat");	
	fout << g.number_of_nodes << endl;
	fout << len_samp;
	fout.close();		
}

void indexed_sampler(Graph g, const string dir, const uint64_t len_samp, const int ind) {
	/*
		len_samp: the lenght of the time interval of samples
	*/
	string samp_file = dir+"indexed/S-"+to_string(len_samp)+"-"+to_string(ind);
	g.gen_sample(len_samp, samp_file);
	
	ofstream fout(samp_file + ".stat");	
	fout << g.number_of_nodes << endl;
	fout << len_samp;
	fout.close();		
}

// ************************************************
//	Functins for Solver task
// ************************************************

void update(vector<double>& p, vector<double>& W) {
	for (int i=0; i < p.size(); ++i) {
		p[i] *= W[i];
	}
	double sum = accumulate(p.begin(), p.end(), 0.0);
	for (int i=0; i < p.size(); ++i) {
		p[i] /= sum;
	}
}

void solver(const string samp_file, 
			const double theta, 
			const int num_iter, const int probes) {

	int number_of_nodes, len_samp;
	ifstream fin(samp_file + ".stat");
	
	fin >> number_of_nodes >> len_samp;	
	fin.close();

	vector<double> p(number_of_nodes, 1.0/number_of_nodes);	//uniform schedule	
	vector<double> cost(num_iter, 0);
	vector<clock_t> iteration_time(num_iter);
	

	int k, u;
	double pS, tmp;
	ofstream fout_p(samp_file+ "_probes-" + to_string(probes) + ".p");	
	for (int r=0; r < num_iter; ++r) {				
		// saving the schedule:
		for (double d : p) {
			fout_p << d << "\t";
		}
		fout_p << endl;		
		clock_t round_time = clock();
		fin.open(samp_file);
		vector<double> W(number_of_nodes, 0);			
		while (fin >> k) {
			pS = 0;
			deque<int> S; // TODO: check if it gets new everytime
			for (int i=0; i<k; ++i) {
				fin >> u;				
				pS += p[u];
				S.push_back(u);
			}			
			tmp = 1.0/(1 - theta*pow(1-pS,probes));
			cost[r] += tmp;

			// update W:
			tmp = pow(tmp,2)*pow(1-pS,probes-1);
			for (int v : S) {
				W[v] += tmp;
			}
		}		
		fin.close();		
		cost[r] /= len_samp;
		update(p, W);
		iteration_time[r] = clock()-round_time; 
	}
	fout_p.close();
	
	ofstream fout_conv(samp_file + "_probes-" + to_string(probes) + ".conv");
	for (double c : cost)
		fout_conv << c << "\t";
	fout_conv << endl;

	for (clock_t t : iteration_time)
		fout_conv << t << "\t";
	fout_conv << endl;
	fout_conv.close();

}


// ************************************************
//	Functins for Tester task
// ************************************************

vvdoub_t load_pp(const string pp_file, const int number_of_nodes, 
							   const int num_iter) {
	vvdoub_t pp(num_iter);
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

void tester(const string samp_file, const string test_file,
			const int num_iter, const int probes,
			const double theta) {

	int number_of_nodes, len_samp;
	ifstream fin(samp_file + ".stat");	
	fin >> number_of_nodes >> len_samp;	
	fin.close();

	string pp_file = samp_file + "_probes-" + to_string(probes) + ".p";
	vvdoub_t pp = load_pp(pp_file, number_of_nodes, num_iter);

	fin.open(test_file);
	if (! fin.is_open()) {
		cerr << "cost: file not found: " << test_file << endl;		
	}
	

	vdoub_t cost_val(pp.size(),0);	
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
			cost_val[r] += 1.0/(1 - theta*pow(1-pS[r],probes));
		}
	}

	for (int r=0; r < pp.size(); ++r) {		
		cost_val[r] /= len_samp;	
	}

	
	ofstream fout_test(samp_file + "_probes-" + to_string(probes) +  ".test");
	for (double c : cost_val) {			
			fout_test << c << "\t";
	}
	fout_test << endl;
}






// ************************************************
// Functions for Simulator
// ************************************************
void simulator(const string filename, const int num_simuls, const int num_iter, const double theta, const double epsilon) {
	/*
		In PKDD submission, we let epsilon = 1, theta=0.75, and num_iter = 15;
	*/

	for (int idx=0; idx < num_simuls; ++idx) {
		Graph g(filename, theta);
		vdoub_t cost = g.simul_process(epsilon, num_iter,""+to_string(idx));
		ofstream fout;
		fout.open(my_get_dir_name(filename) + "simul-"+to_string(idx));
		for (double c : cost)
			fout << c << "\t";
		fout << endl;
		fout.close();
	}
}

// ************************************************
// Functions for Comparison
// ************************************************

vdoub_t load_p(const string pp_file, const int number_of_nodes, 
							   const int num_iter) {
	
	ifstream fin;
	fin.open(pp_file);
	if (! fin.is_open()) {
		cerr << "load_pp: file not found: " << pp_file << endl;
	} else {
		cout << "opened: " <<pp_file << endl;
	}
	double p_val;
	vdoub_t p(number_of_nodes);
	for (int r=0; r < num_iter-1; ++r) {
		for (int i=0; i<number_of_nodes; ++i) {
			fin >> p_val;
		}
	}
	for (int i=0; i<number_of_nodes; ++i) {
		fin >> p_val;
		p[i] = p_val;
	}

	cout << "number of nodes : " << p.size() << endl;
	return p;
}



void compare(const string filename, const string samp_file, const string test_file,
			const int num_iter, const int probes,
			const double theta) {

	Graph g(filename, theta);

	int number_of_nodes, len_samp;
	// ifstream fin(samp_file + ".stat");
	ifstream fin(test_file + ".stat");	
	fin >> number_of_nodes >> len_samp;	
	fin.close();

	string pp_file = samp_file + "_probes-" + to_string(probes)+ ".p";
	
	// getting p
	vdoub_t p = load_p(pp_file, number_of_nodes, num_iter);

	fin.open(test_file);
	if (! fin.is_open()) {
		cerr << "cost: file not found: " << test_file << endl;		
	}

	// getting uniform:
	vdoub_t unif(number_of_nodes, 1.0/number_of_nodes);

	// getting out_degree:
	vdoub_t out_degree(number_of_nodes);
	for (int u = 0; u<number_of_nodes; ++u) {
		out_degree[u] = (double) g.outdeg[u]/(g.number_of_edges);
	}
	
	// getting in_degree:
	vdoub_t in_degree(number_of_nodes);
	for (int u = 0; u<number_of_nodes; ++u) {
		in_degree[u] = (double) g.indeg[u]/(g.number_of_edges);
	}


	vdoub_t cost_val(4,0);	
	int k, u;
	double pS_p, pS_unif, pS_out, pS_in;	
	double cost_p = 0, cost_unif = 0, cost_out = 0, cost_in = 0;
	while (fin >> k) {
		// vector<double> pS(pp.size(), 0);
		pS_p = 0;
		pS_unif = 0;
		pS_out = 0;
		pS_in = 0;

		for (int i=0; i<k; ++i) {			
			fin >> u;
			pS_p += p[u];
			pS_unif += unif[u];
			pS_out += out_degree[u];
			pS_in  += in_degree[u];
			
		}

		cost_p += 1.0/(1 - theta*pow(1-pS_p,probes));
		cost_unif += 1.0/(1 - theta*pow(1-pS_unif,probes));
		cost_out += 1.0/(1 - theta*pow(1-pS_out,probes));
		cost_in += 1.0/(1 - theta*pow(1-pS_in,probes));

	}
	fin.close();

	cost_p /= len_samp;
	cost_unif /= len_samp;
	cost_out /= len_samp;
	cost_in /= len_samp;
	
	ofstream fout_comp(samp_file + "_probes-" + to_string(probes) + ".compare");
	fout_comp << "p\t" << cost_p << endl;
	fout_comp << "unif\t" << cost_unif << endl;
	fout_comp << "out_degree\t" << cost_out << endl;
	fout_comp << "in_degree\t" << cost_in << endl;

	fout_comp.close();
}

void compare_with_file(const string samp_file, const string test_file,
			const int num_iter, const int probes,
			const double theta) {

	

	int number_of_nodes, len_samp;
	
	ifstream fin(test_file + ".stat");	
	fin >> number_of_nodes >> len_samp;	
	fin.close();

	string pp_file = samp_file + "_probes-" + to_string(probes)+ ".p";
	
	// getting p
	vdoub_t p = load_p(pp_file, number_of_nodes, num_iter);

	fin.open(test_file);
	if (! fin.is_open()) {
		cerr << "cost: file not found: " << test_file << endl;		
	}

	// getting uniform:
	vdoub_t unif(number_of_nodes, 1.0/number_of_nodes);

	


	vdoub_t cost_val(2,0);	
	int k, u;
	double pS_p, pS_unif;
	double cost_p = 0, cost_unif = 0;
	while (fin >> k) {
		// vector<double> pS(pp.size(), 0);
		pS_p = 0;
		pS_unif = 0;

		for (int i=0; i<k; ++i) {			
			fin >> u;
			pS_p += p[u];
			pS_unif += unif[u];
			
		}

		cost_p += 1.0/(1 - theta*pow(1-pS_p,probes));
		cost_unif += 1.0/(1 - theta*pow(1-pS_unif,probes));


	}
	fin.close();

	cost_p /= len_samp;
	cost_unif /= len_samp;

	
	ofstream fout_comp(samp_file + "_probes-" + to_string(probes) + ".compare");
	fout_comp << "p\t" << cost_p << endl;
	fout_comp << "unif\t" << cost_unif << endl;


	fout_comp.close();
}

	



int main(int argc, char *argv[]) {
	/*
		Parsing the arguments
	*/
	string TASK, filename, SFILE, TFILE;
	uint64_t len_samp, len_test;

	int num_iter = 51, probes, num_simuls = 1, ind; 
	double epsilon = 0.1, theta = 0.75;

	for (int i=0; i<argc; ++i) {			
		if (! strcmp(argv[i],"-TASK")) {			
			TASK = argv[i+1];
		} else if (! strcmp(argv[i],"-SAMPLE")) {
			len_samp = stoi(argv[i+1]);
		} else if (! strcmp(argv[i],"-TEST")) {
			len_test = stoi(argv[i+1]);
		} else if (! strcmp(argv[i],"-DIR")) {
			filename = argv[i+1];
			filename +=  "graph";			
		} else if (! strcmp(argv[i],"-EPSILON")) {
			epsilon = stod(argv[i+1]);
		} else if (! strcmp(argv[i],"-THETA")) {
			theta = stod(argv[i+1]);
		} else if (! strcmp(argv[i],"-NUM_ITER")) {
			num_iter = stod(argv[i+1]);
		} else if (! strcmp(argv[i],"-PROBES")) {
			probes = stod(argv[i+1]);
		} else if (! strcmp(argv[i],"-NUM_SIMULS")) {
			num_simuls = stod(argv[i+1]);
		} else if (! strcmp(argv[i],"-IND")) {
			ind = stod(argv[i+1]);
		} else if (! strcmp(argv[i],"-SAMPLE_FILE")) {
			SFILE = argv[i+1];
		} else if (! strcmp(argv[i],"-TEST_FILE")) {
			TFILE = argv[i+1];
		}

	}

	/*
		Running the tasks
	*/



	if (TASK == "sampler") {
		Graph g(filename, theta);
		sampler(g, my_get_dir_name(filename), len_samp);
	} else if (TASK == "indexed_sampler") {
		Graph g(filename, theta);
		indexed_sampler(g, my_get_dir_name(filename), len_samp, ind);

	} else if (TASK == "solver") {
		string samp_file = my_get_dir_name(filename)  +"samples/S-"+to_string(len_samp);
		solver(samp_file, theta, num_iter, probes);	

	} else if (TASK == "solver-with-file") {
		string samp_file = SFILE;
		solver(samp_file, theta, num_iter, probes);	

	} else if (TASK == "indexed_solver") {
		string samp_file = my_get_dir_name(filename)  +"indexed/S-"+to_string(len_samp)+"-"+to_string(ind);
		solver(samp_file, theta, num_iter, probes);	

	} else if (TASK == "tester") {
		string samp_file = my_get_dir_name(filename)  +"samples/S-"+to_string(len_samp);
		string test_file = my_get_dir_name(filename)  +"samples/S-"+to_string(len_test);
		tester(samp_file, test_file, num_iter, probes, theta);

	} else if (TASK == "indexed_tester") {
		string samp_file = my_get_dir_name(filename)  +"indexed/S-"+to_string(len_samp)+"-"+to_string(ind);
		string test_file = my_get_dir_name(filename)  +"indexed/S-"+to_string(len_test)+"-"+to_string(ind);
		tester(samp_file, test_file, num_iter, probes, theta);
		
	} else if (TASK == "simulator") {
		simulator(filename, num_simuls, num_iter, theta, epsilon);
	
	} else if (TASK == "compare") {
		string samp_file = my_get_dir_name(filename)  +"samples/S-"+to_string(len_samp);
		string test_file = my_get_dir_name(filename)  +"samples/S-"+to_string(len_test);
		compare(filename, samp_file, test_file, num_iter, probes, theta);	
	} else if (TASK == "compare-with-file") {
		string samp_file = SFILE;
		string test_file = TFILE;
		compare_with_file(samp_file, test_file, num_iter, probes, theta);	
	
	} else if (TASK == "indexed_compare") {
		string samp_file = my_get_dir_name(filename)  +"indexed/S-"+to_string(len_samp)+"-"+to_string(ind);
		string test_file = my_get_dir_name(filename)  +"indexed/S-"+to_string(len_test)+"-"+to_string(ind);
		compare(filename, samp_file, test_file, num_iter, probes, theta);
	}
	
	return 0;
}
