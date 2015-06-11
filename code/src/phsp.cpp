#include <cmath>
#include "graph.h"
// #include <direct.h>

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
	ofstream fout_p(samp_file+".p");	
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
	
	ofstream fout_conv(samp_file + ".conv");
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

	string pp_file = samp_file +".p";
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

	
	ofstream fout_test(samp_file + ".test");
	for (double c : cost_val) {			
			fout_test << c << "\t";
	}
	fout_test << endl;
}






// ************************************************
// ************************************************



int main(int argc, char *argv[]) {
	/*
		Parsing the arguments
	*/
	string TASK, filename;
	uint64_t len_samp, len_test;

	int num_iter = 51, probes; // arm counts the number of probes
	double epsilon = 0.1, theta = 0.75;

	for (int i=0; i<argc; ++i) {			
		if (! strcmp(argv[i],"-TASK")) {			
			TASK = argv[i+1];
		} else if (! strcmp(argv[i],"-SAMPLE")) {
			len_samp = stoi(argv[i+1]);
		} else if (! strcmp(argv[i],"-TEST")) {
			len_test = stoi(argv[i+1]);
		} else if (! strcmp(argv[i],"-FILE")) {
			filename = argv[i+1];
		} else if (! strcmp(argv[i],"-EPSILON")) {
			epsilon = stod(argv[i+1]);
		} else if (! strcmp(argv[i],"-THETA")) {
			theta = stod(argv[i+1]);
		} else if (! strcmp(argv[i],"-NUM_ITER")) {
			num_iter = stod(argv[i+1]);
		} else if (! strcmp(argv[i],"-PROBES")) {
			probes = stod(argv[i+1]);
		}
	}

	/*
		Running the tasks
	*/



	if (TASK == "sampler") {
		Graph g(filename, theta);
		sampler(g, my_get_dir_name(filename), len_samp);

	} else if (TASK == "solver") {
		string samp_file = my_get_dir_name(filename)  +"samples/S-"+to_string(len_samp);
		solver(samp_file, theta, num_iter, probes);		

	} else if (TASK == "tester") {
		string samp_file = my_get_dir_name(filename)  +"samples/S-"+to_string(len_samp);
		string test_file = my_get_dir_name(filename)  +"samples/S-"+to_string(len_test);
		tester(samp_file, test_file, num_iter, probes, theta);
	}
	
	return 0;
}