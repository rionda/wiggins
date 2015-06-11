#include "graph.h"
#include <random>
#include <cmath>
// #include <set>

random_device rd;
mt19937 eng(rd());




// double pr_group(const int degree) {	
// 	// return 2*atan(degree)/3.14159265;
// 	return 1;
// }

double pr_group(const int degree) {	
	if (degree >= 1000) {		
		return 0.1;
	} else if (degree >= 500) {
		return 0.05;
	} else if (degree >= 100) {
		return 0.01;
	} else {
		return 0.0;
	}
}


void update(vector<double>& p, const vector<double>& W) {
	for (int i=0; i < p.size(); ++i) {
		p[i] *= W[i];
	}
	double sum = accumulate(p.begin(), p.end(), 0.0);
	for (int i=0; i < p.size(); ++i) {
		p[i] /= sum;
	}
}

Graph::Graph(const string input_file, const double theta_input) {
	theta = theta_input;
	graph_name = input_file;
	// influence_type: suppose you have an edge (a,b). This is =0 if b influences a, and = 1
	// if a influences b.
	const string direction =  input_file.substr(input_file.find_last_of(".") + 1);
	cout << "The graph is " << direction << endl;
	
	ifstream fin;
	fin.open(input_file);
	if (! fin.is_open()) {
		cout << "file not found: " << input_file << endl;
		return;
	} 
	unset_t V;
	int u, v;
	double p;
	while (fin >> u >> v >> p) {
		V.insert(u);
		V.insert(v);
	}
	fin.close();
	
	number_of_nodes = V.size();
	nodes = vvint_t(number_of_nodes);
	distr = uniform_int_distribution<>(0,number_of_nodes-1);
	bern = vector<vector<bernoulli_distribution>>(number_of_nodes);

	// loading the edges:
	number_of_edges = 0;
	fin.open(input_file);
	while (fin >> u >> v >> p) {
		++number_of_edges;		
		nodes[u].push_back(v);
		bern[u].push_back(bernoulli_distribution(p));
	}	
	fin.close();

	// keep adjs sorted
	for(int i=0; i < number_of_nodes; ++i) {
		sort(nodes[i].begin(), nodes[i].end());		
	}	
}


set<int> Graph::random_spot(int source) {	
	unset_t Q;
	set<int> spot;

	Q.insert(source);
	int u, v;
	while(!Q.empty()) {
		u = *Q.begin();		
		Q.erase(u);
		spot.insert(u);
		for (int v_idx=0; v_idx < nodes[u].size(); ++v_idx) {
			v = nodes[u][v_idx];
			if(bern[u][v_idx](eng)) {
				if (!spot.count(v)) {
					Q.insert(v);
				}  
			}
		}
	}
	return spot;
}

void Graph::gen_sample(const int time_interval, const string output_file) {
	ofstream fout;
	fout.open(output_file);
	for (int u=0; u < number_of_nodes; ++u) {
		if (u%1000 == 0) {
			cout << u << ", "; cout.flush();
		}
		bernoulli_distribution b(pr_group(nodes[u].size()));
		for (int i=0; i < time_interval; ++i) {
			if (b(eng)) {
				set<int> spot = random_spot(u);
				fout << spot.size() << "\t";
				for (int n : spot) {
					fout << n << "\t";
				}
				fout << endl;
			}
		}
	}
	fout.close();
	cout << endl;
}

void Graph::gen_sample_with_time(const int time_interval, const string output_file, const int time_idx) {	
	ofstream fout;
	fout.open(output_file);
	for (int u=0; u < number_of_nodes; ++u) {
		
		// cout << u << ", "; cout.flush();
		
		bernoulli_distribution b(pr_group(nodes[u].size()));
		for (int i=0; i < time_interval; ++i) {
			if (b(eng)) {
				set<int> spot = random_spot(u);
				fout << spot.size() << "\t" << i+time_idx << "\t";
				for (int n : spot) {
					fout << n << "\t";
				}
				fout << endl;
			}
		}
	}
	fout.close();
}


vdoub_t Graph::simul_process(double eps, int num_iter, string idxstr) {
	int u;
	vdoub_t p0, p1, p2;
	// int R = ceil(log(g.number_of_nodes) * pow(theta,2) /(2*pow(eps,2)*pow(1-theta,2)));
	int R = ceil(3*(log(number_of_nodes)+log(2))/ ((1-theta)*pow(eps,2)));
	// int R = 10;
	int a = 10*R;
	int b = 2*R;
	// int a = 100*R;
	// int b = 20*R;
	cout << "R: " << R << ", a: " << a << ", b: " << b << endl;
	
	p0 = simul_with_file(R, graph_name+".tmp0"+idxstr, num_iter, 0);
	p1 = simul_with_file(R, graph_name+".tmp1"+idxstr, num_iter, a+b-R);
	p2 = simul_with_file(R, graph_name+".tmp2"+idxstr, num_iter, 2*(a+b)-R);

	discrete_distribution<int> disc0(p0.begin(), p0.end());
	discrete_distribution<int> disc1(p1.begin(), p1.end());
	discrete_distribution<int> disc2(p2.begin(), p2.end());

	// vvint_t schedule(3*a + 2*b);
	vvint_t schedule(number_of_nodes);
	//
	cout << "1" << endl; cout.flush();
	cout << "schedule size: " << schedule.size() << endl; cout.flush();
	cout << "p0 size: " << p0.size() << endl; cout.flush();
	for (int t=0; t< a; ++t) {
		u = disc0(eng);
		schedule[u].push_back(t);
	}
	cout << "1 (shuffle)" << endl; cout.flush();
	shuffle(p0.begin(), p0.end(), eng);
	disc0 = discrete_distribution<int>(p0.begin(), p0.end());
	for (int t=a; t< (a+b); ++t) {
		u = disc0(eng);
		schedule[u].push_back(t);
	}
	//
	cout << "2" << endl; cout.flush();
	for (int t=(a+b); t< 2*a+b; ++t) {
		u = disc1(eng);
		schedule[u].push_back(t);
	}
	cout << "2(shuffle)" << endl; cout.flush();
	shuffle(p1.begin(), p1.end(), eng);
	disc1 = discrete_distribution<int>(p1.begin(), p1.end());
	for (int t=2*a+b; t< 2*(a+b); ++t) {
		u = disc1(eng);
		schedule[u].push_back(t);
	}
	//
	cout << "3" << endl; cout.flush();
	for (int t=2*(a+b); t< 3*a+2*b; ++t) {
		u = disc2(eng);
		schedule[u].push_back(t);
	}


	// go for cost:
	vdoub_t cost(3*a+2*b, 0);
	cout << "online(1)" << endl; cout.flush();
	go_online(schedule, cost, 0, a+b-R);
	cout << "offline (1)" << endl; cout.flush();
	go_offline(schedule, cost, graph_name+".tmp1"+idxstr);
	cout << "online(2)" << endl;cout.flush();
	go_online(schedule, cost, a+b, a+b-R);
	cout << "offline (2)" << endl;cout.flush();
	go_offline(schedule, cost, graph_name+".tmp2"+idxstr);
	cout << "online(3)" << endl;cout.flush();
	go_online(schedule, cost, 2*(a+b), a);

	return cost;
}


vector<double> Graph::simul_with_file(const int time_interval, const string sample_file, const int num_iter, const int time_idx) {
	// cout << ">>>>>>> 1 \n"; cout.flush();
	// clock_t t=clock();
	// gen_sample(time_interval, sample_file+"XXX");
	// cout << "took: " << clock() - t << " milsec\n";flush();
	// cout << ">>>>>>> 2 \n"; cout.flush();
	// t = clock();
	gen_sample_with_time(time_interval, sample_file, time_idx);
	// cout << "took: " << clock() - t << " milsec\n";cpitflush();
	
	return solver(sample_file, time_interval, num_iter);
}

vector<double> Graph::solver(const string sample_file, const int time_length, const int num_iter) {
	vector<double> p(number_of_nodes, 1.0/number_of_nodes);	
	int k, u;
	double pS, tmp;
	ifstream fin;
	
	double time_tmp;

	for (int r=0; r < num_iter; ++r) {		
		fin.open(sample_file);
		vector<double> W(number_of_nodes, 0);		
		while (fin >> k) {
			fin >> time_tmp;
			pS = 0;
			deque<int> S; // TODO: check if it gets new everytime
			for (int i=0; i<k; ++i) {
				fin >> u;
				pS += p[u];
				S.push_back(u);
			}
			tmp = 1.0/(1 - theta*(1-pS));
			
			// update W:
			tmp = pow(tmp,2);
			for (int v : S) {
				W[v] += tmp;
			}
		}
		fin.close();
		update(p, W);		
	}	
	return p;
}

void Graph::go_online(vvint_t & schedule, vector<double>& cost, const int time_idx, const int time_interval) {
	set<int> spot;
	for (int u=0; u < number_of_nodes; ++u) {
		bernoulli_distribution b(pr_group(nodes[u].size()));
		for (int i=0; i < time_interval; ++i) {
			if (b(eng)) {
				spot = random_spot(u);
				impose(schedule, spot, i+time_idx, cost);				
			}
		}
	}
}

void Graph::go_offline(vvint_t & schedule, vector<double>& cost, const string sample_file) {	
	ifstream fin;
	fin.open(sample_file);
	int k, u, current_time;

	while (fin >> k) {
		set<int> spot;
		fin >> current_time;
		for (int i=0; i<k; ++i) {
			fin >> u;
			spot.insert(u);
		}
		impose(schedule, spot, current_time, cost);
	}
}


void Graph::impose(vvint_t& schedule, set<int>& spot, const int current_time, vector<double>& cost) {
	int catch_time = cost.size();
	for (int u : spot) {
		auto res = lower_bound(schedule[u].begin(),schedule[u].end(), current_time);
		if (res != schedule[u].end()) {
			if (*res < catch_time) {
				catch_time = *res;
			}
		}
	}	
	for (int i=current_time; i<= catch_time; ++i) {
		cost[i] += pow(theta, i-current_time);
	}
}


// void Graph::impose(vvint_t& schedule, set<int>& spot, const int current_time, vector<double>& cost) {
// 	int catch_time = cost.size();
// 	for (int u : spot) {
// 		for (int t : schedule[u]) {
// 			if (t >= current_time && t < catch_time) {
// 				catch_time = t;
// 				break;
// 			}
// 		}
// 	}	
// 	for (int i=current_time; i<= catch_time; ++i) {
// 		cost[i] += pow(theta, i-current_time);
// 	}
// }

// void Graph::simulate(vvint_t & schedule, vector<double>& cost) {
// 	int time_interval = cost.size();	
// 	set<int> spot;

// 	for (int u=0; u < number_of_nodes; ++u) {
// 		bernoulli_distribution b(pr_group(nodes[u].size()));
// 		for (int i=0; i < time_interval; ++i) {
// 			if (b(eng)) {
// 				spot = random_spot(u);
// 				impose(schedule, spot, i, cost);				
// 			}
// 		}
// 	}
// }

void Graph::inf_gen_sample(const int num_samples, const string output_file) {
	ofstream fout;
	fout.open(output_file);

	int u;
	set<int> spot;
	// cout << "sample counter: ";
	for (int j=0; j < num_samples; ++j) {
		// cout << j << ", ";
		u = distr(eng);
		spot = random_spot(u);
		fout << spot.size() << "\t";
		for (int n : spot) {
			fout << n << "\t";
		}
		fout << endl;
	}
	cout << endl;

	fout.close();
	cout << endl;
}
