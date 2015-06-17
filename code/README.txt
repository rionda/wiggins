To compile:
	g++ -std=c++11 -Ofast src/graph.cpp src/phsp.cpp -o phsp


To Run:
	-TASK: a flag determines the task to be done
	-SAMPLE: The length of the time interval (for sampling)
	-TEST: the length of the time interval for test-sample
	-DIR: The directory of the database: the graph file is always called "graph" in that folder
	-EPSILON: for the epsilon parameter; default value is 0.1
	-THETA: the theta parameter; default value is 0.75
	-NUM_ITER: number of iteration in computing the schedules (in our iterative method); default 51;
	-PROBES: the number of probes at each time steps; default =1


--------------------------------------------
Examples:
***for generating a sample during a time interval of length 5000 we use:
	./phsp -TASK sampler -DIR data/twitter/ -SAMPLE 5000


***Now having tow samples one for learning (time interval of length 2000) and the other for testing (time interval of length 10000) we run:
	./phsp -TASK tester -TEST 10000 -SAMPLE 2000 -DIR data/twitter/ -PROBES 3


*** For comparing the cost of our algorithm trained by a sample of length (of time interval) 2000, and 3 other algorithms (uniform schedule, and propoertional to in/out-degrees), on a test sample of length 10000, we run:
	./phsp -TASK compare -TEST 10000 -SAMPLE 2000 -DIR data/twitter/ -PROBES 5

*** For running the simulation of applying the algorithm (the last figure of the graph), run the following commands which outputs the "load of the system" at many time steps. (using one probe at a time)

	./phsp -TASK simulator -DIR data/twitter/ -NUM_SIMULS 1

--------------------------------------------
Note that in order to use tester, sampler and compare, you need to generate the sample files in advance. 
