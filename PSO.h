#ifndef __PSO_h
#define __PSO_h

#include <vector>
#include <string>
#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

struct Particle {
    vector <double> velocity;
    vector <double> position;
    vector <double> p_best_position, n_best_position;
    double p_best_value, n_best_value;
    vector <int> neighborhood_indices;
};

enum Topology {GLOBAL, RING, VON_NEUMANN, RANDOM};
enum Function {ACKLEY, ROSENBROCK, RASTRIGIN};

class PSO {
public:
    PSO(string neighborhoodTopology, int swarmSize, int numIterations, string function, int numDimensions);
    vector<double> runPSO();

private:
    Topology neighborhood_topo;
    Function function_to_optimize;
    int swarm_size, num_iterations, num_dimensions, k_neighbors;
    double g_best_value, phi1, phi2, constriction_factor;
    double min_position, max_position, min_velocity, max_velocity;
    vector <double> g_best_position;
    vector <Particle> swarm;

    void initialize_swarm();
    void evaluate_neighborhoods();
    void create_neighborhoods();
    void create_ring();
    void create_von_neumann();
    void create_random();
    void update_velocities();
    void update_positions();
    void evaluate_swarm();
    void iterate();
    double function_value(vector<double> position);
    double rosenbrock_function(vector<double> position);
    double ackley_function(vector<double> position);
    double rastrigin_function(vector<double> position);

    // utility functions
    vector<double> vector_subtraction(vector<double> vector1, vector<double> vector2);
    double rand_in_range(double min, double max);

};


#endif
