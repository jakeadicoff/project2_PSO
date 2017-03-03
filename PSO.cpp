#include "PSO.h"


double BIG_DOUBLE = 9999999999;
double PI = 3.1415926535897;

PSO::PSO(string neighborhoodTopology, int swarmSize, int numIterations,
	 string function, int numDimensions) {
  this->swarm_size = swarmSize;
  this->num_iterations = numIterations;
  this->num_dimensions = numDimensions;
  this->g_best_value = BIG_DOUBLE;

  this->phi1 = 2.05;
  this->phi2 = 2.05;
  this->constriction_factor = 0.7298;
  //HEY
  this->k_neighbors = 5;

  vector <double> gBestPos;
  for(int i = 0; i < num_dimensions; i++){
    gBestPos.push_back(BIG_DOUBLE);
  }
  this->g_best_position = gBestPos;

  if(neighborhoodTopology == "gl") {
    this->neighborhood_topo = GLOBAL;
  }
  if(neighborhoodTopology == "vn") {
    this->neighborhood_topo = VON_NEUMANN;
  }
  if(neighborhoodTopology == "ri") {
    this->neighborhood_topo = RING;
  }
  if(neighborhoodTopology == "ra") {
    this->neighborhood_topo = RANDOM;
  }

  if(function == "rok") {
    this->min_position = 15;
    this->max_position = 30;
    this->min_velocity = -2;
    this->max_velocity = 2;
    this->function_to_optimize = ROSENBROCK;
  } else if(function == "ack") {
    this->min_position = 16;
    this->max_position = 32;
    this->min_velocity = -2;
    this->max_velocity = 4;
    this->function_to_optimize = ACKLEY;
  } else {
    this->min_position = 2.56;
    this->max_position = 5.12;
    this->min_velocity = -2;
    this->max_velocity = 4;
    this->function_to_optimize = RASTRIGIN;
  }
  vector <Particle> s;
  this->swarm = s;
  swarm.clear();
  srand(clock());
  initialize_swarm();

}

void PSO::runPSO() {
    double start_time = clock();
    for(int i = 0; i < num_iterations; ++i) {
        iterate();
    }
    double end_time = clock();
    cout << "Runtime: " << end_time - start_time << endl;
    cout << "Best function value: " << g_best_value << endl;
    cout << "Best position: ";
    for(int j = 0; j < num_dimensions; ++j) {
        cout << g_best_position[j] << " ";
    }
    cout << endl;

}
void PSO::initialize_swarm() {

  for(int i = 0; i < swarm_size; ++i) {
    Particle p;

    // initialize randomly in range
    for(int n = 0; n < num_dimensions; ++ n) {
      p.position.push_back(rand_in_range(min_position, max_position));
      p.velocity.push_back(rand_in_range(min_velocity, max_velocity));
    }// every dimension
    p.p_best_value = function_value(p.position);
    p.p_best_position = p.position;
    p.n_best_value = p.p_best_value;
    p.n_best_position = p.position;
    // update gbest if necessary
    if(p.p_best_value < g_best_value) {
      g_best_position = p.position;
    }
    swarm.push_back(p);
  }// every particle

  // build topology
  create_neighborhoods();
  evaluate_neighborhoods();
}


// get a random double between min and max
double PSO::rand_in_range(double min, double max) {
  return double(rand()/double(RAND_MAX)) * (max - min) + min;
}

// evaluate the neighborhoods and assign neighborhood best
void PSO::evaluate_neighborhoods() {
  switch(neighborhood_topo) {

  case GLOBAL:
    for(int i = 0; i < swarm_size; ++i) {
      swarm[i].n_best_value = g_best_value;
      swarm[i].n_best_position = g_best_position;
    }

  default: //needs fixing
    for(int i = 0; i < swarm_size; ++i) {
      for(unsigned int j = 0; j < swarm[i].neighborhood_indices.size(); ++j) {
	if(swarm[j].p_best_value < swarm[i].p_best_value) {
	  swarm[i].n_best_value = swarm[j].p_best_value;
	  swarm[i].n_best_position = swarm[j].p_best_position;
	}//end if
      }//end for
    }//end for
  }//end case
}//end func


void PSO::create_random () {
  for(int i = 0; i < swarm_size; ++i) {
    swarm[i].neighborhood_indices.push_back(i);
    for(int j = 0; j < k_neighbors - 1; ++j) {
      int random_index = i;

      // find should search vector for the value j
      while(find(swarm[i].neighborhood_indices.begin(), swarm[i].neighborhood_indices.end(), random_index) != swarm[i].neighborhood_indices.end()) {
	random_index = rand() % swarm_size;
      }
      swarm[i].neighborhood_indices.push_back(j);
    } // for every neighbor
  } // for every particle
}



void PSO::create_neighborhoods() {
  switch(neighborhood_topo) {
  case GLOBAL: break; //special logic, handled in each relevant
  case RING:
    create_ring();
    break;
  case VON_NEUMANN:
    create_von_neumann();
    break;
  case RANDOM:
    create_random();
    break;
  }
}

// give particles neighbors ahead and behind, wrap around
void PSO::create_ring() {
  for(int i = 0; i < swarm_size; ++i) {
    swarm[i].neighborhood_indices.push_back(i);
    swarm[i].neighborhood_indices.push_back((i+1)%swarm_size);
    swarm[i].neighborhood_indices.push_back((i-1) % swarm_size);
  }
}


void PSO::create_von_neumann() {
  int width = floor(sqrt(swarm_size));

  // imagine rectangle with width "width", wrap around
  for(int i = 0; i < swarm_size; ++i) {
    swarm[i].neighborhood_indices.push_back(i);
    swarm[i].neighborhood_indices.push_back((i+1)%swarm_size);
    swarm[i].neighborhood_indices.push_back((i-1) % swarm_size);
    swarm[i].neighborhood_indices.push_back((i - width) % swarm_size);
    swarm[i].neighborhood_indices.push_back((i + width) % swarm_size);
  }
}

void PSO::update_velocities() {
  for(int i = 0; i < swarm_size; ++i) {
    // HOT-SPOT FOR THE BUGS
    vector<double> p_best_difference = vector_subtraction(swarm[i].p_best_position, swarm[i].position);
    vector<double> n_best_difference = vector_subtraction(swarm[i].n_best_position, swarm[i].position);
    for(int j = 0; j < num_dimensions; ++j) {
      swarm[i].velocity[j] = constriction_factor * (swarm[i].velocity[j] + rand_in_range(0,phi1) * p_best_difference[j] + rand_in_range(0,phi2) * n_best_difference[j]);
    }
  }
}



void PSO::update_positions() {
  for(int i = 0; i < swarm_size; ++i) {
    for(int j = 0; j < num_dimensions; ++j) {
      swarm[i].position[j] = swarm[i].position[j] + swarm[i].velocity[j];
    }
  }
}

void PSO::evaluate_swarm() {
  for(int i = 0; i < swarm_size; ++i) {
    double p_value = function_value(swarm[i].position);
    if(p_value < swarm[i].p_best_value) {
      swarm[i].p_best_value = p_value;
      swarm[i].p_best_position = swarm[i].position;
      if(p_value < g_best_value) {
	g_best_value = p_value;
	g_best_position = swarm[i].position;
      }
    }
  }
}

double PSO::function_value(vector<double> position) {
  switch(function_to_optimize) {
  case ACKLEY:
    return ackley_function(position);
  case ROSENBROCK:
    return rosenbrock_function(position);
  case RASTRIGIN:
    return rastrigin_function(position);
  }
  // whats the right way to do this?
  cout << "ERROR ERROR ERROR\n" << endl;
  return BIG_DOUBLE;
}

double PSO::rosenbrock_function(vector<double> position) {
  double value = 0;
  for(int i = 0; i < num_dimensions-1; ++i) {
    value += 100 * pow(-pow(position[i], 2) + position[i+1], 2) + pow(position[i - 1] - 1, 2);
  }
  return value;
}

double PSO::ackley_function(vector<double> position) {

  double sum_part_one = 0;
  double sum_part_two = 0;

  for(int i = 0; i < num_dimensions; ++i) {
    sum_part_one += pow(position[i], 2);
    sum_part_two += (position[i], 2) * cos(2 * PI * position[i]);  // define PI in file
  }

  return -20 * exp(-0.2 * sqrt(sum_part_one / num_dimensions)) - exp(sum_part_two / num_dimensions) + 20 + exp(1);
}

double PSO::rastrigin_function(vector<double> position) {
  double value = 0;
  for(int i = 0; i < num_dimensions; ++i) {
    value += pow(position[i], 2) - 10 * cos(2 * PI * position[i]) + 10;
  }
  return value;
}




// element-wise vector subtraction
vector<double> PSO::vector_subtraction(vector<double> vector1, vector<double> vector2) {
  vector<double> vector_difference;
  for(unsigned int i = 0; i < vector1.size(); ++i) {
    vector_difference.push_back(vector1[i] - vector2[i]);
  }
  return vector_difference;
}


void PSO::iterate() {
  update_velocities();
  update_positions();
  evaluate_swarm();
  evaluate_neighborhoods();
}
