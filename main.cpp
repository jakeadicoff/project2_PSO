#include "main.h"

using namespace std;

int main(int argc, char** argv) {
  string test;
  string neighborhoodTopology = argv[1];
  int swarmSize = stoi(argv[2]);
  int numIterations = stoi(argv[3]);
  string function = argv[4];
  int numDimensions = stoi(argv[5]);
  string output = "";
  cout << "run all 720 tests?! (y or n)" << endl;
  cin >> test;
  if (test == "n") {
    PSO pso_alg(neighborhoodTopology, swarmSize, numIterations, function, numDimensions);
    pso_alg.runPSO();
    return 0;
  } 
 
  string outputFileName;
  cout << "name the output file with .csv" << endl;
  cin >> outputFileName;
  ofstream outputFile;
  output += "Neighborhood Topology,Particles,Function,G_best,1000,2000,3000,4000,5000,6000,7000,8000,9000\n";
  vector <double> values;
  for(int i = 0; i < 4; i++) {
    for(int j = 0; j < 3; j++) {
      for(int k = 0; k < 3; k++) {
	if (i == 0) { 
	  neighborhoodTopology = "gl";
	} else if (i == 1) {
	  neighborhoodTopology = "ri";
	} else if (i == 2) {
	  neighborhoodTopology = "vn";
	} else {
	  neighborhoodTopology = "ra";
	}

	if (j == 0) {
	  swarmSize = 16;
	} else if (j == 1) {
	  swarmSize = 30;
	} else {
	  swarmSize = 49;
	}

	if (k == 0) {
	  function = "rok";
	} else if (k == 1) {
	  function = "ack";
	} else {
	  function = "ras";
	}
	for(int l = 0; l < 20; l++) {
	  PSO pso_alg(neighborhoodTopology, swarmSize, numIterations, function, numDimensions);
	  values.clear();
	  values = pso_alg.runPSO();
	  output += neighborhoodTopology + "," + to_string(swarmSize) + "," + function + "," + to_string(values[10]);
	  cout << "size = " << values.size() << endl;
	  for(int m = 0; m < 10; m++) {
	    output += "," + to_string(values[i]);
	  }
	  output += "\n";
	}
      }
    }
  }
  outputFile.open(outputFileName);
  outputFile << output;
  outputFile.close();
  cout << output << endl;
return 0;
}
