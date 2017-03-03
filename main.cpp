#include "main.h"
#include <stdlib.h>
#include <stdio.h>

using namespace std;


int main(int argc, char** argv) {
    string neighborhoodTopology = argv[1];
    int swarmSize = stoi(argv[2]);
    int numIterations = stoi(argv[3]);
    string function = argv[4];
    int numDimensions = stoi(argv[5]);
    PSO pso_alg(neighborhoodTopology, swarmSize, numIterations, function, numDimensions);
    pso_alg.runPSO();
    return 0;
}
