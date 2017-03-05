# project2_PSO
To run: 
   "make"
   "./PSO neighborhoodTopology  swarmSize  numIterations  function  numDimensions"
for example:
    "./PSO gl 100 10000 ack 30"

You will be prompted if you want to run 720 test of the PSO. type "n" and enter if you do not want to. If you do want to, you will be prompted again to give a filename for the ouput file. Add a .csv extension. This method will test each neighborhood topology with each function with numParticles = 16,30, and 49. The output file will appear in the current directory.  