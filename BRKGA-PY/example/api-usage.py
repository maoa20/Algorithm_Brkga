# Simple example of how to instantiate and use the BRKGA API.

#Join the directory above (BRKGA's directory) to import the BRKGA module:
import sys
import os.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

from BRKGA import * #Import BRKGA module

from SampleDecoder import * #Import the Decoder implementation
from PyRand import * #Import the Random Number Generator implementation

print "Welcome to the BRKGA API sample driver.\nFinding a (heuristic) minimizer for f(x) = sum_i (x_i * i) where x \\in [0,1)^n.\n"

n = 10 #size of chromosomes
p = 100 #size of population
pe = 0.10 #fraction of population to be the elite-set
pm = 0.10 #fraction of population to be replaced by mutants
rhoe = 0.7 #probability that offspring inherit an allele from elite parent
K = 3 #number of independent populations
MAXT = 1 #number of threads for parallel decoding

decoder = SampleDecoder () #initialize the decoder
rngSeed = 12345 #seed to the random number generator
rng = PyRand (rngSeed) #initialize the random number generator

#initialize the BRKGA-based heuristic
b = BRKGA (n, p, pe, pm, rhoe, decoder, rng, K, MAXT)

generation = 0 #current generation
X_INTVL = 10 #exchange best individuals at every 100 generations
X_NUMBER = 2 #exchange top 2 best
MAX_GENS = 100 #run for 100 gens

print "Running for %d generations...\n"%(MAX_GENS)

while generation < MAX_GENS:
    b.evolve() #evolve the population for one generation
    if (generation + 1)%X_INTVL == 0:
        b.exchangeElite(X_NUMBER) #exchange top individuals
    generation += 1
	
#print the fitness of the top 10 individuals of each population:
print "Fitness of the top 10 individuals of each population: \n"

bound = min(p,20) #makes sure we have 10 individuals

for i in range (K):
    print "Population #%d\n"%(i)
    for j in range (bound):
        print "\t%d) %f\n"%(j, b.getPopulation(i).getFitness(j))

print "Best solution found has objective value = %f\n" %(b.getBestFitness())

