from multiprocessing import Pool 

from Population import *
from IDecoder import *
from IRNG import *

class BRKGA:
	# '''
	# This class encapsulates a Biased Random-key Genetic Algorithm for minimization problems
	# with K  independent Populations stored in two vectors of Population, current and previous.
	# It supports multi-threading via the multiprocessing lib, and implements the following key methods:
	# BRKGA() constructor: initializes the populations with parameters described below.
	# - evolve() operator: evolve each Population following the BRKGA methodology. This method
                       # supports the multiprocessing lib to evolve up to K independent Populations in parallel.
                       # Please note that double Decoder.decode(...) MUST be thread-safe.
					   
	# Required parameters:
	# - n: number of genes in each chromosome
	# - p: number of elements in each population
	# - pe: pct of elite items into each population
	# - pm: pct of mutants introduced at each generation into the population
	# - rhoe: probability that an offspring inherits the allele of its elite parent
	
	# Optional parameters:
	# - K: number of independent Populations (set to 1 if not supplied)
	# - MAX_THREADS: number of threads to perform parallel decoding (set to 1 if not supplied)
                 # WARNING: Decoder.decode(...) MUST be thread-safe if MAX_THREADS > 1!
 
 	# The following objects are required:
	# RNG: random number generator that implements the methods below.
      # - RNG(seed) to initialize a new RNG with 'seed'
      # - rand() to return a double precision random deviate in range [0,1)
      # - randInt() to return a >=32-bit unsigned random deviate in range [0,2^32-1)
      # - randInt(N) to return a unsigned random deviate in range [0, N] with N < 2^32
	  
	# Decoder: problem-specific decoder that implements any of the decode methods outlined below. When
           # compiling and linking BRKGA with -fopenmp (i.e., with multithreading support via
           # multiprocessing), the method must be thread-safe.
      # - decode(chromosome)
	  # '''

    def __init__ (self, n, p, pe, pm, rhoe, refDecoder, refRNG, K=1, MAX_THREADS = 1):
		#'''Error check:'''
        if n == 0:
            raise IndexError, "Chromosome size n cannot be zero."
        if p == 0:
            raise IndexError, "Population size p cannot be zero."
        if pe == 0:
            raise IndexError, "Elite-set size pe cannot be zero."
        if pe > p:
            raise IndexError, "Elite-set size cannot be greater than population size (pe > p)."
        if pe > p:
            raise IndexError, "Mutant-set size (pm) cannot be greater than population size (p)."
        if pm + pe > p:
            raise IndexError, "Elite + mutant sets cannot be greater than population size (p)."
        if K == 0:
            raise IndexError, "Number of parallel populations cannot be zero."
        self.__n = n
        self.__p = p
        self.__pe = int (pe*p)
        self.__pm = int (pm*p)
        self.__rhoe = rhoe
        self.__refRNG = refRNG
        self.__refDecoder = refDecoder
        self.__K = K
        self.__MAX_THREADS = MAX_THREADS
        self.__current = [] #'''current populations'''
        self.__previous = [] #'''previous populations'''
		
		#'''Initialize and decode each chromosome of the current population, then copy to previous:'''
        for i in range (K):
			
			#'''Allocate:'''
            self.__current.append (Population (None, n,p))
			
			#'''Initialize:'''
            self.initialize(i);
			
			#'''Then just copy to previous:'''
            self.__previous.append (Population (self.__current[i]))
			

    def initialize (self, i):
		#'''Initialize current population 'i' with random keys'''
        for j in range (self.__p):
            for k in range (self.__n):
                self.__current[i].setAllele (j,k,self.__refRNG.rand())

		#'''Decode:'''
        if self.__MAX_THREADS > 1:

            pool = Pool (self.__MAX_THREADS)
            tasks = [self.__current[i].getChromosome(j) for j in range (self.__p)]
            results = pool.map (self.__refDecoder.decode, tasks)
            pool.close()

            for j in range (self.__p):
                self.__current[i].setFitness (j, results[j])

        else:
            for j in range (self.__p):
                self.__current[i].setFitness (j, self.__refDecoder.decode(self.__current[i].getChromosome(j)[:]))

		#'''Sort:'''
        self.__current[i].sortFitness()

    def evolution (self, curr__pop, next_pop):
		#'''We now will set every chromosome of 'current', iterating with 'i':'''
        i = 0

        #'''The 'pe' best chromosomes are maintained, so we just copy these into 'current':'''
        while i < self.__pe: #'''Iterate chromosome by chromosome'''
            for j in range (self.__n): #'''Iterate allele by allele'''
                next_pop.setAllele (i,j,curr__pop.getAllele(curr__pop.getFitnessPair(i)[1],j))

            next_pop.setFitness (i,curr__pop.getFitness(i))
            i+=1

		#'''We'll mate 'p - pe - pm' pairs; initially, i = pe, so we need to iterate until i < p - pm:'''
        while i < self.__p-self.__pm:
			
			#'''Select an elite parent:'''
            eliteParent = self.__refRNG.randInt (self.__pe-1)
			
			#'''Select a non-elite parent:'''
            noneliteParent = self.__pe + self.__refRNG.randInt (self.__p - self.__pe - 1)

			#'''Mate:'''
            for j in range (self.__n):
                if self.__refRNG.rand() < self.__rhoe:
                    sourceParent = eliteParent
                else:
                    sourceParent = noneliteParent

                next_pop.setAllele (i,j,curr__pop.getAllele(curr__pop.getFitnessPair(sourceParent)[1],j))

            i+=1

		#'''We'll introduce 'pm' mutants:'''
        while i<self.__p:
            for j in range (self.__n):
                next_pop.setAllele (i,j,self.__refRNG.rand())
            i+=1

		#'''Time to compute fitness, in parallel:'''
        if self.__MAX_THREADS > 1:

            pool = Pool (self.__MAX_THREADS)
            tasks = [next_pop.getChromosome(i) for i in range (self.__pe, self.__p)]
            results = pool.map (self.__refDecoder.decode, tasks)
            pool.close()
            
            for i in range (self.__pe, self.__p):
                next_pop.setFitness (i, results[i-self.__pe])
            
        else:

            for i in range (self.__pe, self.__p):
                next_pop.setFitness (i,self.__refDecoder.decode(next_pop.getChromosome(i)[:]))
                
        #'''Now we must sort 'current' by fitness, since things might have changed:'''
		next_pop.sortFitness()


    def reset (self):
		#'''Resets all populations with brand new keys'''
        for i in range (self.__K):
            self.initialize(i)

    def evolve (self, generations = 1):
		# '''Evolve the current populations following the guidelines of BRKGAs
        # @param generations number of generations (must be even and nonzero)
        # @param J interval to exchange elite chromosomes (must be even; 0 ==> no synchronization)
        # @param M number of elite chromosomes to select from each population in order to exchange'''
        if generations == 0:
            raise IndexError, "Cannot evolve for 0 generations."


        for i in range (generations):
            for j in range (self.__K):
                self.evolution (self.__current[j], self.__previous[j]) #'''First evolve the population (curr, next)'''
				
				#'''Update (prev = curr; curr = prev == next)'''
                self.__current[j], self.__previous[j] = self.__previous[j], self.__current[j]
                

    def exchangeElite (self, M):
		# '''Exchange elite-solutions between the populations
        # @param M number of elite chromosomes to select from each population'''
        if M == 0 or M >= self.__p:
            raise IndexError, "M cannot be zero or >= p."
		
		#'''Population i will receive some elite members from each Population j below:'''
        for i in range (self.__K):
            dest = self.__p - 1 #'''Last chromosome of i (will be updated below)'''
            for j in range (self.__K):
                if j == i:
                    continue
					
				#''' Copy the M best of Population j into Population i:'''
                for m in range (M):
					#'''Copy the m-th best of Population j into the 'dest'-th position of Population i:'''
                    bestOfJ = self.__current[j].getChromosome(m)
                    self.__current[i].setChromosome (dest, bestOfJ[:])
                    self.__current[i].setFitnessValue (dest, self.__current[j].getFitness(m))
                    dest-=1

        for j in range (self.__K):
            self.__current[j].sortFitness()

    def getPopulation (self, K=0):
		#'''Returns the current population'''
        if(K >= self.__K):
            raise IndexError, "Invalid population identifier."

        return self.__current[K]

    def getBestChromosome (self):
		#'''Returns the chromosome with best fitness so far among all populations'''
        bestK = 0
        for i in range (1, self.__K):
            if self.__current[i].getBestFitness() < self.__current[bestK].getBestFitness():
                bestK = i

        return self.__current[bestK].getChromosome(0)

    def getBestFitness (self):
		#'''Returns the best fitness found so far among all populations'''
        best = self.__current[0].getFitness (0)
        for i in range (1,self.__K):
            if self.__current[i].getFitness (0) < best:
                best = self.__current[i].getFitness(0)

        return best

	#'''Return copies to the internal parameters:'''
    def getN (self):
        return self.__n

    def getP (self):
        return self.__p
    
    def getPe (self):
        return self.__pe
    
    def getPm (self):
        return self.__pm

    def getPo (self):
        return self.__p - self.__pe - self.__pm

    def getRhoe (self):
        return self.__rhoe

    def getK (self):
        return self.__K

    def getMAX_THREADS (self):
        return self.__MAX_THREADS
