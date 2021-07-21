class Population:
    # '''Encapsulates a population of chromosomes represented by a vector of doubles. We don't decode
    # nor deal with random numbers here; instead, we provide support methods to set the fitness of a 
    # specific chromosome as well as access methods to each allele.'''
    def __init__ (self, pop=None, n=None, p=None):
        if n != None and p != None:
            if p == 0:
                raise IndexError, "Population size p cannot be zero."
            if n == 0:
                raise IndexError, "Chromosome size n cannot be zero."
            self.__population = [] #'''Population as vectors of prob.'''
            self.__fitness = [] #'''Fitness (double) of a each chromosome'''
            for i in range (p):
                self.__population.append ([0.0]*n)
                self.__fitness.append ([0.0, i])
        elif pop != None:
            self.__population = pop.getPopulation()
            self.__fitness = pop.getF()
        
    def getN(self):
        #'''Size of each chromosome'''
        return len (self.__population [0])

    def getP(self):
        #'''Size of population'''
        return len (self.__population)

    def getBestFitness (self):
        # '''These methods REQUIRE fitness to be sorted, and thus a call to sortFitness() beforehand
        # (this is done by BRKGA, so rest assured: everything will work just fine with BRKGA).
        # Returns the best fitness in this population:'''
        return self.getFitness (0)

    def getFitness (self, i):
        # '''Returns the fitness of chromosome i \in {0, ..., getP() - 1}'''
        if i >= self.getP():
            raise IndexError, "Invalid individual identifier."
        return self.__fitness[i][0]

    def getFitnessPair (self, i):
        return self.__fitness[i][:]
    
    def getF (self):
        return self.__fitness[:]
	
    def getPopulation (self):
        return self.__population[:]

    def getChromosome (self, i):
        # '''Returns a chromosome'''
        if i >= self.getP():
            raise IndexError, "Invalid individual identifier."
        return self.__population [self.__fitness[i][1]]

    def setChromosome (self, i, new_chromosome):
        if i >= self.getP():
            raise IndexError, "Invalid individual identifier."
        self.__population [i] = new_chromosome[:]

    def getAllele (self, chromossome, allele):
        return self.__population[chromossome][allele]

    def setAllele (self, chromossome, allele, value):
        self.__population[chromossome][allele] = value

    def sortFitness (self):
        # '''Sorts 'fitness' by its first parameter'''
        self.__fitness.sort()

    def setFitness (self, i, f):
        # '''Sets the fitness of chromosome i'''
        self.__fitness[i][0] = f
        self.__fitness[i][1] = i
        
    def setFitnessValue (self, i, f):
        self.__fitness[i][0] = f
    
