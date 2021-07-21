from IDecoder import *

class SampleDecoder (IDecoder):
    
    def decode (self, chromosome):
        myFitness = 0
        
        for i in range (len(chromosome)-1):
            myFitness+=(i+1)*chromosome[i]
			
        return myFitness
