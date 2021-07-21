from IRNG import *
from random import *

class PyRand (IRNG):

    def __init__ (self, rnd_seed = None):
        seed (rnd_seed)

    def rand (self):
        return random ()

    def randInt (self, N):
        return randint (0, N)
