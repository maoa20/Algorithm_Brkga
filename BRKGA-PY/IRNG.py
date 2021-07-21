from abc import ABCMeta, abstractmethod

class IRNG:
    __metaclass__ = ABCMeta
    
    @abstractmethod
    def rand(self):
        pass

    @abstractmethod
    def randInt(self, N):
        pass

