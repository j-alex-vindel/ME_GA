
class Population(object):
    """ Population representation - a group of individuals
        populations can merge together"""
    def __init__(self):
        self.population = []
        self.fronts = []

    def __len__(self):
        return len(self.population)

    def __iter__(self):
        """Allows for individual's iteration"""
        return self.population.__iter__()

    def extend(self,new_individuals):
        """Creates a new population that consists of
        old and new individuals"""
        self.population.extend(new_individuals)