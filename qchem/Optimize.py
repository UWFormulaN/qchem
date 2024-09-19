from scipy import optimize

# Implement energy minimization algorithms here (eventually convert to C).

class GeoOpt:
    def __init__(self, molecule, method):
        self.molecule = molecule
        self.method = method
