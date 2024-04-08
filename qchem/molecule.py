import pandas as pd
from scipy import optimize

# We will create molecule objects which will store information about the molecule
# Includes coordinates, atom types, how optimization was performed, how energy calculations were performed, etc.
# We can take segmented properties from the output file, and store them as attributes of the molecule object


class Molecule:
    def __init__(
        self,
        name,
        multiplicity,
        coordinates,
        atom_types,
        optimization_method,
        energy_method,
        energy,
    ):
        self.name = name
        self.multiplicity = multiplicity
        self.coordinates = coordinates
        self.atom_types = atom_types
        self.optimization_method = optimization_method
        self.energy_method = energy_method
        self.energy = energy
