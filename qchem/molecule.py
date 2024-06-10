from ast import List
from os import read
import string
import pandas as pd
from scipy import optimize
from .Data.constants import CovalentRadiiConstants

# We will create molecule objects which will store information about the molecule
# Includes coordinates, atom types, how optimization was performed, how energy calculations were performed, etc.
# We can take segmented properties from the output file, and store them as attributes of the molecule object


class Molecule:

    # The Following are Variables that don't need to be initialized immediately or specified by the user, but belong to the Molecule class
    # For Some reason the Description needs to be defined after?

    name: str = ""
    """Name of the Molecule"""

    #  Atom Symbol, X, Y, Z ?
    XYZCoordinates: pd.core.frame.DataFrame
    """List of the XYZ Coordinates of the Atoms in the Molecule"""

    energy = None
    """The Energy of the Molecule"""

    multiplicity = None
    """Multiplicity of the Molecule"""

    atom_types = None
    """Atom types of the Molecule"""

    optimization_method = None
    """Optimization Method of the Molecule"""

    energy_method = None
    """Energy Method used for the Molecule"""

    def ReadXYZ(self, path: str) -> pd.core.frame.DataFrame:
        return pd.read_csv(
            path, sep=r"\s+", skiprows=1, names=["Atom", "X", "Y", "Z"], engine="python"
        )

    def __init__(self, name: str, XYZFilePath: str):
        self.name = name

        # Load the XYZ File from XYZ File
        self.XYZCoordinates = self.ReadXYZ(path=XYZFilePath)

