from ast import List
from math import sqrt
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
            path, sep=r"\s+", skiprows=2, names=["Atom", "X", "Y", "Z"], engine="python"
        )

    def GetGeometry(self):
        for atom in self.XYZCoordinates.itertuples():
            print(atom)

    def GetRadius(self, coord1, coord2):
        return sqrt(
            pow(coord1[0] - coord2[0], 2)
            + pow(coord1[1] - coord2[1], 2)
            + pow(coord1[2] - coord2[2], 2)
        )

    def GetBondGraph(self):

        at_types = self.XYZCoordinates["Atom"].values
        coords = [
            (
                self.XYZCoordinates["X"][i],
                self.XYZCoordinates["Y"][i],
                self.XYZCoordinates["Z"][i],
            )
            for i in range(0, len(self.XYZCoordinates["X"].values))
        ]
        bond_graphs = [[] for i in range(len(at_types))]

        for i in range(len(at_types)):
            radii1 = CovalentRadiiConstants[at_types[i]]
            for j in range(i + 1, len(at_types)):
                radii2 = CovalentRadiiConstants[at_types[j]]
                thresh = 1.1 * (radii1 + radii2)
                dist = self.GetRadius(coords[i], coords[j])
                if dist < thresh:
                    bond_graphs[i].append(j)
                    bond_graphs[j].append(i)

        print("%s\n" % (self.name), end="")
        for i in range(len(at_types)):
            print(" %4i %-2s -" % (i + 1, at_types[i]), end="")
            for j in range(len(bond_graphs[i])):
                print(" %i" % (bond_graphs[i][j] + 1), end="")
            print("\n", end="")
        print("\n", end="")

    def __init__(self, name: str, XYZFilePath: str):
        self.name = name

        # Load the XYZ File from XYZ File
        self.XYZCoordinates = self.ReadXYZ(path=XYZFilePath)
