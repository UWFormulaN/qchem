from ast import List
import dis
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

    AtomCount: int = 0

    # Atom Index, Atom Symbol, Bond Index 1,Bond Index 2 , Bond Index 3, Bond Index 4
    Bonds: pd.core.frame.DataFrame

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

    def GetRadius(self, atom1, atom2):
        """Gets the Radius Between 2 Atoms """
        return sqrt(
            pow(atom1[0] - atom2[0], 2)
            + pow(atom1[1] - atom2[1], 2)
            + pow(atom1[2] - atom2[2], 2)
        )
    

    def GetBonds(self):

        at_types = self.XYZCoordinates["Atom"].values
        coords = [
            (
                self.XYZCoordinates["X"][i],
                self.XYZCoordinates["Y"][i],
                self.XYZCoordinates["Z"][i],
            )
            for i in range(0, len(self.XYZCoordinates["X"].values))
        ]
        bonds = [[] for i in range(self.AtomCount)]
        bonds_distance = [[] for i in range(self.AtomCount)]

        for i in range(self.AtomCount):
            radii1 = CovalentRadiiConstants[at_types[i]]
            for j in range(i + 1, self.AtomCount):
                radii2 = CovalentRadiiConstants[at_types[j]]
                thresh = 1.1 * (radii1 + radii2)
                dist = self.GetRadius(coords[i], coords[j])
                if dist < thresh:
                    bonds_distance[i].append(dist)
                    bonds_distance[j].append(dist)
                    bonds[i].append(j)
                    bonds[j].append(i)

        index = [i+1 for i in range(self.AtomCount)] 
        bondDataFrame = pd.DataFrame({
            "Index": index,
            "Atom" : at_types,
            "Bonds": bonds,
            "Bond Distance": bonds_distance
            
        })
        self.Bonds = bondDataFrame

        
    def DisplayBondGraph (self):

        print("%s\n" % (self.name), end="")

        for i in range(self.AtomCount):

            index = self.Bonds["Index"][i]

            atom = self.Bonds["Atom"][i]

            # Create the String for Bonds 
            bonds = ""
            for i in self.Bonds["Bonds"][i]:
                bonds += str(i) + " "

            print(" %4i   %-2s - %s" % (index, atom, bonds))

    def __init__(self, name: str, XYZFilePath: str):
        self.name = name

        # Load the XYZ File from XYZ File
        self.XYZCoordinates = self.ReadXYZ(path=XYZFilePath)

        self.AtomCount = len(self.XYZCoordinates["Atom"].values)

        self.GetBonds()
