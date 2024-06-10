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
    """Data Frame of the XYZ Coordinates of the Atoms in the Molecule"""

    AtomCount: int = 0

    # Atom Index, Atom Symbol 
    Bonds: pd.core.frame.DataFrame
    """Data Frame of the Bonds connected to each Atom, as well as the Bond Lengths"""

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
        """Reads the XYZ file and Giving a Data Frame with the Position of Each Atom
        
        Returns : Pandas Data Frame with Columns for Atom Symbol and X, Y, Z Position
        """
        return pd.read_csv(
            path, sep=r"\s+", skiprows=2, names=["Atom", "X", "Y", "Z"], engine="python"
        )

    def GetGeometry(self):
        """Displays the Geometry of the Molecule in the Terminal"""
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
        """Generates a Data Frame with all Bond related Information"""

        # Pre initialize most variables
        at_types = self.XYZCoordinates["Atom"].values
        index = [i for i in range(self.AtomCount)] 
        bonds = [[] for i in range(self.AtomCount)]
        bonds_distance = [[] for i in range(self.AtomCount)]
        coords = [
            (
                self.XYZCoordinates["X"][i],
                self.XYZCoordinates["Y"][i],
                self.XYZCoordinates["Z"][i],
            )
            for i in range(0, len(self.XYZCoordinates["X"].values))
        ]

        # Get the Bonds and save to Array
        for i in range(self.AtomCount):
            radii1 = CovalentRadiiConstants[at_types[i]]
            for j in range(i + 1, self.AtomCount):
                radii2 = CovalentRadiiConstants[at_types[j]]
                thresh = 1.1 * (radii1 + radii2)
                dist = self.GetRadius(coords[i], coords[j])
                if dist < thresh:
                    bonds[i].append(j)
                    bonds[j].append(i)

        #Get Bond Distances
        for i in range(self.AtomCount):
            for j in bonds[i]:
                bonds_distance[i].append(self.GetRadius(coords[i], coords[j]))
    
        # Save new Bonds Data Frame to Bonds Variable
        self.Bonds = pd.DataFrame({
            "Index": index,
            "Atom" : at_types,
            "Bonds": bonds,
            "Bond Distance": bonds_distance
        })

        
    def DisplayBondGraph (self):
        """Displays the Bond Graph in Terminal"""

        print("   %s\n" % (self.name), end="")

        for i in range(self.AtomCount):
            # Get Index and Atom Symbol
            index = self.Bonds["Index"][i]
            atom = self.Bonds["Atom"][i]

            # Create the String for Bonds 
            bonds = ""
            for j in self.Bonds["Bonds"][i]:
                bonds += str(j+1) + " "

            # Create Distance for Bond Distance
            bond_dist = ""
            for j in self.Bonds["Bond Distance"][i]:
                bond_dist += "%.3fÅ " % j 

            print(" %4i   %-2s - %s          %4s" % (index + 1, atom, bonds, bond_dist))

    def __init__(self, name: str, XYZFilePath: str):
        """Initializes a New Molecule Object"""
        self.name = name

        # Load the XYZ File from XYZ File
        self.XYZCoordinates = self.ReadXYZ(path=XYZFilePath)

        self.AtomCount = len(self.XYZCoordinates["Atom"].values)

        self.GetBonds()
