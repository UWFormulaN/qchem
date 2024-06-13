from array import array
from ast import List
import dis
from math import pi, sqrt
import math
from os import read
import random
import string
from unittest.util import sorted_list_difference
import pandas as pd
import numpy as np
from scipy import optimize

from .Data.constants import CovalentRadiiConstants

# We will create molecule objects which will store information about the molecule
# Includes coordinates, atom types, how optimization was performed, how energy calculations were performed, etc.
# We can take segmented properties from the output file, and store them as attributes of the molecule object


# May need to Redesign all the math to use dependency injection to a Class called Internal Molecular Math

class Molecule:
    """Class that represents an Entire Molecule. Stores information aboput Bonds, Atomic Positions unique Atom Properties and Allows for Specific Simulation Calculations"""

    # The Following are Variables that don't need to be initialized immediately or specified by the user, but belong to the Molecule class
    # For Some reason the Description needs to be defined after?

    PositionSlice = slice(1, 4)

    Name: str = ""
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

    def GetRadius(
        self, atom1: tuple[float, float, float], atom2: tuple[float, float, float]
    ):
        """Gets the Radius Between 2 Atoms"""
        return sqrt(
            pow(atom1[0] - atom2[0], 2)
            + pow(atom1[1] - atom2[1], 2)
            + pow(atom1[2] - atom2[2], 2)
        )
    
    def GetRadiusByIndex(
        self, atomIndex1: int, atomIndex2: int
    ):
        atom1 = self.GetPositionTuple(atomIndex1)
        atom2 = self.GetPositionTuple(atomIndex2)

        """Gets the Radius Between 2 Atoms"""
        return sqrt(
            pow(atom1[0] - atom2[0], 2)
            + pow(atom1[1] - atom2[1], 2)
            + pow(atom1[2] - atom2[2], 2)
        )

    def GetBonds(self):
        """Generates a Data Frame with all Bond related Information"""

        # Pre initialize variables
        at_types = self.XYZCoordinates["Atom"].values
        index = [i for i in range(self.AtomCount)]
        bonds = [[] for i in range(self.AtomCount)]
        bonds_distance = [[] for i in range(self.AtomCount)]
        coords = [
            self.GetPositionTuple(i)
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

        # Get Bond Distances
        for i in range(self.AtomCount):
            for j in bonds[i]:
                bonds_distance[i].append(self.GetRadius(coords[i], coords[j]))

        # Save new Bonds Data Frame to Bonds Variable
        self.Bonds = pd.DataFrame(
            {
                "Index": index,
                "Atom": at_types,
                "Bonds": bonds,
                "Bond Distance": bonds_distance,
            }
        )

    def GetPositionTuple(self, atomIndex):
        return (
            self.XYZCoordinates["X"][atomIndex],
            self.XYZCoordinates["Y"][atomIndex],
            self.XYZCoordinates["Z"][atomIndex],
        )

    def GetMagnitude(self, vector):
        return sqrt( pow(vector[0], 2) + pow(vector[1], 2) + pow(vector[2], 2))

    def GetDotProductAngle(self, vector1, vector2):
        dotProduct = 0
        for i in range(len(vector1)):
            dotProduct += vector1[i] * vector2[i]
        return math.acos(
            dotProduct / (self.GetMagnitude(vector1) * self.GetMagnitude(vector2))
        )

    def GetAtomPosition(self, atomIndex):
        """Returns a Numpy Array of the Atoms Position"""
        return np.array(self.XYZCoordinates.iloc[atomIndex, self.PositionSlice], dtype=float)

    def GetDihedralAngle (self, atomIndex1, atomIndex2, atomIndex3, atomIndex4):
        atom1Pos = self.GetAtomPosition(atomIndex1)
        atom2Pos = self.GetAtomPosition(atomIndex2)
        atom3Pos = self.GetAtomPosition(atomIndex3)
        atom4Pos = self.GetAtomPosition(atomIndex4)

        v21 = atom2Pos - atom1Pos
        v32 = atom3Pos - atom2Pos
        v43 = atom4Pos - atom3Pos

        v1 = np.cross(v21, v32)
        v1 = v1 / np.linalg.norm(v1)
        v2 = np.cross(v43, v32)
        v2 = v2 / np.linalg.norm(v2)
        m1 = np.cross(v1, v32)
        m1 = m1 / np.linalg.norm(m1)
        x = np.dot(v1, v2)
        y = np.dot(m1, v2)
        chi = np.arctan2(y, x)
        chi = -180.0 - 180.0 * chi / np.pi
        if (chi < -180.0):
            chi = chi + 360.0

        return chi

    def GetVector(self, pos1: tuple[float, float, float], pos2: tuple[float, float, float]) -> tuple[float, float, float]: 
        return (pos2[0] - pos1[0], pos2[1] - pos1[1], pos2[2] - pos1[2])

    def GetAngleBetweenAtoms(self, atomIndex1, atomIndex2, atomIndex3):
        # Atom Index 2 is in the Middle  2 -> 1 and 2 -> 3
        vec1 = self.GetVector(self.GetPositionTuple(atomIndex2), self.GetPositionTuple(atomIndex1))
        vec2 = self.GetVector(self.GetPositionTuple(atomIndex2), self.GetPositionTuple(atomIndex3))
        return (180 / pi) * self.GetDotProductAngle(vec1, vec2)

    def PrintAngles(self):
        angles = []
        skipList = []
        for j in range(self.AtomCount):
            numOfBonds = len(self.Bonds["Bonds"][j])
            for a in range(numOfBonds):
                
                if (len(angles) > 0):
                    if skipList.__contains__(angles[:][0]): 
                        next()
                atomIndex = self.Bonds["Bonds"][j][a]
                # i = bond_graph[j][a]
                for b in range(a + 1, numOfBonds):
                    atomIndex2 = self.Bonds["Bonds"][j][b]
                    skipList.append(atomIndex2)
                    # k = bond_graph[j][b]
                    bondAngle = self.GetAngleBetweenAtoms(atomIndex, j, atomIndex2)
                    #a123 = get_a123(coords[i], coords[j], coords[k])
                    angles.append([atomIndex, j, atomIndex2, bondAngle])

        for i in range(len(angles)):
            print(angles[i])
        return angles

    def DisplayBondGraph(self):
        """Displays the Bond Graph in Terminal"""
        # Display Title Header
        print("   %s\n" % (self.name), end="")

        for i in range(self.AtomCount):
            # Get Index and Atom Symbol
            index = self.Bonds["Index"][i]
            atom = self.Bonds["Atom"][i]

            # Create the String for Bonds
            bonds = ""
            for j in self.Bonds["Bonds"][i]:
                bonds += str(j + 1) + " "

            # Create Distance for Bond Distance
            bond_dist = ""
            for j in self.Bonds["Bond Distance"][i]:
                bond_dist += "%.3fÅ " % j

            # Print Line to Screen
            print(" %4i   %-2s - %s    %4s" % (index + 1, atom, bonds, bond_dist))

    def RecursiveBondSearch (self, depth, atomIndex, previousVisits: list[int],  previousChain: tuple[int, int, int]):
        for i in self.Bonds["Bonds"][atomIndex]:
            if i in previousVisits:
                continue
            if depth >= self.AtomCount:
                return previousChain
           
            newChain = (previousChain[1], previousChain[2], i+1)
            previousVisits.append(i)
            print(newChain)
            result = self.RecursiveBondSearch( depth + 1, i, previousVisits, newChain )

            if result:
                return result

    def GetAtomChain(self, atomChain: list[int], depth = 0):

        # Grab Bonds and Sort by incrementing index
        bonds : list[int] = self.Bonds["Bonds"][atomChain[depth]]
        bonds.sort()

        # Loop through the Sorted Bonds
        for i in bonds:

            # Return the Result if we are already at a length of 4
            if (len(atomChain) == 4):
                continue

            # Skip if Atom is already in chain
            if i in atomChain:
                continue

            # Get this Atoms Bonds
            currentAtomBonds : list[int] = self.Bonds["Bonds"][i]

            # Skip if the Number of Bonds that the Current Atom has is Only one, Unless it is going to be the Last Atom in the Chain
            if (len(currentAtomBonds) == 1 and (not len(atomChain) == 3)):
                continue

            # Chosen to be a valid chain so we add it to the list
            atomChain.append(i)

            # Go one layer deeper in recursion to add an Atom to the chain
            self.GetAtomChain(atomChain, depth + 1)

        return atomChain

    def DisplayZMatrix(self):
        z_matrix = self.CreateZMatrixFile()
        for i in z_matrix:
            print(i)

    def CreateZMatrixFile (self):
        """Creates an Array of Strings that Reprensents the ZMatrix File"""

        # Create an Array Reprensenting the File, First line is Number of Atoms, Second is the Name of the Molecule
        z_matrix = []
        z_matrix.append(f"{self.AtomCount}")
        z_matrix.append(self.name)

        # Loop through all the Atoms
        for i in range(self.AtomCount):

            # Create an Atom Chain starting with the Current Atom
            atomChain = self.GetAtomChain([i])

            # If the Chain hasn't reached 4 add Random Atoms until we reach 4
            while (len(atomChain) < 4):
                randomAtom = random.randint(0, self.AtomCount - 1)
                if randomAtom not in atomChain:
                    atomChain.append(randomAtom)

            j = atomChain[1]
            k = atomChain[2]
            l = atomChain[3]

            atomSymbol = self.XYZCoordinates["Atom"][i]

            # Check if First Atom, Only Append symbol
            if (i == 0):
                z_matrix.append(f"{atomSymbol}")
                continue

            # Check if Second Atom, Append Symbol and Distance from Next Closest Atom
            if (i == 1):
                z_matrix.append(f"{atomSymbol} {j+1} {self.GetRadiusByIndex(i, j)}")
                continue

            # Check if Third Atom, Append Symbol, Distance and Angle
            if (i == 2):
                z_matrix.append(f"{atomSymbol} {j+1} {self.GetRadiusByIndex(i, j)} {k+1} {self.GetAngleBetweenAtoms(i, j, k)}")
                continue

            # Append Symbol, Distance, Angle and Dihedral Angle
            z_matrix.append(f"{atomSymbol} {j+1} {self.GetRadiusByIndex(i, j)} {k+1} {self.GetAngleBetweenAtoms(i, j, k)} {l+1} {self.GetDihedralAngle(i, j, k, l)}")

        return z_matrix
    

    def __init__(self, name: str, XYZFilePath: str):
        """Initializes a New Molecule Object"""
        self.name = name

        # Load the XYZ File from XYZ File
        self.XYZCoordinates = self.ReadXYZ(path=XYZFilePath)

        self.AtomCount = len(self.XYZCoordinates["Atom"].values)

        self.GetBonds()
