from array import array
from ast import List
import dis
from math import pi, sqrt
import math
from os import read
import string
import pandas as pd
from scipy import optimize
from .Data.constants import CovalentRadiiConstants

# We will create molecule objects which will store information about the molecule
# Includes coordinates, atom types, how optimization was performed, how energy calculations were performed, etc.
# We can take segmented properties from the output file, and store them as attributes of the molecule object


class Molecule:
    """Class that represents an Entire Molecule. Stores information aboput Bonds, Atomic Positions unique Atom Properties and Allows for Specific Simulation Calculations"""

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

    def GetRadius(
        self, atom1: tuple[float, float, float], atom2: tuple[float, float, float]
    ):
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

    # def get_u12(coords1, coords2):
    #     r12 = get_r12(coords1, coords2)
    #     u12 = [0.0 for p in range(3)]
    #     for p in range(3):
    #         u12[p] = (coords2[p] - coords1[p]) / r12
    #     return u12

    def get_udp(uvec1, uvec2):
        udp = 0.0
        for p in range(3):
            udp += uvec1[p] * uvec2[p]
        udp = max(min(udp, 1.0), -1.0)
        return udp

    def GetMagnitude(self, vector):
        return sqrt( pow(vector[0], 2) + pow(vector[1], 2) + pow(vector[2], 2))

    def GetDotProductAngle(self, vector1, vector2):
        dotProduct = 0
        for i in range(len(vector1)):
            dotProduct += vector1[i] * vector2[i]
        return math.acos(
            dotProduct / (self.GetMagnitude(vector1) * self.GetMagnitude(vector2))
        )

    # # calculate unit cross product between two unit vectors
    # def get_ucp(uvec1, uvec2):
    #     ucp = [0.0 for p in range(3)]
    #     cos_12 = get_udp(uvec1, uvec2)
    #     sin_12 = math.sqrt(1 - cos_12**2)
    #     ucp[0] = (uvec1[1] * uvec2[2] - uvec1[2] * uvec2[1]) / sin_12
    #     ucp[1] = (uvec1[2] * uvec2[0] - uvec1[0] * uvec2[2]) / sin_12
    #     ucp[2] = (uvec1[0] * uvec2[1] - uvec1[1] * uvec2[0]) / sin_12
    #     return ucp

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
            print(" %4i   %-2s - %s          %4s" % (index + 1, atom, bonds, bond_dist))

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

            

    def CreateInternalFile (self):

        z_matrix = [None] * self.AtomCount
        z_matrix[0] = (self.XYZCoordinates["Atom"][0])
        z_matrix[1] = ((self.XYZCoordinates["Atom"][1], 1, self.GetRadius(self.GetPositionTuple(0), self.GetPositionTuple(1))))
        z_matrix[2] = ((self.XYZCoordinates["Atom"][2], 2, self.GetRadius(self.GetPositionTuple(1), self.GetPositionTuple(2)), 1,   self.GetAngleBetweenAtoms(0, 1, 2)))


        for i in range(3, self.AtomCount):
            for j in range(len(self.Bonds["Bonds"][i])):
                for k in range(len(self.Bonds["Bonds"][j])):

                    if (i == j or i == k or j == k):
                        continue
                    atom = self.XYZCoordinates["Atom"][k]
                    dist = self.GetRadius(self.GetPositionTuple(j), self.GetPositionTuple(k))
                    angle = self.GetAngleBetweenAtoms(i, j, k)
                    dihedral = self.GetAngleBetweenAtoms(i, j, k) # Fix this to dihedral angle

                    z_matrix[i] = (atom, k, dist, j, angle, i, dihedral)

                    #z_matrix.append((atom, k, dist, j, angle, i, dihedral))
                # z_matrix.append()

   
        # z_matrix.append((self.XYZCoordinates["Atom"][0],)) 
        # z_matrix.append((self.XYZCoordinates["Atom"][1], 1, self.GetRadius(self.GetPositionTuple(0), self.GetPositionTuple(1))))
        # z_matrix.append((self.XYZCoordinates["Atom"][2], 2, self.GetRadius(self.GetPositionTuple(1), self.GetPositionTuple(2)), 1,   self.GetAngleBetweenAtoms(0, 1, 2)))

        
        # for i in range(3, self.AtomCount):
        #     atom = self.XYZCoordinates["Atom"][i]
        #     dist = self.GetRadius(self.GetPositionTuple(i), self.GetPositionTuple(i-1))
        #     angle = self.GetAngleBetweenAtoms(i-2, i-1, i)
        #     dihedral = self.GetAngleBetweenAtoms(i-2, i-1, i) # Fix this to dihedral angle
        #    # z_matrix.append()

        #     # dist = calculate_distance(coordinates[i], coordinates[i-1])
        #     # angle = calculate_angle(coordinates[i-2], coordinates[i-1], coordinates[i])
        #     # dihedral = calculate_dihedral(coordinates[i-3], coordinates[i-2], coordinates[i-1], coordinates[i])
        #     z_matrix.append((atom, i, dist, i-1, angle, i-2, dihedral))
        
        return z_matrix


    def __init__(self, name: str, XYZFilePath: str):
        """Initializes a New Molecule Object"""
        self.name = name

        # Load the XYZ File from XYZ File
        self.XYZCoordinates = self.ReadXYZ(path=XYZFilePath)

        self.AtomCount = len(self.XYZCoordinates["Atom"].values)

        self.GetBonds()
