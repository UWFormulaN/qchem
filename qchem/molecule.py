import dis
from math import pi
import math
import random
from uu import Error
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

    PositionSlice = slice(1, 4)
    """Constant that Grabs a Array/List of the Atom Position. Used for the Function GetAtomPosition. Is global so that you can get the position using self.XYZCoordinates.iloc[atomIndex, self.PositionSlice]"""

    Name: str = ""
    """Name of the Molecule"""

    #  Atom Symbol, X, Y, Z 
    XYZCoordinates: pd.core.frame.DataFrame
    """Data Frame of the XYZ Coordinates of the Atoms in the Molecule"""

    AtomCount: int = 0

    # Atom Index, Atom Symbol, Array of Index of other Atoms it's bonded to, Array of Booleans determining if Bond is Rotatable
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

    def GetRadius(self, atomIndex1: int, atomIndex2: int):
        """Gets the Radius Between 2 Atoms"""
        vector = self.GetAtomPosition(atomIndex1) - self.GetAtomPosition(atomIndex2)
        return np.linalg.norm(vector)


    def GetAllAtomsAfterBond(self, atomIndex1, atomIndex2) -> list[int]:
        """Starts a Recursive Search to find all the Atoms Present after a Bond"""

        if (atomIndex2 not in self.Bonds["Bonds"][atomIndex1]):
            Error("These Atoms are not Bonded Together")

        return self.BranchSearch(atomIndex2, [] , atomIndex1)

    def BranchSearch (self, currentIndex: int, atoms: list[int] = None,  ignoreIndex: int = None) -> list[int]:
        """Recursively Branch Searches through all the Atoms in a Molecule"""
        bonds: list[int] = self.Bonds["Bonds"][currentIndex]

        # Loop through all Bonds
        for i in bonds:

            # Ignore the Atom That we just came from
            if i == ignoreIndex:
                continue
            
            # If the Atom is 
            if i not in atoms:
                atoms.append(i)
                self.BranchSearch(i,  atoms, currentIndex)
          
        return atoms

    def FindRotatableBonds (self):
        """Goes through all Bonds and Determine if it's rotatable"""
        #TODO: Will need to Factor in Double Bonds in the Future

        rotatableBonds: list[list[bool]] =[]

        # Go through all Bonds in the Molecule. Do a Chain search 
        for i in range(self.AtomCount):

            # Cache the Bonds for the current Atom and initialize 
            bonds = self.Bonds["Bonds"][i]
            rotatableList :list[bool] = [] 

            # Loop through All Atoms that are Bonded
            for j in range(bonds):
                # Get the List of Atoms that come after the Bond of i - j -> 
                atoms: list[int] = self.GetAllAtomsAfterBond(i, j)

                # Rules that determine if the Bond is Rotatable TODO: Expand to factor in double bonds and more rules
                if i in atoms:
                    rotatableList.append(False)
                else:
                    rotatableList.append(True)

            # Adds the results of rotatable Bond 
            rotatableBonds.append(rotatableList)

        # Add to Data Frame
        self.Bonds["Rotatable"] = rotatableBonds
            
    def GetBonds(self):
        """Generates a Data Frame with all Bond related Information"""
        # Pre initialize variables
        at_types = self.XYZCoordinates["Atom"].values
        index = [i for i in range(self.AtomCount)]
        bonds = [[] for i in range(self.AtomCount)]
        bonds_distance = [[] for i in range(self.AtomCount)]

        # Get the Bonds and save to Array
        for i in range(self.AtomCount):
            radii1 = CovalentRadiiConstants[at_types[i]]
            for j in range(i + 1, self.AtomCount):
                radii2 = CovalentRadiiConstants[at_types[j]]
                thresh = 1.1 * (radii1 + radii2)
                dist = self.GetRadius(i, j)
                if dist < thresh:
                    bonds[i].append(j)
                    bonds[j].append(i)

        # Get Bond Distances
        for i in range(self.AtomCount):
            for j in bonds[i]:
                bonds_distance[i].append(self.GetRadius(i, j))

        # Save new Bonds Data Frame to Bonds Variable
        self.Bonds = pd.DataFrame(
            {
                "Index": index,
                "Atom": at_types,
                "Bonds": bonds,
                "Bond Distance": bonds_distance,
            }
        )

    def GetAtomPosition(self, atomIndex):
        """Returns a Numpy Array of the Atoms Position"""
        return np.array(self.XYZCoordinates.iloc[atomIndex, self.PositionSlice], dtype=float)

    def GetDihedralAngle(self, atomIndex1, atomIndex2, atomIndex3, atomIndex4):
        """Returns the Dihedral Angle of between Atoms"""
        # Get the position of the Atoms
        atom1Pos = self.GetAtomPosition(atomIndex1)
        atom2Pos = self.GetAtomPosition(atomIndex2)
        atom3Pos = self.GetAtomPosition(atomIndex3)
        atom4Pos = self.GetAtomPosition(atomIndex4)

        # Get the Vectors between each atom
        v21 = atom2Pos - atom1Pos
        v32 = atom3Pos - atom2Pos
        v43 = atom4Pos - atom3Pos

        # Calculate the Dihedral Angles
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
        if chi < -180.0:
            chi = chi + 360.0
        return chi

    def GetAngleBetweenAtoms(self, atomIndex1, atomIndex2, atomIndex3):
        """Returns the Angle between Atoms in Degrees"""
        # Get Atom Positions
        atom1Pos = self.GetAtomPosition(atomIndex1)
        atom2Pos = self.GetAtomPosition(atomIndex2)
        atom3Pos = self.GetAtomPosition(atomIndex3)

        # Get Normalized Vectors
        v1 = atom1Pos - atom2Pos
        v1 = v1 / np.linalg.norm(v1)
        v2 = atom3Pos - atom2Pos
        v2 = v2 / np.linalg.norm(v2)

        # Return the Angle between the Vectors aka between Atoms
        return (180 / pi) * math.acos(np.dot(v1, v2))

    def DisplayBondGraph(self):
        """Displays the Bond Graph in Terminal"""
        # Display Title Header
        print("   %s\n" % (self.name), end="")

        for i in range(self.AtomCount):
            # Get Index and Atom Symbol
            index = self.Bonds["Index"][i]
            atom = self.Bonds["Atom"][i]
            rotatable = self.Bonds["Rotatable"]

            # Create the String for Bonds
            bonds = ""
            for j in self.Bonds["Bonds"][i]:
                bonds += str(j + 1) + " "

            # Create String for Bond Distance
            bond_dist = ""
            for j in self.Bonds["Bond Distance"][i]:
                bond_dist += "%.3fÅ " % j

            # Create String for Rotatable Bonds
            rotatable = ""
            for j in self.Bonds["Rotatable"][i]:
                if j:
                    rotatable += "T "
                else:
                    rotatable += "F "

            # Print Line to Screen
            print(" %4i   %-2s - %s    %4s  %2s" % (index + 1, atom, bonds, bond_dist, rotatable))

    def GetDihedralAtomChain(self, atomChain: list[int], depth=0):
        """Recursively Searches for a Viable Atom Chain of length 4"""
        # Grab Bonds and Sort by incrementing index
        bonds: list[int] = self.Bonds["Bonds"][atomChain[depth]]
        bonds.sort()

        # Loop through the Sorted Bonds
        for i in bonds:

            # Return the Result if we are already at a length of 4
            if len(atomChain) == 4:
                continue

            # Skip if Atom is already in chain
            if i in atomChain:
                continue

            # Get this Atoms Bonds
            currentAtomBonds: list[int] = self.Bonds["Bonds"][i]

            # Skip if the Number of Bonds that the Current Atom has is Only one, Unless it is going to be the Last Atom in the Chain
            if len(currentAtomBonds) == 1 and (not len(atomChain) == 3):
                continue

            # Chosen to be a valid chain so we add it to the list
            atomChain.append(i)

            # Go one layer deeper in recursion to add an Atom to the chain
            self.GetDihedralAtomChain(atomChain, depth + 1)

        return atomChain

    def DisplayZMatrix(self):
        """Displays the Z Matrix in the Terminal"""
        z_matrix = self.CreateZMatrixFile()
        for i in z_matrix:
            print(i)

    def CreateZMatrixFile(self):
        """Creates an Array of Strings that Reprensents the ZMatrix File"""
        # Create an Array Reprensenting the File, First line is Number of Atoms, Second is the Name of the Molecule
        z_matrix = []
        z_matrix.append(f"{self.AtomCount}")
        z_matrix.append(self.name)

        # Loop through all the Atoms
        for i in range(self.AtomCount):

            # Create an Atom Chain starting with the Current Atom
            atomChain = self.GetDihedralAtomChain([i])

            # If the Chain hasn't reached 4 add Random Atoms until we reach 4
            while len(atomChain) < 4:
                randomAtom = random.randint(0, self.AtomCount - 1)
                if randomAtom not in atomChain:
                    atomChain.append(randomAtom)

            # Save the Indexes to Variables
            j = atomChain[1]
            k = atomChain[2]
            l = atomChain[3]

            # Cache the Z Matrix Variables
            atomSymbol = self.XYZCoordinates["Atom"][i]
            radius = self.GetRadius(i, j)
            angle = self.GetAngleBetweenAtoms(i, j, k)
            dihedralAngle = self.GetDihedralAngle(i, j, k, l)

            # Check if First Atom, Only Append symbol
            if i == 0:
                z_matrix.append(f"{atomSymbol}")
                continue

            # Check if Second Atom, Append Symbol and Distance from Next Closest Atom
            if i == 1:
                z_matrix.append(f"{atomSymbol} {j+1} {radius}")
                continue

            # Check if Third Atom, Append Symbol, Distance and Angle
            if i == 2:
                z_matrix.append(f"{atomSymbol} {j+1} {radius} {k+1} {angle}")
                continue

            # Append Symbol, Distance, Angle and Dihedral Angle
            z_matrix.append(
                f"{atomSymbol} {j+1} {radius} {k+1} {angle} {l+1} {dihedralAngle}"
            )

        return z_matrix

    def __init__(self, name: str, XYZFilePath: str):
        """Initializes a New Molecule Object"""
        self.name = name

        # Load the XYZ File from XYZ File
        self.XYZCoordinates = self.ReadXYZ(path=XYZFilePath)

        self.AtomCount = len(self.XYZCoordinates["Atom"].values)

        self.GetBonds()
        self.FindRotatableBonds()