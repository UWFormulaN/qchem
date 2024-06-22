import copy
from math import pi
import math
import random
from uu import Error
import pandas as pd
import numpy as np


from .Data.constants import CovalentRadiiConstants

# We will create molecule objects which will store information about the molecule
# Includes coordinates, atom types, how optimization was performed, how energy calculations were performed, etc.
# We can take segmented properties from the output file, and store them as attributes of the molecule object
# May need to Redesign all the math to use dependency injection to a Class called Internal Molecular Math

class Molecule:
    """Class that represents an Entire Molecule. Stores information aboput Bonds, Atomic Positions unique Atom Properties and Allows for Specific Simulation Calculations"""

    Name: str = ""
    """Name of the Molecule"""

    #  Atom Symbol, X, Y, Z 
    XYZCoordinates: pd.core.frame.DataFrame
    """Data Frame of the XYZ Coordinates of the Atoms in the Molecule"""

    AtomCount: int = 0
    """The Number of Atoms present in the Molecule"""

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

    def __init__(self, name: str, XYZFilePath: str):
        """Initializes a New Molecule Object"""
        self.Name = name

        # Load the XYZ File from XYZ File
        self.XYZCoordinates = self.ReadXYZ(path=XYZFilePath)
        self.AtomCount = len(self.XYZCoordinates["Atom"].values)

        self.GetBonds()
        self.FindRotatableBonds()

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
            raise Error("These Atoms are not Bonded Together")

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
    
    def CopyPositions (self) -> pd.core.frame.DataFrame:
        return copy.deepcopy(self.XYZCoordinates)

    def Copy (self):
        return copy.deepcopy(self)
    
    def GeneratePerpendicularVector(self, v):
        # Check if the input vector is not zero
        if np.linalg.norm(v) == 0:
            raise ValueError("The input vector cannot be the zero vector")
        
        # Create a second vector which is not parallel to the input vector
        if v[0] == 0:
            w = np.array([1, 0, 0])  # Choose a simple vector if the first element is zero
        else:
            w = np.array([v[1], -v[0], 0])  # This ensures the vector is not parallel

        # Compute the cross product
        cross_product = np.cross(v, w)
        
        # Normalize the cross product
        norm = np.linalg.norm(cross_product)
        if norm == 0:
            raise ValueError("The cross product resulted in a zero vector, which should not happen")
        
        normalized_vector = cross_product / norm
        
        return normalized_vector

    def RotateBond (self, atomIndex1: int, atomIndex2: int, radians: float):
        """Rotates a Bond in the Molecule by the Specified Radians"""

        atomIndexes: list[int] = self.GetAllAtomsAfterBond(atomIndex1, atomIndex2)
        atomPositions: np.ndarray = np.array(self.GetAtomPosition(atomIndexes[0]) - self.GetAtomPosition(atomIndex2))

        for i in range(1, len(atomIndexes)):
            atomPositions = np.vstack((atomPositions, self.GetAtomPosition(atomIndexes[i]) - self.GetAtomPosition(atomIndex2)))

        z_vector = self.GetAtomPosition(atomIndex2) - self.GetAtomPosition(atomIndex1)
        z_vector = z_vector / np.linalg.norm(z_vector)

        x_vector = self.GeneratePerpendicularVector(z_vector)
        
        y_vector = np.cross(z_vector, x_vector)
        y_vector = y_vector / np.linalg.norm(y_vector)

        original_axes = np.array([
            x_vector,  # Original X-axis
            y_vector,  # Original Y-axis
            z_vector   # Original Z-axis
        ])

        # Define rotation matrix around Z-axis by given radians
        cos_theta = np.cos(radians)
        sin_theta = np.sin(radians)
        z_rotation_matrix = np.array([
            [cos_theta, -sin_theta, 0],
            [sin_theta,  cos_theta, 0],
            [0,          0,         1]
        ])

        # Combine the transformations to find the full rotation matrix
        rotation_basis = original_axes.T
        rotation_matrix = rotation_basis @ z_rotation_matrix @ np.linalg.inv(rotation_basis)

        # Apply the rotation to all atom positions
        rotatedAtoms = atomPositions @ rotation_matrix.T

        # Update atom coordinates in the original structure
        for i, atomIndex in enumerate(atomIndexes):
            newPosition = rotatedAtoms[i] + self.GetAtomPosition(atomIndex2)
            self.XYZCoordinates.loc[atomIndex, "X"] = newPosition[0]
            self.XYZCoordinates.loc[atomIndex, "Y"] = newPosition[1]
            self.XYZCoordinates.loc[atomIndex, "Z"] = newPosition[2]

    def GetConformers (self, atomIndex1, atomIndex2, steps):
        """Returns a list of Conformers Molecules where a bond is rotated by (2pi / steps) radians """
        stepSizeRad = (2 * pi)/steps
        conformers: list[Molecule] = []

        # Loops and creates a new Rotated Molecule to add to the List of Conformers
        for i in range(steps):
            rotation = stepSizeRad * i
            newMolecule = self.Copy()
            newMolecule.Name = f"{self.Name}_rot_{(180 / pi) * rotation}"
            newMolecule.RotateBond(atomIndex1, atomIndex2, rotation)
            conformers.append(newMolecule)

        return conformers

    def DisplayXYZ (self):
        """Prints the XYZ File to the Terminal"""
        from qchem.XYZFile import XYZFile
        print(XYZFile(self).GetFileAsString())

    def SaveAsXYZ (self, filePath):
        """Saves the Molecule to a XYZ File"""
        from qchem.XYZFile import XYZFile
        XYZFile(self).SaveToFile(filePath)

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
            for j in bonds:
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
        return np.array(self.XYZCoordinates.iloc[atomIndex, slice(1, 4)], dtype=float)

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
        print("   %s\n" % (self.Name), end="")

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
        z_matrix.append(self.Name)

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