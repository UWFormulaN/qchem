import copy
from math import pi
import math
import random
import pandas as pd
import numpy as np
from .Data.Constants import AtomicMassConstants, CovalentRadiiConstants

# We will create molecule objects which will store information about the molecule
# Includes coordinates, atom types, how optimization was performed, how energy calculations were performed, etc.
# We can take segmented properties from the output file, and store them as attributes of the molecule object
# May need to Redesign all the math to use dependency injection to a Class called Internal Molecular Math

class Molecule:
    """Class that represents an Entire Molecule. Stores information aboput Bonds, Atomic Positions unique Atom Properties and Allows for Specific Simulation Calculations"""

    name: str = ""
    """Name of the Molecule"""

    #  Atom Symbol, X, Y, Z 
    XYZCoordinates: pd.core.frame.DataFrame
    """Data Frame of the XYZ Coordinates of the Atoms in the Molecule"""

    atomCount: int = 0
    """The Number of Atoms present in the Molecule"""

    # Atom Index, Atom Symbol, Array of Index of other Atoms it's bonded to, Array of Booleans determining if Bond is Rotatable
    bonds: pd.core.frame.DataFrame
    """Data Frame of the Bonds connected to each Atom, as well as the Bond Lengths"""

    energy = None
    """The Energy of the Molecule"""

    multiplicity = None
    """Multiplicity of the Molecule"""

    atomTypes = None
    """Atom types of the Molecule"""

    optimizationMethod = None
    """Optimization Method of the Molecule"""

    energyMethod = None
    """Energy Method used for the Molecule"""

    def __init__(self, name: str, XYZ):
        """Initializes a New Molecule Object\n
            name : str
            XYZ : str | XYZFile
        """
        from qchem.XYZFile import XYZFile # <-- To Avoid Circular Imports
        
        if (not isinstance(name, (str))):
            raise ValueError("Name must be a String")
        
        if (name == ""):
            raise ValueError("Name cannot be Empty")
        
        self.name = name
        
        if (isinstance(XYZ, (str))):
            # Load the XYZ File from XYZ File
            self.XYZCoordinates = self.sortAtomDataFrame(self.readXYZ(path=XYZ))
            self.atomCount = len(self.XYZCoordinates["Atom"].values)
        elif (isinstance(XYZ, (XYZFile))):
            # Load from the XYZ File Class
            self.atomCount = XYZ.atomCount
            self.XYZCoordinates = XYZ.atomPositions

        self.getBonds()
        self.findRotatableBonds()

    def readXYZ(self, path: str) -> pd.core.frame.DataFrame:
        """Reads the XYZ file and Giving a Data Frame with the Position of Each Atom

        Returns : Pandas Data Frame with Columns for Atom Symbol and X, Y, Z Position
        """
        return pd.read_csv(
            path, sep=r"\s+", skiprows=2, names=["Atom", "X", "Y", "Z"], engine="python"
        )

    def getGeometry(self):
        """Displays the Geometry of the Molecule in the Terminal"""
        for atom in self.XYZCoordinates.itertuples():
            print(atom)

    def getRadius(self, atomIndex1: int, atomIndex2: int):
        """Gets the Radius Between 2 Atoms"""
        vector = self.getAtomPosition(atomIndex1) - self.getAtomPosition(atomIndex2)
        return np.linalg.norm(vector)

    def getAllAtomsAfterBond(self, atomIndex1, atomIndex2) -> list[int]:
        """Starts a Recursive Search to find all the Atoms Present after a Bond"""
        if (atomIndex2 not in self.bonds["Bonds"][atomIndex1]):
            raise Exception("These Atoms are not Bonded Together")

        return self.branchSearch(atomIndex2, [] , atomIndex1)

    def branchSearch (self, currentIndex: int, atoms: list[int] = None,  ignoreIndex: int = None) -> list[int]:
        """Recursively Branch Searches through all the Atoms in a Molecule"""
        bonds: list[int] = self.bonds["Bonds"][currentIndex]

        # Loop through all Bonds
        for i in bonds:

            # Ignore the Atom That we just came from
            if i == ignoreIndex:
                continue
            
            # If the Atom is 
            if i not in atoms:
                atoms.append(i)
                self.branchSearch(i,  atoms, currentIndex)
          
        return atoms
    
    def copyPositions (self) -> pd.core.frame.DataFrame:
        return copy.deepcopy(self.XYZCoordinates)

    def copy (self):
        return copy.deepcopy(self)
    
    def generatePerpendicularVector(self, v):
        # Check if the input vector is not zero
        if np.linalg.norm(v) == 0:
            raise ValueError("The input vector cannot be the zero vector")
        
        # Create a second vector which is not parallel to the input vector
        if v[0] == 0:
            w = np.array([1, 0, 0])  # Choose a simple vector if the first element is zero
        else:
            w = np.array([v[1], -v[0], 0])  # This ensures the vector is not parallel

        # Compute the cross product
        crossProduct = np.cross(v, w)
        
        # Normalize the cross product
        norm = np.linalg.norm(crossProduct)
        if norm == 0:
            raise ValueError("The cross product resulted in a zero vector, which should not happen")
        
        normalizedVector = crossProduct / norm
        
        return normalizedVector

    def rotateBond (self, atomIndex1: int, atomIndex2: int, radians: float):
        """Rotates a Bond in the Molecule by the Specified Radians"""

        atomIndexes: list[int] = self.getAllAtomsAfterBond(atomIndex1, atomIndex2)
        atomPositions: np.ndarray = np.array(self.getAtomPosition(atomIndexes[0]) - self.getAtomPosition(atomIndex2))

        for i in range(1, len(atomIndexes)):
            atomPositions = np.vstack((atomPositions, self.getAtomPosition(atomIndexes[i]) - self.getAtomPosition(atomIndex2)))

        zVector = self.getAtomPosition(atomIndex2) - self.getAtomPosition(atomIndex1)
        zVector = zVector / np.linalg.norm(zVector)

        xVector = self.generatePerpendicularVector(zVector)
        
        yVector = np.cross(zVector, xVector)
        yVector = yVector / np.linalg.norm(yVector)

        originalAxes = np.array([
            xVector,  # Original X-axis
            yVector,  # Original Y-axis
            zVector   # Original Z-axis
        ])

        # Define rotation matrix around Z-axis by given radians
        cosTheta = np.cos(radians)
        sinTheta = np.sin(radians)
        zRotationMatrix = np.array([
            [cosTheta, -sinTheta, 0],
            [sinTheta,  cosTheta, 0],
            [0,          0,         1]
        ])

        # Combine the transformations to find the full rotation matrix
        rotationBasis = originalAxes.T
        rotationMatrix = rotationBasis @ zRotationMatrix @ np.linalg.inv(rotationBasis)

        # Apply the rotation to all atom positions
        rotatedAtoms = atomPositions @ rotationMatrix.T

        # Update atom coordinates in the original structure
        for i, atomIndex in enumerate(atomIndexes):
            newPosition = rotatedAtoms[i] + self.getAtomPosition(atomIndex2)
            self.XYZCoordinates.loc[atomIndex, "X"] = newPosition[0]
            self.XYZCoordinates.loc[atomIndex, "Y"] = newPosition[1]
            self.XYZCoordinates.loc[atomIndex, "Z"] = newPosition[2]

    def getConformers (self, atomIndex1, atomIndex2, steps):
        """Returns a list of Conformers Molecules where a bond is rotated by (2pi / steps) radians """
        stepSizeRad = (2 * pi)/steps
        conformers: list[Molecule] = []

        # Loops and creates a new Rotated Molecule to add to the List of Conformers
        for i in range(steps):
            rotation = stepSizeRad * i
            newMolecule = self.Copy()
            newMolecule.name = f"{self.Name}_rot_{(180 / pi) * rotation}"
            newMolecule.rotateBond(atomIndex1, atomIndex2, rotation)
            conformers.append(newMolecule)

        return conformers

    def XYZ (self):
        """Prints the XYZ File to the Terminal"""
        from qchem.XYZFile import XYZFile
        return XYZFile(self).getFileAsString()
        
    def XYZBody (self) -> str:
        """Prints the XYZ File to the Terminal"""
        from qchem.XYZFile import XYZFile
        return XYZFile(self).getXYZBody()

    def saveAsXYZ (self, filePath):
        """Saves the Molecule to a XYZ File"""
        from qchem.XYZFile import XYZFile
        XYZFile(self).saveToFile(filePath)

    def findRotatableBonds (self):
        """Goes through all Bonds and Determine if it's rotatable"""
        #TODO: Will need to Factor in Double Bonds in the Future

        rotatableBonds: list[list[bool]] =[]

        # Go through all Bonds in the Molecule. Do a Chain search 
        for i in range(self.atomCount):

            # Cache the Bonds for the current Atom and initialize 
            bonds = self.bonds["Bonds"][i]
            rotatableList :list[bool] = [] 

            # Loop through All Atoms that are Bonded
            for j in bonds:
                # Get the List of Atoms that come after the Bond of i - j -> 
                atoms: list[int] = self.getAllAtomsAfterBond(i, j)

                # Rules that determine if the Bond is Rotatable TODO: Expand to factor in double bonds and more rules
                if i in atoms:
                    rotatableList.append(False)
                else:
                    rotatableList.append(True)

            # Adds the results of rotatable Bond 
            rotatableBonds.append(rotatableList)

        # Add to Data Frame
        self.bonds["Rotatable"] = rotatableBonds
            
    def getBonds(self):
        """Generates a Data Frame with all Bond related Information"""
        # Pre initialize variables
        atTypes = self.XYZCoordinates["Atom"].values
        index = [i for i in range(self.atomCount)]
        bonds = [[] for i in range(self.atomCount)]
        bondsDistance = [[] for i in range(self.atomCount)]

        # Get the Bonds and save to Array
        for i in range(self.atomCount):
            radii1 = CovalentRadiiConstants[atTypes[i]]
            for j in range(i + 1, self.atomCount):
                radii2 = CovalentRadiiConstants[atTypes[j]]
                thresh = 1.1 * (radii1 + radii2)
                dist = self.getRadius(i, j)
                if dist < thresh:
                    bonds[i].append(j)
                    bonds[j].append(i)

        # Get Bond Distances
        for i in range(self.atomCount):
            for j in bonds[i]:
                bondsDistance[i].append(self.getRadius(i, j))

        # Save new Bonds Data Frame to Bonds Variable
        self.bonds = pd.DataFrame(
            {
                "Index": index,
                "Atom": atTypes,
                "Bonds": bonds,
                "Bond Distance": bondsDistance,
            }
        )

    def getAtomPosition(self, atomIndex):
        """Returns a Numpy Array of the Atoms Position"""
        return np.array(self.XYZCoordinates.iloc[atomIndex, slice(1, 4)], dtype=float)

    def getDihedralAngle(self, atomIndex1, atomIndex2, atomIndex3, atomIndex4):
        """Returns the Dihedral Angle of between Atoms"""
        # Get the position of the Atoms
        atom1Pos = self.getAtomPosition(atomIndex1)
        atom2Pos = self.getAtomPosition(atomIndex2)
        atom3Pos = self.getAtomPosition(atomIndex3)
        atom4Pos = self.getAtomPosition(atomIndex4)

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

    def getAngleBetweenAtoms(self, atomIndex1, atomIndex2, atomIndex3):
        """Returns the Angle between Atoms in Degrees"""
        # Get Atom Positions
        atom1Pos = self.getAtomPosition(atomIndex1)
        atom2Pos = self.getAtomPosition(atomIndex2)
        atom3Pos = self.getAtomPosition(atomIndex3)

        # Get Normalized Vectors
        v1 = atom1Pos - atom2Pos
        v1 = v1 / np.linalg.norm(v1)
        v2 = atom3Pos - atom2Pos
        v2 = v2 / np.linalg.norm(v2)

        # Return the Angle between the Vectors aka between Atoms
        return (180 / pi) * math.acos(np.dot(v1, v2))

    def displayBondGraph(self):
        """Displays the Bond Graph in Terminal"""
        # Display Title Header
        print("   %s\n" % (self.Name), end="")

        for i in range(self.atomCount):
            # Get Index and Atom Symbol
            index = self.bonds["Index"][i]
            atom = self.bonds["Atom"][i]
            rotatable = self.bonds["Rotatable"]

            # Create the String for Bonds
            bonds = ""
            for j in self.bonds["Bonds"][i]:
                bonds += str(j + 1) + " "

            # Create String for Bond Distance
            bondDist = ""
            for j in self.bonds["Bond Distance"][i]:
                bondDist += "%.3fÅ " % j

            # Create String for Rotatable Bonds
            rotatable = ""
            for j in self.bonds["Rotatable"][i]:
                if j:
                    rotatable += "T "
                else:
                    rotatable += "F "

            # Print Line to Screen
            print(" %4i   %-2s - %s    %4s  %2s" % (index + 1, atom, bonds, bondDist, rotatable))

    def getDihedralAtomChain(self, atomChain: list[int], depth=0):
        """Recursively Searches for a Viable Atom Chain of length 4"""
        # Grab Bonds and Sort by incrementing index
        bonds: list[int] = self.bonds["Bonds"][atomChain[depth]]
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
            currentAtomBonds: list[int] = self.bonds["Bonds"][i]

            # Skip if the Number of Bonds that the Current Atom has is Only one, Unless it is going to be the Last Atom in the Chain
            if len(currentAtomBonds) == 1 and (not len(atomChain) == 3):
                continue

            # Chosen to be a valid chain so we add it to the list
            atomChain.append(i)

            # Go one layer deeper in recursion to add an Atom to the chain
            self.getDihedralAtomChain(atomChain, depth + 1)

        return atomChain

    def displayZMatrix(self):
        """Displays the Z Matrix in the Terminal"""
        zMatrix = self.createZMatrixFile()
        for i in zMatrix:
            print(i)

    def createZMatrixFile(self):
        """Creates an Array of Strings that Reprensents the ZMatrix File"""
        # Create an Array Reprensenting the File, First line is Number of Atoms, Second is the Name of the Molecule
        zMatrix = []
        zMatrix.append(f"{self.atomCount}")
        zMatrix.append(self.name)

        # Loop through all the Atoms
        for i in range(self.atomCount):

            # Create an Atom Chain starting with the Current Atom
            atomChain = self.getDihedralAtomChain([i])

            # If the Chain hasn't reached 4 add Random Atoms until we reach 4
            while len(atomChain) < 4:
                randomAtom = random.randint(0, self.atomCount - 1)
                if randomAtom not in atomChain:
                    atomChain.append(randomAtom)

            # Save the Indexes to Variables
            j = atomChain[1]
            k = atomChain[2]
            l = atomChain[3]

            # Cache the Z Matrix Variables
            atomSymbol = self.XYZCoordinates["Atom"][i]
            radius = self.getRadius(i, j)
            angle = self.getAngleBetweenAtoms(i, j, k)
            dihedralAngle = self.getDihedralAngle(i, j, k, l)

            # Check if First Atom, Only Append symbol
            if i == 0:
                zMatrix.append(f"{atomSymbol}")
                continue

            # Check if Second Atom, Append Symbol and Distance from Next Closest Atom
            if i == 1:
                zMatrix.append(f"{atomSymbol} {j+1} {radius}")
                continue

            # Check if Third Atom, Append Symbol, Distance and Angle
            if i == 2:
                zMatrix.append(f"{atomSymbol} {j+1} {radius} {k+1} {angle}")
                continue

            # Append Symbol, Distance, Angle and Dihedral Angle
            zMatrix.append(
                f"{atomSymbol} {j+1} {radius} {k+1} {angle} {l+1} {dihedralAngle}"
            )

        return zMatrix
    
    def sortAtomDataFrame (self, dataFrame : pd.core.frame.DataFrame):
        """Sorts the Molecule Data Frame when initializing the Molecule to make sure the Heaviest Atoms are at the Top of the Data Frame"""
        # Create a duplicate of the DataFrame
        sortedDataFrame = dataFrame.copy()

        atoms = []

        for i in range(len(dataFrame["Atom"].values)):
            atoms.append([i, AtomicMassConstants[dataFrame["Atom"][i]]])

        # Sort atoms by atomic mass in descending order
        sortedAtoms = sorted(atoms, key=lambda x: x[1], reverse=True)

        # Extract the sorted indices from the sorted_atoms list
        sortedIndices = [atom[0] for atom in sortedAtoms]

        # Use the sorted indices to reindex the original DataFrame
        sortedDataFrame = dataFrame.iloc[sortedIndices]

        # Replace rows in the duplicate DataFrame using the sorted indices
        for newIndex, (originalIndex, _) in enumerate(sortedAtoms):
            sortedDataFrame.iloc[newIndex] = dataFrame.iloc[originalIndex]

        # Resets the Now Sorted indexes so that they properly increment
        sortedDataFrame.reset_index(drop=True, inplace=True)

        return sortedDataFrame
    
    def getMolecularWeight (self) -> float:
        """Returns the Molecular Weight of the Molecule in g/mol"""
        mw = 0

        # Loop through all Atoms and Add their Individual Atomic Mass
        for i in range(self.atomCount):
            mw += AtomicMassConstants[self.XYZCoordinates["Atom"][i]]

        return mw
