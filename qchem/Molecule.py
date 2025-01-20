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
    """Object Representing a Molecule. Exposes basic information and methods for the Molecule"""

    name: str = ""
    """Name of the Molecule (Often the Scientific name)"""

    #  Atom Symbol, X, Y, Z
    XYZCoordinates: pd.DataFrame
    """Pandas DataFrame storing each atoms Atomic Symbol and XYZ position in rows"""

    atomCount: int = 0
    """Number of Atoms present in the Molecule"""

    # Atom Index, Atom Symbol, Array of Index of other Atoms it's bonded to, Array of Booleans determining if Bond is Rotatable
    bonds: pd.DataFrame
    """Pandas DataFrame storing the index of each Atom, which Atoms they are Bonded to, their Bond Distances and if the Bond is rotatable """

    energy: float = None
    """Energy of the Molecule, the sum of the bond energies and the energy needed in it’s conformation"""

    # multiplicity: int = None
    # """Multiplicity of the Molecule"""
    #
    # atomTypes = None
    # """Atom types of the Molecule"""
    #
    # optimizationMethod = None
    # """Optimization Method of the Molecule"""
    #
    # energyMethod = None
    # """Energy Method used for the Molecule"""

    def __init__(self, name: str, XYZ):
        """Initializes a New Molecule Object\n
        name : str
        XYZ : str | XYZFile
        """
        from qchem.XYZFile import XYZFile  # <-- To Avoid Circular Imports

        if not isinstance(name, (str)):
            raise ValueError("Name must be a String")

        if name == "":
            raise ValueError("Name cannot be Empty")

        self.name = name

        if isinstance(XYZ, (str)):
            # Load the XYZ File from XYZ File
            self.XYZCoordinates = self.sortAtomDataFrame(self.readXYZ(path=XYZ))
            self.atomCount = len(self.XYZCoordinates["Atom"].values)
        elif isinstance(XYZ, (XYZFile)):
            # Load from the XYZ File Class
            self.atomCount = XYZ.atomCount
            self.XYZCoordinates = XYZ.atomPositions

        self.getBonds()
        self.findRotatableBonds()

    def readXYZ(self, path: str) -> pd.DataFrame:
        """Reads the Provided XYZ Files and returns the XYZ Format in a DataFrame

        ## Parameters : \n
            self : Molecule - Default Parameter for the Class Instance
            path : str - Path to the XYZ File

        ## Returns : \n
            pd.DataFrame - Pandas DataFrame with each row being an Atoms XYZ Position
        """
        return pd.read_csv(
            path, sep=r"\s+", skiprows=2, names=["Atom", "X", "Y", "Z"], engine="python"
        )

    def getGeometry(self):
        """Displays the Molecules Geometry in the Terminal

        ## Parameters : \n
            self : Molecule - Default Parameter for the Class Instance

        ## Returns : \n
            None - No Return Value
        """
        for atom in self.XYZCoordinates.itertuples():
            print(atom)

    def getRadius(self, atomIndex1: int, atomIndex2: int):
        """Gets the Radius between 2 Atoms in Angstroms

        ## Parameters : \n
            self : Molecule - Default Parameter for the Class Instance
            atomIndex1 : int - Index of the First Atom
            atomIndex2 : int - Index of the Second Atom

        ## Returns : \n
            float - Distance between the 2 Atoms in Angstroms
        """
        vector = self.getAtomPosition(atomIndex1) - self.getAtomPosition(atomIndex2)
        return np.linalg.norm(vector)

    def getAllAtomsAfterBond(self, atomIndex1: int, atomIndex2: int) -> list[int]:
        """Recusrsively Searches through all Atoms in the Molecule present after the specified bond

        ## Parameters : \n
            self : Molecule - Default Parameter for the Class Instance
            atomIndex1 : int - Index of the First Atom
            atomIndex2 : int - Index of the Second Atom

        ## Returns : \n
            list[int] - List of the index of each Atom present after the bond
        """
        if atomIndex2 not in self.bonds["Bonds"][atomIndex1]:
            raise Exception("These Atoms are not Bonded Together")

        return self.branchSearch(atomIndex2, [], atomIndex1)

    def branchSearch(
        self, currentIndex: int, atoms: list[int] = None, ignoreIndex: int = None
    ) -> list[int]:
        """Recursively Searches through all Atoms in the Molecule

        ## Parameters : \n
            self : Molecule - Default Parameter for the Class Instance
            currentIndex : int - Index of the Current Atom in the recursive search
            atoms : list[int] - List of Atoms previously visited in the Branch Search
            ignoreIndex : int - Index of the Molecule to ignore in the Search (Often just the Current index

        ## Returns : \n
            list[int] - List of the Atoms in the Molecule present after the specified bond
        """
        bonds: list[int] = self.bonds["Bonds"][currentIndex]

        # Loop through all Bonds
        for i in bonds:

            # Ignore the Atom That we just came from
            if i == ignoreIndex:
                continue

            # If the Atom is
            if i not in atoms:
                atoms.append(i)
                self.branchSearch(i, atoms, currentIndex)

        return atoms

    def copyPositions(self) -> pd.DataFrame:
        """Creates a Deep Copy of the XYZ Coordinates DataFrame

        ## Parameters : \n
            self : Molecule - Default Parameter for the Class Instance

        ## Returns : \n
            pd.DataFrame - Deep Copy of the XYZ Coordinates DataFrame
        """
        return copy.deepcopy(self.XYZCoordinates)

    def copy(self):
        """Creates a Deep Copy of the Molecule Object

        ## Parameters : \n
            self : Molecule - Default Parameter for the Class Instance

        ## Returns : \n
            Molecule - Deep Copy of the Molecule Object
        """
        return copy.deepcopy(self)

    def generatePerpendicularVector(self, vector: np.ndarray) -> np.ndarray:
        """Generates a Perpendicular Vector to the Vector Provided

        ## Parameters : \n
            self : Molecule - Default Parameter for the Class Instance \n
            vector : np.ndarray - Numpy Array of the Vector to Generate the Perpendicular Vector from

        ## Returns : \n
            np.ndarray - Vector that is Perpendicular to the one inputted
        """
        # Check if the input vector is not zero
        if np.linalg.norm(vector) == 0:
            raise ValueError("The input vector cannot be the zero vector")

        # Create a second vector which is not parallel to the input vector
        if vector[0] == 0:
            w = np.array(
                [1, 0, 0]
            )  # Choose a simple vector if the first element is zero
        else:
            w = np.array(
                [vector[1], -vector[0], 0]
            )  # This ensures the vector is not parallel

        # Compute the cross product
        crossProduct = np.cross(vector, w)

        # Normalize the cross product
        norm = np.linalg.norm(crossProduct)
        if norm == 0:
            raise ValueError(
                "The cross product resulted in a zero vector, which should not happen"
            )

        normalizedVector = crossProduct / norm

        return normalizedVector

    def rotateBond(self, atomIndex1: int, atomIndex2: int, radians: float):
        """Rotates a Bond and all Molecules proceeding it by the specified number of Radians

        ## Parameters : \n
            self : Molecule - Default Parameter for the Class Instance \n
            atomIndex1 : int - Index of the First Atom in the Bond \n
            atomIndex2 : int - Index of the Second Atom in the Bond \n
            radians : float - Radians to Rotate the Bond and all Atoms proceeding it

        ## Returns : \n
            None - No Return Value
        """

        atomIndexes: list[int] = self.getAllAtomsAfterBond(atomIndex1, atomIndex2)
        atomPositions: np.ndarray = np.array(
            self.getAtomPosition(atomIndexes[0]) - self.getAtomPosition(atomIndex2)
        )

        for i in range(1, len(atomIndexes)):
            atomPositions = np.vstack(
                (
                    atomPositions,
                    self.getAtomPosition(atomIndexes[i])
                    - self.getAtomPosition(atomIndex2),
                )
            )

        zVector = self.getAtomPosition(atomIndex2) - self.getAtomPosition(atomIndex1)
        zVector = zVector / np.linalg.norm(zVector)

        xVector = self.generatePerpendicularVector(zVector)

        yVector = np.cross(zVector, xVector)
        yVector = yVector / np.linalg.norm(yVector)

        originalAxes = np.array(
            [
                xVector,  # Original X-axis
                yVector,  # Original Y-axis
                zVector,  # Original Z-axis
            ]
        )

        # Define rotation matrix around Z-axis by given radians
        cosTheta = np.cos(radians)
        sinTheta = np.sin(radians)
        zRotationMatrix = np.array(
            [[cosTheta, -sinTheta, 0], [sinTheta, cosTheta, 0], [0, 0, 1]]
        )

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

    def getConformers(self, atomIndex1: int, atomIndex2: int, steps: int):
        """Generates a List of Conformer Molecules with each Conformers specified bond rotated by (2pi / steps) radians

        ## Parameters : \n
            self : Molecule - Default Parameter for the Class Instance \n
            atomIndex1 : int - Index of the First Atom in the Bond \n
            atomIndex2 : int - Index of the Second Atom in the Bond \n
            steps : int - Number of Steps / Conformer Molecules to be outputted

        ## Returns : \n
            list[Molecule] - List of Conformer Molecules along the Specified bond
        """
        stepSizeRad = (2 * pi) / steps
        conformers: list[Molecule] = []

        # Loops and creates a new Rotated Molecule to add to the List of Conformers
        for i in range(steps):
            rotation = stepSizeRad * i
            newMolecule = self.copy()
            newMolecule.name = f"{self.name}_rot_{(180 / pi) * rotation}"
            newMolecule.rotateBond(atomIndex1, atomIndex2, rotation)
            conformers.append(newMolecule)

        return conformers

    def XYZ(self):
        """Returns the XYZ Files content as a string generated from the Molecules

        ## Parameters : \n
            self : Molecule - Default Parameter for the Class Instance

        ## Returns : \n
            str - Content of XYZ File generated from the Molecule
        """
        from qchem.XYZFile import XYZFile

        return XYZFile(self).getFileAsString()

    def XYZBody(self) -> str:
        """Returns the Body of the XYZ File as a string generated from the Molecules

        ## Parameters : \n
            self : Molecule - Default Parameter for the Class Instance

        ## Returns : \n
            str - Body of XYZ File generated from the Molecule
        """
        from qchem.XYZFile import XYZFile

        return XYZFile(self).getXYZBody()

    def saveAsXYZ(self, fileDir: str):
        """Saves the Molecule as a XYZ File to the specified file directory. Uses the Molcules name as the File Name and the file directory as the folders path

        ## Parameters : \n
            self : Molecule - Default Parameter for the Class Instance \n
            fileDir : str - Path of the Directory the XYZ file is saved (Does not include name or file extension)

        ## Returns : \n
            None - No Return Value
        """
        from qchem.XYZFile import XYZFile

        XYZFile(self).saveToFile(fileDir)

    def findRotatableBonds(self):
        """Finds all rotatable bonds in the Molecule and adds it to the Bonds DataFrame under the column “Rotatable”. Rotatable bonds are not completely accurate, relies on finding a loop in the molecule (Double and Tripple bonds are still considered rotatable in this case so beware)

        ## Parameters : \n
            self : Molecule - Default Parameter for the Class Instance

        ## Returns : \n
            None - No Return Value
        """
        # TODO: Will need to Factor in Double Bonds in the Future

        rotatableBonds: list[list[bool]] = []

        # Go through all Bonds in the Molecule. Do a Chain search
        for i in range(self.atomCount):

            # Cache the Bonds for the current Atom and initialize
            bonds = self.bonds["Bonds"][i]
            rotatableList: list[bool] = []

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
        """Generates the Bonds DataFrame for the Molecule, fills the DataFrame with Atomic indices in the Molecule, their atomic symbol, indices of the atoms they are bonded to and the distances of those bonds

        ## Parameters : \n
            self : Molecule - Default Parameter for the Class Instance

        ## Returns : \n
            None - No Return Value
        """
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

    def getAtomPosition(self, atomIndex: int):
        """Returns a Numpy Array with the XYZ Position of the Specified Atom

        ## Parameters : \n
            self : Molecule - Default Parameter for the Class Instance \n
            atomIndex : int - Index of the Atom in the Molecule

        ## Returns : \n
            NDArray[Any] - Numpy Array of the Atoms XYZ Position
        """
        return np.array(self.XYZCoordinates.iloc[atomIndex, slice(1, 4)], dtype=float)

    def getDihedralAngle(
        self, atomIndex1: int, atomIndex2: int, atomIndex3: int, atomIndex4: int
    ) -> float:
        """Gets the Dihedral Angle in Degrees between 2 Bonds in the Molecule. (Gets the Angle between 2 vectors represented by bonds)

        ## Parameters : \n
            self : Molecule - Default Parameter for the Class Instance \n
            atomIndex1 : int - Index of the First Atom \n
            atomIndex2 : int - Index of the Second Atom \n
            atomIndex3 : int - Index of the Third Atom \n
            atomIndex4 : int - Index of the Fourth Atom

        ## Returns : \n
            float - Angle between 2 Bonds in Degrees
        """
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

    def getAngleBetweenAtoms(self, atomIndex1: int, atomIndex2: int, atomIndex3: int):
        """Determines the Angle in Degrees between 2 atoms bonded with a common middle atom

        ## Parameters : \n
            self : Molecule - Default Parameter for the Class Instance \n
            atomIndex1 : int - Index of the First Atom \n
            atomIndex2 : int - Index of the Second Atom, is the commonly shared Atom \n
            atomIndex3 : int - Index of the Third Atom \n

        ## Returns : \n
            float - Angle between the bonded atoms in Degrees
        """
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
        """Prints the Bond graph to the Terminal. Displays the index of atoms that are bonded to the displayed index

        ## Parameters : \n
            self : Molecule - Default Parameter for the Class Instance \n

        ## Returns : \n
            None - No Return Value
        """
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
            print(
                " %4i   %-2s - %s    %4s  %2s"
                % (index + 1, atom, bonds, bondDist, rotatable)
            )

    def getDihedralAtomChain(self, atomChain: list[int], depth: int = 0):
        """Recursively searches for a atom chain of length 4 to calculate the Dihedral angle

        ## Parameters : \n
            self : Molecule - Default Parameter for the Class Instance \n
            atomChain : list[int] - List of Atoms indices that are in the Chain \n
            depth : int - The Depth of the Recursive Search

        ## Returns : \n
            list[int] - The Dihedral Atom chain list. Is of length 4
        """
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
        """Displays the Molecules Atoms info in Z Matrix Format

        ## Parameters : \n
            self : Molecule - Default Parameter for the Class Instance

        ## Returns : \n
            None - No Return Value
        """
        zMatrix = self.createZMatrixFile()
        for i in zMatrix:
            print(i)

    def createZMatrixFile(self):
        """Converts the Molecules info from XYZ Format to Z Matrix Format as a list of strings

        ## Parameters : \n
            self : Molecule - Default Parameter for the Class Instance

        ## Returns : \n
            list[str] - List of Lines in the Z Matrix File 
        """
        # Create an Array Reprensenting the File, First line is Number of Atoms, Second is the Name of the Molecule
        zMatrix : list[str] = []
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

    def sortAtomDataFrame(self, dataFrame: pd.core.frame.DataFrame):
        """Sorts the Molecules Atoms Data Frame by ordering the Atoms from Heaviest Molecular Weight to Lowest

        ## Parameters : \n
            self : Molecule - Default Parameter for the Class Instance \n
            dataFrame : pd.DataFrame - The XYZ Position DataFrame to Sort

        ## Returns : \n
            pd.DataFrame - The DataFrame with Atom Positions sorted from Highest to Lowest Molecular Weight
        """
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

    def getMolecularWeight(self) -> float:
        """Returns the Molecular Weight of the Molecule in g/mol (grams per mol)

        ## Parameters : \n
            self : Molecule - Default Parameter for the Class Instance

        ## Returns : \n
            float - The Molecular Weight of the Molecule in g/mol
        """
        mw = 0

        # Loop through all Atoms and Add their Individual Atomic Mass
        for i in range(self.atomCount):
            mw += AtomicMassConstants[self.XYZCoordinates["Atom"][i]]

        return mw
