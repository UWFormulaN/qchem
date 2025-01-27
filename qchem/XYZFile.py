import os
import re
import pandas as pd
from .Data.Constants import AtomicMassConstants


class XYZFile:
    """Describes the Structure of a XYZ file. Allows for loading and saving to the file format"""

    atomCount: int
    """Number of Atoms Present in the File"""

    moleculeName: str
    """Name of the Molecule being stored or loaded in the file"""

    atomPositions: pd.DataFrame
    """Pandas DataFrame storing the Position of each Atom in the Molecule, stored in the Rows"""

    def __init__(self, molecule: str | list[str]):

        # Check if it's a String path
        if isinstance(molecule, list) and all(isinstance(item, str) for item in molecule):
            self.moleculeName = molecule[1].strip()
            
            atoms = []
            columnHeaders = ["Atom", "X", "Y", "Z"]
            
            for i in range(len(molecule)):
                if (self.isValidXYZLine(molecule[i].strip())):
                    atoms.append(re.sub(r'\s+', " ", molecule[i].strip()).split(" "))
                
            self.atomCount = len(atoms)
            self.atomPositions = pd.DataFrame(atoms, columns=columnHeaders)
        
        elif isinstance(molecule, (str)) and os.path.exists(molecule):
            # Open file and Extract all Lines
            xyzFileLines = open(molecule).readlines()

            # Check if first line is a Integer, likely to be Atom Count
            if xyzFileLines[0].strip().isdigit():
                self.atomCount = int(xyzFileLines[0].strip())

            # Set Second line as Molecule Name
            self.moleculeName = xyzFileLines[1].strip()
            
            # Load the Positions of the Atoms
            self.atomPositions = self.sortAtomDataFrame(self.readXYZ(molecule))
        
        else:
            raise ValueError("Invalid Input, must be a Path to a File or a List of Strings")

    def isValidXYZLine(self, line: str):
        """Checks if the string provided follows the format of the Body of a XYZ File (Atom Symbol - X Y Z)

        ## Parameters : \n
            self : XYZFile - Default Parameter for the Class Instance \n
            line : str - The Line to be Checked

        ## Returns : \n
            bool - True if the line follows (Atom Symbol - X Y Z), False Otherwise
        """
        splitLine = re.sub(r"\s+", " ", line).split(" ")

        if len(splitLine) == 4:
            if (
                self.isValidFloat(splitLine[1])
                and self.isValidFloat(splitLine[2])
                and self.isValidFloat(splitLine[3])
            ):
                return True

        return False

    def isValidFloat(self, value: str):
        """Checks if the String provided can be converted to a float

        ## Parameters : \n
            self : XYZFile - Default Parameter for the Class Instance \n
            value : str - The String that will be checked if it can be converted to a float

        ## Returns : \n
            bool - True if the String can be converted to a Float, False Otherwise
        """
        try:
            float(value)
            return True
        except ValueError:
            return False

    def getFileAsString(self):
        """Formats the Classes content into a XYZFile format in a single String

        ## Parameters : \n
            self : XYZFile - Default Parameter for the Class Instance

        ## Returns : \n
            str - Content of the XYZFile formatted properly in a single String
        """
        fileString = f"{self.atomCount}"
        fileString += f"\n{self.moleculeName}\n"
        fileString += self.atomPositions.to_string(header=False, index=False)
        return fileString

    def saveToFile(self, directory: str = ""):
        """Saves the XYZFile Class Object to a XYZFile

        ## Parameters : \n
            self : XYZFile - Default Parameter for the Class Instance \n
            directory : str - The Directory the XYZ File is saved to (Do not include Name and File Extension)

        ## Returns : \n
            None - No Return Value
        """
        with open(os.path.join(directory, f"{self.moleculeName}.xyz"), "w") as file:
            file.write(self.getFileAsString())

    def getXYZBody(self):
        """Retrieves the Body of the XYZ File, The Atomic Symbols and their XYZ Positions. Excludes the Name of the Molecule and the number of Atoms

        ## Parameters : \n
            self : XYZFile - Default Parameter for the Class Instance \n

        ## Returns : \n
            str - The Body of the XYZ File as a String (Rows of A - X Y Z)
        """
        return self.atomPositions.to_string(header=False, index=False)

    def sortAtomDataFrame(self, dataFrame: pd.DataFrame):
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

    def readXYZ(self, path: str) -> pd.DataFrame:
        """Reads the Provided XYZ Files and returns the XYZ Format in a DataFrame

        ## Parameters : \n
            self : Molecule - Default Parameter for the Class Instance
            path : str - Path to the XYZ File

        ## Returns : \n
            pd.DataFrame - Pandas DataFrame with each row being an Atoms XYZ Position
        """
        if not os.path.exists(path):
            raise FileNotFoundError(f"File {path} not found")
    
        return pd.read_csv(
            path, sep=r"\s+", skiprows=2, names=["Atom", "X", "Y", "Z"], engine="python"
        )
