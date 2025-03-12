import os
import pandas as pd
import re
from qchem.Molecule import Molecule


class XYZFile:
    """Class that Describes a XYZ Molecule File"""

    AtomCount: int
    """Number of Atoms in the Molecule"""

    MoleculeName: str
    """Name of the Molecule"""

    AtomPositions : pd.DataFrame
    """Positions of each Atom and the Atom Atomic Element in a DataFrame"""

    def __init__(self, molecule : Molecule | str | list[str], name : str = "Molecule"):
        # Check if it's already a Molecule
        if (isinstance(molecule, (Molecule))):
            self.AtomCount = molecule.AtomCount
            self.MoleculeName = molecule.Name
            self.AtomPositions = molecule.XYZCoordinates
            
        # Extract if it is a list of Strings 
        elif isinstance(molecule, list) and all(isinstance(item, str) for item in molecule):
            self.MoleculeName = name
            
            atoms = []
            columnHeaders = ["Atom", "X", "Y", "Z"]
            
            for i in range(len(molecule)):
                if (self.isValidXYZLine(molecule[i].strip())):
                    atoms.append(re.sub(r'\s+', " ", molecule[i].strip()).split(" "))
                
            self.AtomCount = len(atoms)
            self.AtomPositions = pd.DataFrame(atoms, columns=columnHeaders)
        # Check if it's a String path
        elif isinstance(molecule, (str)):
            
            # Open file and Extract all Lines
            xyzFileLines = open(molecule).readlines()
            
            # Check if first line is a Integer, likely to be Atom Count
            if (xyzFileLines[0].strip().isdigit()):
                self.AtomCount = int(xyzFileLines[0].strip())
            
            # Set Second line as Molecule Name
            self.MoleculeName = xyzFileLines[1].strip()
            
            # Default Arrays
            atoms = []
            columnHeaders = ["Atom", "X", "Y", "Z"]
            
            for i in range(len(xyzFileLines)):
                if (self.isValidXYZLine(xyzFileLines[i].strip())):
                    atoms.append(re.sub(r'\s+', " ", xyzFileLines[i].strip()).split(" "))
                    
            self.AtomPositions = pd.DataFrame(atoms, columns=columnHeaders)
                 
    def isValidXYZLine (self, line:str):
        """Checks if a String / Line from a File is a Valid XYZ File Format Line"""
        splitLine = re.sub(r'\s+', " ", line).split(" ")
        
        if (len(splitLine) == 4):
            if (self.isValidFloat(splitLine[1]) and self.isValidFloat(splitLine[2]) and self.isValidFloat(splitLine[3])):
                return True
        
        return False
            
    def isValidFloat(self, value):
        """Checks if the String provided can be converted to a Float"""
        try:
            float(value)
            return True
        except ValueError:
            return False        
        
    def GetFileAsString(self):
        """Returns the Entire XYZ File as a String"""
        fileString = f"{self.AtomCount}"
        fileString += f"\n{self.MoleculeName}\n"
        fileString += self.AtomPositions.to_string(header=False, index=False)
        return fileString
    
    def SaveToFile (self, directory: str = ""):
        """Saves the XYZ File to the Specified Path"""
        with open(os.path.join(directory, f"{self.MoleculeName}.xyz"), "w") as file:
            file.write(self.GetFileAsString())

    def GetXYZBody (self):
        """Gets the XYZ Atom Position Body of the file as a String"""
        return self.AtomPositions.to_string(header=False, index=False)