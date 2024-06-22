import os
import pandas as pd
from qchem.Molecule import Molecule

class XYZFile:
    """Class that Describes a XYZ Molecule File"""

    AtomCount: int
    """Number of Atoms in the Molecule"""

    MoleculeName: str
    """Name of the Molecule"""

    AtomPositions: list[str]
    """Positions of each Atom in the Molecule in String Arrays"""

    def __init__(self, molecule : Molecule):
        self.AtomCount = molecule.AtomCount
        self.MoleculeName = molecule.Name
        self.AtomPositions = []

        positionDF =  molecule.XYZCoordinates

        for i in range(self.AtomCount):
            self.AtomPositions.append(f"{positionDF["Atom"][i]} {positionDF["X"][i]} {positionDF["Y"][i]} {positionDF["Z"][i]}")

    def GetFileAsString(self):
        """Returns the Entire XYZ File as a String"""
        fileString = f"{self.AtomCount}"
        fileString += f"\n{self.MoleculeName}"

        for atom in self.AtomPositions:
            fileString += f"\n{atom}"

        return fileString
    
    def SaveToFile (self, filePath: str):
        """Saves the XYZ File to the Specified Path"""
        with open(os.path.join(filePath, f"{self.MoleculeName}.xyz"), "w") as file:
            file.write(self.GetFileAsString())

    def GetXYZBody (self):
        """Gets the XYZ Atom Position Body of the file as a String"""
        string = self.AtomPositions[0]
        
        for i in range(2, self.AtomCount):
            string += f"\n{self.AtomPositions[i]}"

        return string