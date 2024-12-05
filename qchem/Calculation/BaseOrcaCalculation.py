import os
from .BaseCalculation import BaseCalculation
from ..Molecule import Molecule
from abc import abstractmethod

class BaseOrcaCalculation(BaseCalculation):
    
    molecule: Molecule | str
    """The Molecule Used for Calculation"""
    
    calculationType: str
    """The Calculation Type"""
    
    inputFilePath: str
    """The Path to the Input File on the Device"""

    outputFilePath: str
    """The Path to the Output File on the Device"""
    
    def __init__ (self, name: str, molecule: str | Molecule, index: int = 1, cores: int = 1, isLocal: bool = False, stdout : bool = True):
        
        # Make a Super Call (Use the Base Class Init for some Boilerplate Setup)
        super().__init__(name, index, cores, isLocal, stdout)
        
        # Check if Values are empty or of Wrong Type
        if (not (molecule and isinstance(molecule, (str, Molecule)))):
            raise ValueError("Molecule is not defined! Provide a Path to the XYZ file or a Molecule Object")
    
        if (not isinstance(cores, (int))):
            raise ValueError("Cores must be an integer")
    
        # Set the Values
        self.molecule = molecule
        self.SetName(name)
        
        # Convert to Functions?
        self.outputFilePath = os.path.join(self.orcaCachePath, self.GetOutputFileName())  #f'{self.OrcaCachePath}\\{self.GetOutputFileName()}'
        self.inputFilePath = os.path.join(self.orcaCachePath, self.GetInputFileName()) #f'{self.OrcaCachePath}\\{self.GetInputFileName()}'
    
    
    @abstractmethod
    def RunCalculation(self):
        """Runs the Orca Calculation"""
        pass
    
    def IsFileReference(self):
        """Determines if the Molecule is stored as a File Reference, or is a direct Molecule Object"""
        if (isinstance(self.molecule, (str))):
            return True
        else:
            return False
        
    def GetInputFileName (self):
        """Returns the Input File Name with it's extension"""
        return f"{self.calculationName}.inp"

    def GetOutputFileName (self):
        """Returns the Output File Name with it's extension"""
        return f"{self.calculationName}.out"
    
    def GetOutput (self) -> str:
        # Open the Output File and Grab the Content
        with open(self.outputFilePath, 'r') as file:
            self.CalculationOutput = file.read()
        
    def SetName(self, name:str):
        if (name == ""):
            if (isinstance(self.molecule, (Molecule))):
                self.name = self.molecule.Name
            else:
                print(f"No Name Provided for the Molecule, using Default Name: {self.calculationType}Molecule")
                self.name = f"{self.calculationType}Molecule"
        else:
            self.name = name
    