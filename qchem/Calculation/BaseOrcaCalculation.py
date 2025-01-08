import os
from typing import Any
from qchem.Data.Enums import OrcaInputTemplate
from .OrcaInputFile import OrcaInputFile
from ..Molecule import Molecule
from abc import ABC, abstractmethod

class BaseOrcaCalculation(ABC):
    
    molecule: Molecule | str
    """The Molecule Used for Calculation"""
    
    template: str | OrcaInputTemplate
    """Input File Tenmplate for the Calculation"""
    
    calculationType: str = ""
    """The Orca Calculation Type and Keyword"""
    
    inputFile: OrcaInputFile
    """Input File Object"""
    
    inputFilePath: str
    """The Path to the Input File on the Device"""

    outputFilePath: str
    """The Path to the Output File on the Device"""
    
    variables: dict[str, Any]
    """Additional Variables to be Passed Through to the Input File"""
    
    index: int
    """Container Index so that they can be run in parallel"""

    calculationName: str
    """Name of the Calculation. Saves the Input and Output File as the Name"""

    orcaCachePath: str
    """The Path to the Orca Cache on the Local Device"""

    calculationTime: float
    """The Time it took to run the Calculation (In Seconds)"""
    
    STDOut: bool
    """Boolean Flag Determining if Standard Output Messages will be displayed"""
    
    isLocal: bool
    """Boolean Flag determining if the Calculation should be run locally or inside a Container"""
    
    cores : int
    """Number of Cores to use for the Calculation"""
    
    def __init__ (self, name: str, molecule: str | Molecule, template: str | OrcaInputTemplate, index: int, cores: int, isLocal: bool, stdout : bool, **variables):
        
         # Value Type Checking
        if (not (self.calculationType and isinstance(self.calculationType, (str)))):
            raise ValueError("Calculation Type must be defined and a String")
         
        if (not isinstance(name, (str)) or name == ""):
            raise ValueError("Name of the Calculation must be specified ")
        
        if (not (molecule and isinstance(molecule, (str, Molecule)))):
            raise ValueError("Molecule must be a Molecule Object or String")
        
        if (not isinstance(template, (str, OrcaInputTemplate))):
            raise ValueError("Template must be a String or OrcaInputTemplate")
        
        if (not isinstance(index, (int))):
            raise ValueError("Index must be an integer")
        
        if (not isinstance(cores, (int))):
            raise ValueError("Cores must be an integer")
        
        if (not isinstance(isLocal, (bool))):
            raise ValueError("IsLocal must be a boolean")
        
        if (not isinstance(stdout, (bool))):
            raise ValueError("STDOut must be a boolean")
        
        # Set the Values
        self.name = name
        self.molecule = molecule
        self.template = template
        self.index = index
        self.cores = cores
        self.isLocal = isLocal
        self.STDOut = stdout
        self.variables = variables
        
        # Set Name to the Name of the Molecule
        if (isinstance(molecule, Molecule)):
            self.name = molecule.Name
        
        # Set Appropriate XYZ Format and Default Templates
        if (self.IsFileReference()):
            self.variables["xyzfile"] = self.molecule
            if (self.template == ""):
                if (self.cores == 1):
                    self.template = OrcaInputTemplate.BASIC
                else:
                    self.template = OrcaInputTemplate.BASICPARALLEL
        else:
            self.variables["xyz"] = self.molecule.XYZBody()
            if (self.template == ""):
                if (self.cores == 1):
                    self.template = OrcaInputTemplate.BASICXYZ
                else:
                    self.template = OrcaInputTemplate.BASICXYZPARALLEL
                
        self.variables["calculation"] = self.calculationType
        self.variables["cores"] = self.cores
            
        # Generate Cache Paths
        orcaCache = "OrcaCache" #Change this to QChemCache?
        self.orcaCachePath =  os.path.join(os.getcwd(), orcaCache, self.name)  #f'{os.getcwd()}\\{orcaCache}\\{self.CalculationName}'
        
        # Convert to Functions?
        self.outputFilePath = os.path.join(self.orcaCachePath, self.GetOutputFileName())  #f'{self.OrcaCachePath}\\{self.GetOutputFileName()}'
        self.inputFilePath = os.path.join(self.orcaCachePath, self.GetInputFileName()) #f'{self.OrcaCachePath}\\{self.GetInputFileName()}'

        # Create the Input File
        self.inputFile = OrcaInputFile(self.template, **self.variables)
    
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
        return f"{self.name}.inp"

    def GetOutputFileName (self):
        """Returns the Output File Name with it's extension"""
        return f"{self.name}.out"
    
    def GetOutput (self) -> str:
        # Open the Output File and Grab the Content
        with open(self.outputFilePath, 'r') as file:
            self.CalculationOutput = file.read()
            
    def CreateDirectories (self):
        """Creates the necessary Cache Folders to store the Calculation"""
        # Make Cache Folder if it doesn't Exist
        if not os.path.exists(self.orcaCachePath):
            os.makedirs(self.orcaCachePath)

        # Make a folder for the Specific Calculation
        if not os.path.exists(self.orcaCachePath):
            os.makedirs(self.orcaCachePath)
    
    def ClockTime(self, seconds):
        """Converts Seconds to a Human Readable Time String"""
        # Convert Seconds to Hours, Minutes, and Seconds
        days = seconds // 86400
        hours = (seconds % 86400) // 3600
        minutes = (seconds % 3600) // 60
        remainingSeconds = seconds % 60

        # Generate the Time String
        parts = []
        if days > 0:
            parts.append(f"{int(hours)} day{'s' if days > 1 else ''}")
        if hours > 0:
            parts.append(f"{int(hours)} hour{'s' if hours > 1 else ''}")
        if minutes > 0:
            parts.append(f"{int(minutes)} minute{'s' if minutes > 1 else ''}")
        if remainingSeconds > 0:
            parts.append(
                f"{int(remainingSeconds)} second{'s' if remainingSeconds > 1 else ''}"
            )

        # Return the Time String
        return ", ".join(parts) if parts else "0 seconds"
    
    def BasisSetFunctionalCompliant (self):
        if (not ("basis" in self.variables)):
            raise ValueError("BasisSet not defined! Provide Basis Set Name as a String")
        
        if (not ("functional" in self.variables)):
            raise ValueError("Functional not defined! Provide Functional Name as a String")
        
    def SetCalculationType (self, calculationType):
        """Sets the Calculation Type in the Variables Dictionary"""
        self.variables["calculation"] = calculationType