import os
from abc import ABC, abstractmethod

class BaseCalculation(ABC):
    """Abstract Class to be Inherited to define a new Orca Calculation class"""
    
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
    
    def __init__ (self, name: str, index: int = 1, cores: int = 1, isLocal: bool = False, stdout : bool = True):
        
        # Value Type Checking
        if (not isinstance(name, (str)) or name == ""):
            raise ValueError("Name of the Calculation must be specified ")
        
        if (not isinstance(index, (int))):
            raise ValueError("Index must be an integer")
        
        if (not isinstance(cores, (int))):
            raise ValueError("Cores must be an integer")
        
        if (not isinstance(isLocal, (bool))):
            raise ValueError("IsLocal must be a boolean")
        
        if (not isinstance(stdout, (bool))):
            raise ValueError("STDOut must be a boolean")
        
        #if (not isinstance(inputFile, (OrcaInputFile))):
        #    raise ValueError("Input File must be a OrcaInputFile")
        
        # Set Values
        self.name = name
        self.index = index
        self.cores = cores
        self.STDOut = stdout
        self.isLocal = isLocal
        
        #self.inputFile = inputFile
        
        # Generate Cache Paths
        orcaCache = "OrcaCache" #Change this to QChemCache?
        self.orcaCachePath =  os.path.join(os.getcwd(), orcaCache, self.name)  #f'{os.getcwd()}\\{orcaCache}\\{self.CalculationName}'
        #self.outputFilePath = os.path.join(self.orcaCachePath, self.GetOutputFileName())  #f'{self.OrcaCachePath}\\{self.GetOutputFileName()}'
        #self.inputFilePath = os.path.join(self.orcaCachePath, self.GetInputFileName()) #f'{self.OrcaCachePath}\\{self.GetInputFileName()}'
        
        # Determine how many cores are used by the calculation
        #if ("cores" in self.InputFile.variables):
        #    self.cores = self.InputFile.variables["cores"]
        #else:
        #    self.cores = 1
    
    @abstractmethod
    def RunCalculation(self):
        """Runs the Orca Calculation"""
        pass
    
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