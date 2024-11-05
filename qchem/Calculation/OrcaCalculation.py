import os
import subprocess
import time
from qchem.Data.Enums import OrcaInputTemplate
from qchem.XYZFile import XYZFile
from qchem.Molecule import Molecule
from .OrcaInputFile import OrcaInputFile

class OrcaCalculation:
    """Class capable of running an Orca Calculation"""

    Cores : int
    """Number of Cores to use for the Calculation"""

    CalculationName: str
    """Name of the Calculation. Saves the Input and Output File as the Name"""

    OutputFilePath: str
    """The Path to the Output File on the Device"""

    InputFilePath: str
    """The Path to the Input File on the Device"""
    
    OrcaCachePath: str
    """The Path to the Orca Cache on the Local Device"""

    InputFile: OrcaInputFile
    """The Input File for the Calculation that will be run"""

    Index: int
    """Container Index so that they can be run in parallel"""
    
    IsLocal: bool
    """Boolean Flag determining if the Calculation should be run locally or inside a Container"""
    
    CalculationTime: float
    """The Time it took to run the Calculation (In Seconds)"""
    
    STDOut: bool
    """Boolean Flag Determining if Standard Output Messages will be displayed"""

    def __init__ (self, name: str, inputFile: OrcaInputFile, index: int = 1, isLocal: bool = False, stdout : bool = True):
        
        # Value Type Checking
        if (not isinstance(name, (str)) or name == ""):
            raise ValueError("Name of the Calculation must be specified ")
        
        if (not isinstance(inputFile, (OrcaInputFile))):
            raise ValueError("Input File must be a OrcaInputFile")
        
        if (not isinstance(index, (int))):
            raise ValueError("Index must be an integer")
        
        if (not isinstance(isLocal, (bool))):
            raise ValueError("IsLocal must be a boolean")
        
        if (not isinstance(stdout, (bool))):
            raise ValueError("STDOut must be a boolean")
        
        # Set Values
        self.CalculationName = name
        self.InputFile = inputFile
        self.Index = index
        self.IsLocal = isLocal
        self.STDOut = stdout
        
        # Generate Cache Paths
        orcaCache = "OrcaCache"
        self.OrcaCachePath =  os.path.join(os.getcwd(), orcaCache, self.CalculationName)  #f'{os.getcwd()}\\{orcaCache}\\{self.CalculationName}'
        self.OutputFilePath = os.path.join(self.OrcaCachePath, self.GetOutputFileName())  #f'{self.OrcaCachePath}\\{self.GetOutputFileName()}'
        self.InputFilePath = os.path.join(self.OrcaCachePath, self.GetInputFileName()) #f'{self.OrcaCachePath}\\{self.GetInputFileName()}'
        
        # Determine how many cores are used by the calculation
        if ("cores" in self.InputFile.variables):
            self.Cores = self.InputFile.variables["cores"]
        else:
            self.Cores = 1

    def RunCalculation(self):
        """Runs a Orca Calculation in a Docker Container """

        # Create the Directories to Store Calculation info
        self.CreateDirectories()

        # Save the Input File to the folder
        self.InputFile.SaveInputFile(os.path.join(self.OrcaCachePath, self.CalculationName) + ".inp")

        # Get the Start Time of the Calculation
        startTimer = time.time()
        
        if self.STDOut:
            print(f"Running Calulation : {self.GetInputFileName()}")

        # Determine of running the 
        if (self.IsLocal):
            result = self.RunLocally(self.OrcaCachePath)
        else:
            result = self.RunDockerContainer(self.OrcaCachePath)
        
        # Determine the Elapsed Time
        self.CalculationTime = time.time() - startTimer
        
        # Post a message that an Error may have Occured
        if (result.stderr.__len__() > 0):
            print(f"WARNING Errors Maybe Occured : \n\n{result.stderr}")
            
        if self.STDOut:
            print(f"Calculation Complete ({self.ClockTime(self.CalculationTime)}) : {self.GetInputFileName()}")  

    def RunLocally (self, cachePath: str):
        # Create the Command String
        command = ""
        
        # Windows OS
        if os.name == 'nt':
            orcaPath = subprocess.run("where orca", shell=True, text=True, capture_output=True).stdout.strip(" \n\"").removesuffix(".exe")
            command = f"cd /d {cachePath} && \"{orcaPath}\" \"{self.GetInputFileName()}\" > \"{self.GetOutputFileName()}\""
        else:
            # Unix based OS (Linux, Mac)
            command = f"cd \"{cachePath}\" && /Orca/orca {self.GetInputFileName()} > {self.GetOutputFileName()}"
        
        # Run the Orca Calculation locally
        return subprocess.run(command, shell=True, text=True, capture_output=True)
            
    def RunDockerContainer (self, cachePath):
        # Create the Command String
        command = f'docker run --name qchemorca{self.Index} -v "{cachePath}":/home/orca mrdnalex/orca sh -c "cd /home/orca && /Orca/orca {self.GetInputFileName()} > {self.GetOutputFileName()}"'

        # Kill and Remove qchemorca container if it doesn't exist yet
        subprocess.run(f"docker kill qchemorca{self.Index}" , shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
        subprocess.run(f"docker rm qchemorca{self.Index}" , shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
        
        # Run the Calculation in a Container and wait
        result = subprocess.run(command, shell=True, text=True, capture_output=True)
        
        # Kill and Remove the Container
        subprocess.run(f"docker kill qchemorca{self.Index}" , shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
        subprocess.run(f"docker rm qchemorca{self.Index}" , shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
        
        return result

    def CreateDirectories (self):
        """Creates the necessary Cache Folders to store the Calculation"""
        # Make Cache Folder if it doesn't Exist
        if not os.path.exists(self.OrcaCachePath):
            os.makedirs(self.OrcaCachePath)

        # Make a folder for the Specific Calculation
        if not os.path.exists(self.OrcaCachePath):
            os.makedirs(self.OrcaCachePath)

    def GetInputFileName (self):
        """Returns the Input File Name with it's extension"""
        return f"{self.CalculationName}.inp"

    def GetOutputFileName (self):
        """Returns the Output File Name with it's extension"""
        return f"{self.CalculationName}.out"
    
    def GetOutput (self) -> str:
        # Open the Output File and Grab the Content
        with open(self.OutputFilePath, 'r') as file:
            self.CalculationOutput = file.read()
            
    def ClockTime(self, seconds):
        hours = seconds // 3600
        minutes = (seconds % 3600) // 60
        remaining_seconds = seconds % 60
        
        parts = []
        if hours > 0:
            parts.append(f"{int(hours)} hour{'s' if hours > 1 else ''}")
        if minutes > 0:
            parts.append(f"{int(minutes)} minute{'s' if minutes > 1 else ''}")
        if remaining_seconds > 0:
            parts.append(f"{int(remaining_seconds)} second{'s' if remaining_seconds > 1 else ''}")
        
        return ", ".join(parts) if parts else "0 seconds"


