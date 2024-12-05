import os
import subprocess
import time
from .OrcaInputFile import OrcaInputFile
from .BaseCalculation import BaseCalculation
from ..Molecule import Molecule #from qchem.Molecule import Molecule

class OrcaCalculation(BaseCalculation):
    """Class capable of running an Orca Calculation"""

    #Cores : int
    #"""Number of Cores to use for the Calculation"""
#
    #CalculationName: str
    #"""Name of the Calculation. Saves the Input and Output File as the Name"""
#
    #OutputFilePath: str
    #"""The Path to the Output File on the Device"""
#
    #InputFilePath: str
    #"""The Path to the Input File on the Device"""
    #
    #OrcaCachePath: str
    #"""The Path to the Orca Cache on the Local Device"""
#
    #InputFile: OrcaInputFile
    #"""The Input File for the Calculation that will be run"""
#
    #Index: int
    #"""Container Index so that they can be run in parallel"""
    #
    #IsLocal: bool
    #"""Boolean Flag determining if the Calculation should be run locally or inside a Container"""
    #
    #CalculationTime: float
    #"""The Time it took to run the Calculation (In Seconds)"""
    #
    #STDOut: bool
    #"""Boolean Flag Determining if Standard Output Messages will be displayed"""

    inputFilePath: str
    """The Path to the Input File on the Device"""

    outputFilePath: str
    """The Path to the Output File on the Device"""

    isLocal: bool
    """Boolean Flag determining if the Calculation should be run locally or inside a Container"""
    
    inputFile: OrcaInputFile
    """The Input File for the Calculation that will be run"""
    
    molecule: Molecule | str
    """The Molecule to be Optimized"""

    def __init__ (self, name: str, inputFile: OrcaInputFile, index: int = 1, isLocal: bool = False, stdout : bool = True):
        
        # Make a Super Call (Use the Base Class Init for some Boilerplate Setup)
        super().__init__(name, index, stdout)
        
        # Value Type Checking
        if (not isinstance(inputFile, (OrcaInputFile))):
            raise ValueError("Input File must be a OrcaInputFile")
        
        if (not isinstance(isLocal, (bool))):
            raise ValueError("IsLocal must be a boolean")
        
        #if (not isinstance(name, (str)) or name == ""):
        #    raise ValueError("Name of the Calculation must be specified ")
        #
        #if (not isinstance(index, (int))):
        #    raise ValueError("Index must be an integer")
        #
        #if (not isinstance(isLocal, (bool))):
        #    raise ValueError("IsLocal must be a boolean")
        #
        #if (not isinstance(stdout, (bool))):
        #    raise ValueError("STDOut must be a boolean")
        
        # Set Values
        self.inputFile = inputFile
        self.isLocal = isLocal
        
        #self.calculationName = name
        #self.index = index
        #self.isLocal = isLocal
        #self.STDOut = stdout
        
        # Generate Cache Paths
        #orcaCache = "OrcaCache"
        #self.orcaCachePath =  os.path.join(os.getcwd(), orcaCache, self.calculationName)  #f'{os.getcwd()}\\{orcaCache}\\{self.CalculationName}'
        self.outputFilePath = os.path.join(self.orcaCachePath, self.GetOutputFileName())  #f'{self.OrcaCachePath}\\{self.GetOutputFileName()}'
        self.inputFilePath = os.path.join(self.orcaCachePath, self.GetInputFileName()) #f'{self.OrcaCachePath}\\{self.GetInputFileName()}'
        
        # Determine how many cores are used by the calculation
        if ("cores" in self.inputFile.variables):
            self.cores = self.inputFile.variables["cores"]
        else:
            self.cores = 1

    def RunCalculation(self):
        """Runs a Orca Calculation in a Docker Container"""

        # Create the Directories to Store Calculation info
        self.CreateDirectories()

        # Save the Input File to the folder
        self.inputFile.SaveInputFile(os.path.join(self.orcaCachePath, self.calculationName) + ".inp")

        # Get the Start Time of the Calculation
        startTimer = time.time()
        
        if self.STDOut:
            print(f"Running Calulation : {self.GetInputFileName()}")

        # Run the Calculation Locally or through a Docker Container
        if (self.isLocal):
            result = self.RunLocally(self.orcaCachePath)
        else:
            result = self.RunDockerContainer(self.orcaCachePath)
        
        # Get the Total Calculation time
        self.calculationTime = time.time() - startTimer
        
        # Post a message that an Error may have Occured
        if (result.stderr.__len__() > 0):
            print(f"WARNING Errors Maybe Occured : \n\n{result.stderr}")
        
        # If Standard Output Allowed post the Completion Message
        if self.STDOut:
            print(f"Calculation Complete ({self.ClockTime(self.calculationTime)}) : {self.GetInputFileName()}")  

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
        command = f'docker run --name qchemorca{self.index} -v "{cachePath}":/home/orca mrdnalex/orca sh -c "cd /home/orca && /Orca/orca {self.GetInputFileName()} > {self.GetOutputFileName()}"'

        # Kill and Remove qchemorca container if it doesn't exist yet
        subprocess.run(f"docker kill qchemorca{self.index}" , shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
        subprocess.run(f"docker rm qchemorca{self.index}" , shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
        
        # Run the Calculation in a Container and wait
        result = subprocess.run(command, shell=True, text=True, capture_output=True)
        
        # Kill and Remove the Container
        subprocess.run(f"docker kill qchemorca{self.index}" , shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
        subprocess.run(f"docker rm qchemorca{self.index}" , shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
        
        return result

    #def CreateDirectories (self):
    #    """Creates the necessary Cache Folders to store the Calculation"""
    #    # Make Cache Folder if it doesn't Exist
    #    if not os.path.exists(self.orcaCachePath):
    #        os.makedirs(self.orcaCachePath)
#
    #    # Make a folder for the Specific Calculation
    #    if not os.path.exists(self.orcaCachePath):
    #        os.makedirs(self.orcaCachePath)

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
    
    
    #def ClockTime(self, seconds):
    #    hours = seconds // 3600
    #    minutes = (seconds % 3600) // 60
    #    remaining_seconds = seconds % 60
    #    
    #    parts = []
    #    if hours > 0:
    #        parts.append(f"{int(hours)} hour{'s' if hours > 1 else ''}")
    #    if minutes > 0:
    #        parts.append(f"{int(minutes)} minute{'s' if minutes > 1 else ''}")
    #    if remaining_seconds > 0:
    #        parts.append(f"{int(remaining_seconds)} second{'s' if remaining_seconds > 1 else ''}")
    #    
    #    return ", ".join(parts) if parts else "0 seconds"

    def IsFileReference(self):
        """Determines if the Molecule is stored as a File Reference, or is a direct Molecule Object"""
        if isinstance(self.molecule, (str)):
            return True
        else:
            return False

