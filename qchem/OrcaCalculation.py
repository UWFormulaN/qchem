import os
import subprocess
from qchem.XYZFile import XYZFile
from qchem.Molecule import Molecule
from .OrcaInputFile import OrcaInputFile

class OrcaCalculation:
    """Class capable of running an Orca Calculation"""

    #CalculationMolecule: Molecule
    #"""Molecule that will have an Orca Calculation run on it"""
#
    #CalculationType : str
    #"""The Type of Calculation that will occur on the Molecule"""
    #
    #BasisSet: str
    #"""Basis Set to use for the Calculation"""
#
    #DensityFunctional : str
    #"""The Density Functional to use for the Calculation"""
#
    Cores : int
    """Number of Cores to use for the Calculation"""

    CalculationName: str

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

    def __init__ (self, name: str, inputFile: OrcaInputFile, index: int = 1, isLocal: bool = False):
        
        # Value Type Checking
        if (not isinstance(name, (str)) or name == ""):
            raise ValueError("Name of the Calculation must be specified ")
        
        if (not isinstance(inputFile, (OrcaInputFile))):
            raise ValueError("Input File must be a OrcaInputFile")
        
        if (not isinstance(index, (int))):
            raise ValueError("Index must be an integer")
        
        if (not isinstance(isLocal, (bool))):
            raise ValueError("IsLocal must be a boolean")
        
        self.CalculationName = name
        self.InputFile = inputFile
        self.Index = index
        self.IsLocal = isLocal
        orcaCache = "OrcaCache"
        self.OrcaCachePath = f'{os.getcwd()}\\{orcaCache}\\{self.CalculationName}'
        self.OutputFilePath = f'{self.OrcaCachePath}\\{self.GetOutputFileName()}'
        self.InputFilePath = f'{self.OrcaCachePath}\\{self.GetInputFileName()}'
        
        if ("cores" in self.InputFile.variables):
            self.Cores = self.InputFile.variables["cores"]
        else:
            self.Cores = 1


    #def __init__(self, molecule: Molecule, calculationType: str = "", basisSet: str = "", densityFunctional: str = "", cores : int = 1, index: int = 1, isLocal: bool = False):
    #    
    #    
    #    if (not (molecule and isinstance(molecule, (Molecule)))):
    #        raise ValueError("Molecule is not defined or is Not of Type Molecule")
    #    
    #    if (not isinstance(cores, (int))):
    #        raise ValueError("Cores must be an Integer")
#
    #    if (calculationType and not isinstance(calculationType, (str))):
    #        if (isinstance(calculationType, (str))):
    #            self.CalculationType = calculationType
    #        else:
    #            raise ValueError("Calculation Type must be a String")
    #        
    #    if (basisSet):
    #        if (isinstance(basisSet, (str))):
    #            self.BasisSet = basisSet
    #        else:
    #            raise ValueError("Basis Set Type must be a String")
    #        
    #    if (densityFunctional):
    #        if (isinstance(densityFunctional, (str))):
    #            self.DensityFunctional = densityFunctional
    #        else:
    #            raise ValueError("Density Functional Type must be a String")
    #        
    #    if (not isinstance(index, (int))):
    #        raise ValueError("Index must be an integer")
    #    
    #    if (not isinstance(isLocal, (bool))):
    #        raise ValueError("IsLocal must be a boolean")
    #    
    #    
    #    self.Index = index
    #    self.Cores = cores
    #    self.CalculationMolecule = molecule
    #    self.IsLocal = isLocal
    #    orcaCache = "OrcaCache"
    #    self.OrcaCachePath = f'{os.getcwd()}\\{orcaCache}\\{self.CalculationMolecule.Name.replace('.', '')}'
    #    self.OutputFilePath = f'{self.OrcaCachePath}\\{self.GetOutputFileName()}'
    #    self.InputFilePath = f'{self.OrcaCachePath}\\{self.GetInputFileName()}'

    def RunCalculation(self):
        """Runs a Orca Calculation in a Docker Container """

        self.CreateDirectories()

        # Save the Input File to the folder
        self.InputFile.SaveInputFile(os.path.join(self.OrcaCachePath, self.CalculationName) + ".inp")

        if (self.IsLocal):
            result = self.RunLocally(self.OrcaCachePath)
        else:
            result = self.RunDockerContainer(self.OrcaCachePath)
            
        if (result.stderr.__len__() > 0):
            print(f"WARNING Errors Maybe Occured : \n\n{result.stderr}")

        print(f"Calculation Complete : {self.GetInputFileName()}")  

        # Open the Output File and Grab the Content
        with open(self.OutputFilePath, 'r') as file:
            self.CalculationOutput = file.read()
        
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

        print(f"Running Calulation : {self.GetInputFileName()}")

        # Kill and Remove qchemorca container if it doesn't exist yet
        subprocess.run(f"docker kill qchemorca{self.Index}" , shell=True)
        subprocess.run(f"docker rm qchemorca{self.Index}" , shell=True)
        
        # Run the Calculation in a Container and wait
        result = subprocess.run(command, shell=True, text=True, capture_output=True)

        # Kill and Remove the Container
        subprocess.run(f"docker kill qchemorca{self.Index}" , shell=True)
        subprocess.run(f"docker rm qchemorca{self.Index}" , shell=True)
        
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

    #def GetInputFileName (self):
    #    """Returns the Input File Name with it's extension"""
    #    return f"{self.CalculationMolecule.Name}.inp"
#
    #def GetOutputFileName (self):
    #    """Returns the Output File Name with it's extension"""
    #    return f"{self.CalculationMolecule.Name}.out"

    #def GetInputFile (self):
    #    """Generates a Input file in String format"""
#
    #    firstLine = "!"
#
    #    # Check if the Properties are defined and Add them to the First line of the Input File
    #    if (self.DensityFunctional):
    #        firstLine += f"{self.DensityFunctional}"
#
    #    if (self.BasisSet):
    #        firstLine += f" {self.BasisSet}"
#
    #    if (self.CalculationType):
    #        firstLine += f" {self.CalculationType}"
#
    #    if (self.Cores > 1):
    #        if (self.Cores < 9):
    #            firstLine += " " + f"PAL{self.Cores}"
    #        else:
    #            firstLine += f"\n%PAL NPROCS {self.Cores} END"
#
    #    # Generate XYZ Wrappers
    #    xyzWrapperStart = "* xyz 0 1"
    #    xyzWrapperEnd = "*"
#
    #    # Create a XYZ File for the Molecule
    #    xyz = XYZFile(self.CalculationMolecule)
#
    #    # Generate the Entire Input File as a String
    #    inputFile = f"{firstLine}\n{xyzWrapperStart}\n{xyz.GetXYZBody()}\n{xyzWrapperEnd}"
#
    #    return inputFile

    #def SaveInputFile (self, filePath: str):
    #    """Saves a Input File using the Settings Provided to the Path Specified"""
    #    
    #    
    #    
    #    with open(os.path.join(filePath, self.GetInputFileName()), "w") as file:
    #        file.write(self.InputFile.InputFile)