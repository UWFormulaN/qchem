import os
import subprocess
from textwrap import indent
import comm
from distutils import core
from qchem.XYZFile import XYZFile
from qchem.Molecule import Molecule

class OrcaCalculation:
    """Class capable of running an Orca Calculation"""

    CalculationMolecule: Molecule
    """Molecule that will have an Orca Calculation run on it"""

    CalculationType : str
    """The Type of Calculation that will occur on the Molecule"""
    
    BasisSet: str
    """Basis Set to use for the Calculation"""

    DensityFunctional : str
    """The Density Functional to use for the Calculation"""

    Cores : int
    """Number of Cores to use for the Calculation"""

    OutputFilePath: str
    """The Path to the Output File on the Device"""

    InputFilePath: str
    """The Path to the Input File on the Device"""
    
    OrcaCachePath: str
    """The Path to the Orca Cache on the Local Device"""

    Index: int
    """Container Index so that they can be run in parallel"""
    
    IsLocal: bool
    """Boolean Flag determining if the Calculation should be run locally or inside a Container"""

    def __init__(self, molecule: Molecule, calculationType: str = None, basisSet: str = None, densityFunctional: str = None, cores : int = 1, index: int = 1):
        
        if (molecule and isinstance(molecule, (Molecule))):
            self.CalculationMolecule = molecule
        else:
            raise ValueError("Molecule is not defined or is Not of Type Molecule")
        
        if (isinstance(cores, (int))):
            self.Cores = cores
        else:
            raise ValueError("Cores must be an Integer")

        if (calculationType):
            if (isinstance(calculationType, (str))):
                self.CalculationType = calculationType
            else:
                raise ValueError("Calculation Type must be a String")
            
        if (basisSet):
            if (isinstance(basisSet, (str))):
                self.BasisSet = basisSet
            else:
                raise ValueError("Basis Set Type must be a String")
            
        if (densityFunctional):
            if (isinstance(densityFunctional, (str))):
                self.DensityFunctional = densityFunctional
            else:
                raise ValueError("Density Functional Type must be a String")
            
        if (isinstance(index, (int))):
            self.Index = index
        else:
            raise ValueError("Index must be an integer")
        
        orcaCache = "OrcaCache"
        self.OrcaCachePath = f'{os.getcwd()}\\{orcaCache}\\{self.CalculationMolecule.Name.replace('.', '')}'
        self.OutputFilePath = f'{self.OrcaCachePath}\\{self.GetOutputFileName()}'
        self.InputFilePath = f'{self.OrcaCachePath}\\{self.GetInputFileName()}'

    def RunCalculation(self):
        """Runs a Orca Calculation in a Docker Container """

        # Make Cache Folder if it doesn't Exist
        if not os.path.exists(self.OrcaCachePath):
            os.makedirs(self.OrcaCachePath)

        # Make a folder for the Specific Calculation
        if not os.path.exists(self.OrcaCachePath):
            os.makedirs(self.OrcaCachePath)

        # Save the Input File to the folder
        self.SaveInputFile(self.OrcaCachePath)

        if (self.IsLocal):
            result = self.RunLocally(orcaCachePath)
        else:
            result = self.RunDockerContainer(orcaCachePath)
            
        if (result.stderr.__len__() > 0):
            print(f"WARNING Errors Maybe Occured : \n\n{result.stderr}")

        print(f"Calculation Complete : {self.GetInputFileName()}")

        # Open the Output File and Grab the Content
        with open(self.OutputFilePath, 'r') as file:
            self.CalculationOutput = file.read()
        
    def RunLocally (self, cachePath: str):
        # Create the Command String
        command = ""
        
        if os.name == 'nt':
            orcaPath = subprocess.run("where orca", shell=True, text=True, capture_output=True).stdout.strip(" \n\"").removesuffix(".exe")
            command = f"cd /d {cachePath} && \"{orcaPath}\" \"{self.GetInputFileName()}\" > \"{self.GetOutputFileName()}\""
        else:
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

    def GetInputFileName (self):
        """Returns the Input File Name with it's extension"""
        return f"{self.CalculationMolecule.Name}.inp"

    def GetOutputFileName (self):
        """Returns the Output File Name with it's extension"""
        return f"{self.CalculationMolecule.Name}.out"

    def GetInputFile (self):
        """Generates a Input file in String format"""

        firstLine = "!"

        # Check if the Properties are defined and Add them to the First line of the Input File
        if (self.DensityFunctional):
            firstLine += f"{self.DensityFunctional}"

        if (self.BasisSet):
            firstLine += f" {self.BasisSet}"

        if (self.CalculationType):
            firstLine += f" {self.CalculationType}"

        if (self.Cores > 1):
            if (self.Cores < 9):
                firstLine += " " + f"PAL{self.Cores}"
            else:
                firstLine += f"\n%PAL NPROCS {self.Cores} END"

        # Generate XYZ Wrappers
        xyzWrapperStart = "* xyz 0 1"
        xyzWrapperEnd = "*"

        # Create a XYZ File for the Molecule
        xyz = XYZFile(self.CalculationMolecule)

        # Generate the Entire Input File as a String
        inputFile = f"{firstLine}\n{xyzWrapperStart}\n{xyz.GetXYZBody()}\n{xyzWrapperEnd}"

        return inputFile

    def SaveInputFile (self, filePath: str):
        """Saves a Input File using the Settings Provided to the Path Specified"""
        with open(os.path.join(filePath, self.GetInputFileName()), "w") as file:
            file.write(self.GetInputFile())

class OrcaInputFile:
    # Template is an input file with missing variables in the form of &{variable_name}
    # kwargs is a dictionary with variable names as keys and their values as values
    # Example:
    # !SP &{basis} PBE
    # *xyzfile 0 1 aspirin.xyz
    # with variables={'basis': 'def2-SVP'}
    # will be converted to 
    # !SP def2-SVP PBE
    # *xyzfile 0 1 aspirin.xyz

    # Example with template file:
    # tester = OrcaInputFile(OrcaInputTemplate.BASIC
    #   , calculation='OPT'
    #   , basis='def2-SVP'
    #   , functional='PBE'
    #   , xyzfile='aspirin.xyz'
    #)
    
    def __init__(self, template: str, **variables):
        self.template = template
        self.variables = variables
        self.inputfile = self.GenerateInputFile()

    def GenerateInputFile(self) -> str:
        """Generates the input file content by replacing placeholders with actual values."""
        if isinstance(self.template, OrcaInputTemplate):
            input_content=self.template.value    
        else:
            with open(self.template, 'r') as file:
                input_content = file.read()
        
        for key, value in self.variables.items():
            placeholder = f'&{{{key}}}'
            input_content = input_content.replace(placeholder, str(value))
        return input_content

    def SaveInputFile(self, file_path: str):
        """Saves the generated input file content to a specified path."""
        with open(file_path, 'w') as file:
            file.write(self.inputfile)
