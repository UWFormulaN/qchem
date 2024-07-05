from os import path
import os
import subprocess
from textwrap import indent
import comm
from distutils import core
from qchem.Enums.OrcaCalculationTypes import OrcaCalculationType
from qchem.Enums.OrcaDensityFunctional import OrcaDensityFunctional
from qchem.XYZFile import XYZFile
from qchem.Molecule import Molecule
from qchem.Enums.OrcaBasisSet import OrcaBasisSet 

class OrcaCalculation:
    """Class capable of running an Orca Calculation"""

    CalculationOutput: str
    """The Entire Orca Output File in a single String"""

    CalculationMolecule: Molecule
    """Molecule that will have an Orca Calculation run on it"""

    CalculationType : OrcaCalculationType
    """The Type of Calculation that will occur on the Molecule"""
    
    BasisSet: OrcaBasisSet
    """Basis Set to use for the Calculation"""

    DensityFunctional : OrcaDensityFunctional
    """The Density Functional to use for the Calculation"""

    Cores : int
    """Number of Cores to use for the Calculation"""

    OutputFilePath: str
    """The Path to the Output File on the Device"""

    InputFilePath: str
    """The Path to the Input File on the Device"""

    Index: int
    """Container Index so that they can be run in parallel"""

    def __init__(self, molecule: Molecule, calculationType: OrcaCalculationType = None, basisSet: OrcaBasisSet = None, densityFunctional: OrcaDensityFunctional = None, cores : int = 1):

        self.CalculationMolecule = molecule
        self.CalculationType = calculationType
        self.BasisSet = basisSet
        self.DensityFunctional = densityFunctional
        self.Cores = cores
        self.Index = 1

    def RunCalculation(self):
        """Runs a Orca Calculation in a Docker Container """
        orcaCache = "OrcaCache"
        orcaCachePath = f'{os.getcwd()}\\{orcaCache}\\{self.CalculationMolecule.Name.replace('.', '')}'
        self.OutputFilePath = f'{orcaCachePath}\\{self.GetOutputFileName()}'
        self.InputFilePath = f'{orcaCachePath}\\{self.GetInputFileName()}'

        # Make Cache Folder if it doesn't Exist
        if not os.path.exists(orcaCache):
            os.makedirs(orcaCache)

        # Make a folder for the Specific Calculation
        if not os.path.exists(orcaCachePath):
            os.makedirs(orcaCachePath)

        # Save the Input File to the folder
        self.SaveInputFile(orcaCachePath)

        # Create the Command String
        command = f'docker run --name qchemorca{self.Index} -v "{orcaCachePath}":/home/orca mrdnalex/orca sh -c "cd /home/orca && /Orca/orca {self.GetInputFileName()} > {self.GetOutputFileName()}"'

        print(f"Running Calulation : {self.GetInputFileName()}")

        # Kill and Remove qchemorca container if it doesn't exist yet
        subprocess.run(f"docker kill qchemorca{self.Index}" , shell=True)
        subprocess.run(f"docker rm qchemorca{self.Index}" , shell=True)

        # Run the Calculation in a Container and wait
        subprocess.run(command, shell=True, text=True, capture_output=True)

        # Kill and Remove the Container
        subprocess.run(f"docker kill qchemorca{self.Index}" , shell=True)
        subprocess.run(f"docker rm qchemorca{self.Index}" , shell=True)

        print(f"Calculation Complete : {self.GetInputFileName()}")

        # Open the Output File and Grab the Content
        with open(self.OutputFilePath, 'r') as file:
            self.CalculationOutput = file.read()

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
            firstLine += f"{self.DensityFunctional.value}"

        if (self.BasisSet):
            firstLine += f" {self.BasisSet.value}"

        if (self.CalculationType):
            firstLine += f" {self.CalculationType.value}"

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