from os import path
import os
import subprocess
import comm
from distutils import core
from qchem.Enums.OrcaCalculationTypes import OrcaCalculationType
from qchem.Enums.OrcaDensityFunctional import OrcaDensityFunctional
from qchem.XYZFile import XYZFile
from qchem.molecule import Molecule
from qchem.Enums.OrcaBasisSet import OrcaBasisSet 

class OrcaCalculation:
    """Class capable of running an Orca Calculation"""

    FINISHED = False

    CalculationOutput: str

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

    def __init__(self, molecule: Molecule, calculationType: OrcaCalculationType = None, basisSet: OrcaBasisSet = None, densityFunctional: OrcaDensityFunctional = None, cores : int = 1 ):

        self.CalculationMolecule = molecule
        self.CalculationType = calculationType
        self.BasisSet = basisSet
        self.DensityFunctional = densityFunctional
        self.Cores = cores

    def RunCalculation(self):

        # Create a "Process" like in Java Script
        # Have it create a new Docker Container using the Orca Image
        # Transfer the XYZ File somehow (Save and then transfer?)

        # Density Function, Basis Set, Calculation Type, PAL#, XYZ File

        # process = subprocess.Popen(['docker', 'run', '-it', 'mrdnalex/orca'])

        # print(process.stdout)

        orcaCache = "OrcaCache"

        if not os.path.exists(orcaCache):
            os.makedirs(orcaCache)

        path = f'{os.getcwd()}\\{orcaCache}'

        print(f"Running Calulation {self.GetInputFileName()}")

        self.SaveInputFile(orcaCache)

        command = f'docker run --name qchemorca -v "{path}":/home/orca mrdnalex/orca sh -c "cd /home/orca && /Orca/orca {self.GetInputFileName()} > {self.GetOutputFileName()}"'

        result = subprocess.run(command, shell=True, text=True, capture_output=True)

        subprocess.run("docker kill qchemorca" , shell=True)
        subprocess.run("docker rm qchemorca" , shell=True)

        print(f"Calculation {self.GetInputFileName()} Complete")

        orcaOutput = result.stdout

        with open(f'{os.getcwd()}/{orcaCache}/{self.GetOutputFileName()}', 'r') as file:
            self.CalculationOutput = file.read()

        self.FINISHED = True

        # print(result.stdout)

    def GetInputFileName (self):
        return f"{self.CalculationMolecule.Name}.inp"

    def GetOutputFileName (self):
        return f"{self.CalculationMolecule.Name}.out"

    def GetInputFile (self):
        """Generates a Input file in String format"""

        firstLine = "!"

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

        xyzWrapperStart = "* xyz 0 1"
        xyzWrapperEnd = "*"

        xyz = XYZFile(self.CalculationMolecule)

        inputFile = f"{firstLine}\n{xyzWrapperStart}\n{xyz.GetXYZBody()}\n{xyzWrapperEnd}"

        return inputFile
    

    def SaveInputFile (self, filePath: str):
        """Saves a Input File using the Settings Provided to the Path Specified"""
        with open(os.path.join(filePath, self.GetInputFileName()), "w") as file:
            file.write(self.GetInputFile())
        




    


    


        