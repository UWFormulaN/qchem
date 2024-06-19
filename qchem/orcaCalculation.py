from qchem.Enums.OrcaCalculationTypes import OrcaCalculationType
from qchem.Enums.OrcaDensityFunctional import OrcaDensityFunctional
from qchem.molecule import Molecule
from qchem.Enums.OrcaBasisSet import OrcaBasisSet 

class OrcaCalculation:
    """Class capable of running an Orca Calculation"""

    FINISHED = False

    Molecule: Molecule
    """Molecule that will have an Orca Calculation run on it"""

    CalculationType : OrcaCalculationType
    """The Type of Calculation that will occur on the Molecule"""
    
    BasisSet: OrcaBasisSet
    """Basis Set to use for the Calculation"""

    DensityFunctional : OrcaDensityFunctional
    """The Density Functional to use for the Calculation"""

    Cores : int
    """Number of Cores to use for the Calculation"""

    def __init__(self, calculationType: OrcaCalculationType):

        self.CalculationType = calculationType

    def RunCalculation(self):

        # Create a "Process" like in Java Script
        # Have it create a new Docker Container using the Orca Image
        # Transfer the XYZ File somehow (Save and then transfer?)


        # Density Function, Basis Set, Calculation Type, PAL#, XYZ File

        print("Running Calculation")

    def GenerateInputFile (self):

        firstLine = "!"

        if (self.DensityFunctional):
            firstLine += self.DensityFunctional

        if (self.BasisSet):
            firstLine += " " + self.BasisSet

        if (self.CalculationType):
            firstLine += " " + self.CalculationType

        if (self.BasisSet > 1):
            firstLine += " " + f"PAL{self.Cores}"

        xyzWrapperStart = "* xyz 0 1"
        xyzWrapperEnd = "*"

        




    


    


        