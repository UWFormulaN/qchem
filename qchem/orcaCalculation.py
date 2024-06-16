from qchem.Enums.CalculationTypes import CalculationType
from qchem.molecule import Molecule
from qchem.Enums.BasisSet import BasisSet 


class OrcaCalculation:

    FINISHED = False

    Molecule: Molecule

    Type : CalculationType
    
    BasisSet: BasisSet

    def RunCalculation(self):

        # Create a "Process" like in Java Script
        # Have it create a new Docker Container using the Orca Image
        # Transfer the XYZ File somehow (Save and then transfer?)

        print("Running Calculation")


        