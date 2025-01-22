import time
import pandas as pd
from qchem.Molecule import Molecule
from qchem.Parser import OrcaOutput
from .BaseOrcaCalculation import BaseOrcaCalculation
from qchem.Calculation.OrcaCalculation import runOrcaCalculation
from qchem.Data.Enums import OrcaCalculationType, OrcaInputTemplate


class Frequency(BaseOrcaCalculation):
    """Performs a Frequency Calculation on a Molecule. Will Expose IR and Vibrational Frequencies"""

    #
    # Need to be Set
    #
    calculationType: str = OrcaCalculationType.FREQUENCY.value
    """The Keyword for the Calculation to run on the Molecule (For Pipelines replace with name)"""

    vibrationalFrequencies: pd.DataFrame
    """The resulting Vibrational Frequencies of the Molecule"""

    IRFrequencies: pd.DataFrame
    """The resulting Infra Red Frequencies of the Molecule"""

    def __init__(
        self,
        molecule: str | Molecule,
        template: str | OrcaInputTemplate = "",
        index: int = 1,
        cores: int = 1,
        isLocal: bool = False,
        name: str = "Molecule",
        stdout: bool = True,
        **variables,
    ):
        # Make a Super Call (Use the Base Class Init for some Boilerplate Setup)
        super().__init__(
            name, molecule, template, index, cores, isLocal, stdout, **variables
        )

        # Check if the Calculation has a Basis Set and a Functional Defined (Specific to Certain Calculations)
        self.basisSetFunctionalCompliant()

    def runCalculation(self):
        """Runs the Frequency Calculation and Saves the Vibrational Frequencies and Infra Red Frequencies
        
        ## Parameters: \n
            self - Default Parameter for the Class Instance
            
        ## Returns: \n
            None - No Return Value
        """
        # Start the Clock
        startTime = time.time()

        # Add a Print Statement to say we are running
        print(f"Running FREQ on {self.name}...")

        # Run the Orca Calculation
        calculation = runOrcaCalculation(
            self.name, self.inputFile, self.index, self.isLocal, STDOut=False
        )

        # Get the Calculation Time
        self.calculationTime = time.time() - startTime

        # Save the Output File Path
        self.outputFilePath = calculation.outputFilePath
        self.orcaCachePath = calculation.orcaCachePath

        # Load the Output File
        outputFile = OrcaOutput(calculation.outputFilePath)

        # Extract the Vibrational Frequency from the Output File
        self.vibrationalFrequencies = outputFile.getVibrationalFrequencies()

        # Load the IR Frequency from the
        self.IRFrequencies = outputFile.getIRFrequencies()

        # Display a Print Statement for the Frequency Completion
        print(f"Finished FREQ on {self.name}! ({self.clockTime(self.calculationTime)})")