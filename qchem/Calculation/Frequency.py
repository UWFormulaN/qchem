import time
from typing import Any
import pandas as pd
from qchem.Molecule import Molecule
from qchem.Parser import OrcaOutput
from qchem.Calculation.OrcaInputFile import OrcaInputFile
from qchem.Data.Enums import OrcaCalculationType, OrcaInputTemplate
from .BaseOrcaCalculation import BaseOrcaCalculation
from .OrcaCalcs import RunOrcaCalculation

class Frequency(BaseOrcaCalculation):

    #
    # Need to be Set
    #
    calculationType: str = OrcaCalculationType.FREQUENCY.value

    vibrationalFrequencies: pd.DataFrame
    """Vibrational Frequencies of the Molecule"""

    IRFrequencies: pd.DataFrame
    """IR Frequencies of the Molecule"""
    
    def __init__(
        self,
        molecule: str | Molecule,
        template: str | OrcaInputTemplate = "",
        index: int = 1,
        cores: int = 1,
        isLocal: bool = False,
        name: str = "FREQMolecule",
        stdout: bool = True,
        **variables
    ):
        # Make a Super Call (Use the Base Class Init for some Boilerplate Setup)
        super().__init__(name, molecule, template, index, cores, isLocal, stdout, **variables)

        # Check if the Calculation has a Basis Set and a Functional Defined (Specific to Certain Calculations)
        self.BasisSetFunctionalCompliant()
        
    def RunCalculation(self):
        """Runs through the GOAT Calculation"""

        # Start the Clock
        startTime = time.time()
        
        # Check if we are using a File reference, create the appropriate Input File associated with it
        #if self.IsFileReference():
        #    inputFile = OrcaInputFile(
        #        OrcaInputTemplate.BASIC,
        #        calculation=self.calculationType,
        #        #basis=self.basisSet,
        #        #functional=self.functional,
        #        xyzfile=self.molecule,
        #        self.idk
        #    )
        #else:
        #    inputFile = OrcaInputFile(
        #        OrcaInputTemplate.BASICXYZ,
        #        calculation=self.calculationType,
        #        #basis=self.basisSet,
        #        #functional=self.functional,
        #        xyz=self.molecule.XYZBody(),
        #        variables = self.idk
        #    )

        # Add a Print Statement to say we are running
        print(f"Running Frequency Analysis on {self.name}...")

        # Run the Orca Calculation
        calculation = RunOrcaCalculation(
            self.name, self.inputFile, self.index, self.isLocal, self.STDOut
        )

        # Get the Calculation Time
        self.calculationTime = time.time() - startTime

        # Save the Output File Path
        self.outputFilePath = calculation.outputFilePath
        self.orcaCachePath = calculation.orcaCachePath

        # Load the Output File
        outputFile = OrcaOutput(calculation.outputFilePath)

        # Extract the Vibrational Frequency from the Output File
        self.vibrationalFrequencies = outputFile.get_vibrational_frequencies()

        # Load the IR Frequency from the
        self.IRFrequencies = outputFile.get_ir_frequencies()

        # Display a Print Statement for the Frequency Completion
        print(
            f"Finished {self.calculationType} on {self.name}! ({self.ClockTime(self.calculationTime)})"
        )