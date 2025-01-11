import os
import time
from qchem.XYZFile import XYZFile
from qchem.Molecule import Molecule
from qchem.Parser import OrcaOutput
from qchem.Data.Enums import OrcaInputTemplate
from qchem.Calculation.OrcaCalculation import RunOrcaCalculation
from qchem.Calculation.BaseOrcaCalculation import BaseOrcaCalculation

class GOAT(BaseOrcaCalculation):

    #
    # Need to be Set
    #
    calculationType: str = "GOAT XTB"

    conformers: list[Molecule]
    """List of Conformer Molecules Calculated from GOAT"""

    conformerContribution: list[float]  # Switch to DataFrame?
    """A List of the Contributions each Conformer provides to the Ensemble"""

    def __init__(
        self,
        molecule: str | Molecule,
        template: str | OrcaInputTemplate = "",
        index: int = 1,
        cores: int = 1,
        isLocal: bool = False,
        name: str = "GOATMolecule",
        stdout: bool = True,
        **variables
    ):
        # Set Basis and Functional to Empty
        variables["basis"] = ""
        variables["functional"] = ""

        # Make a Super Call (Use the Base Class Init for some Boilerplate Setup)
        super().__init__(name, molecule, template, index, cores, isLocal, stdout, **variables)

        # Set the Values
        self.conformers = []

    def RunCalculation(self):
        """Runs through the GOAT Calculation"""

        # Start the Clock
        startTime = time.time()

        # Add a Print Statement to say we are running
        print(f"Running GOAT on {self.name}...")

        # Create the Calculation Object
        calculation = RunOrcaCalculation(self.name, self.inputFile, self.index, self.isLocal, STDOut=False)
        
        # Get the Calculation Time
        self.calculationTime = time.time() - startTime

        # Save the Output File Path
        self.outputFilePath = calculation.outputFilePath
        self.orcaCachePath = calculation.orcaCachePath

        # Extract the Conformer Molecules from the Resulting XYZ File
        self.ExtractConformers()

        # Extract the Contributions
        conformerOutput = OrcaOutput(calculation.outputFilePath).conformers
        self.conformerContribution = conformerOutput[conformerOutput.columns[3]].values

        # Display a Print Statement for the GOAT Completion
        print(
            f"Finished GOAT on {self.name}! ({self.ClockTime(self.calculationTime)})"
        )

    def ExtractConformers(self):
        """Extracts all the Conformer Molecules"""
        # Get the Number of Atoms we should Expect
        if self.IsFileReference():
            atomNum = XYZFile(os.path.join(self.orcaCachePath, self.molecule)).AtomCount
        else:
            atomNum = self.molecule.AtomCount

        # Open the File
        ensembleXYZFile = open(
            os.path.join(self.orcaCachePath, f"{self.name}.finalensemble.xyz")
        )

        # Get all the Lines from the File
        allLines = ensembleXYZFile.readlines()

        # Get the Expected Length of a XYZ File
        XYZLength = atomNum + 2

        # Calculate the Number of
        moleculeCount = int((len(allLines)) / (XYZLength))

        for i in range(moleculeCount):
            # Get the Lines for a Single XYZ File
            molLines = allLines[i * XYZLength : (i + 1) * XYZLength :]

            # Make the Name
            moleculeName = f"{self.name}_Conf_{i}"

            # Load as a XYZ File
            xyz = XYZFile(molecule=molLines, name=moleculeName)

            # Convert to a Molecule
            molecule = Molecule(moleculeName, xyz)

            # Add to the Conformer List
            self.conformers.append(molecule)