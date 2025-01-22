import os
import time
from qchem.XYZFile import XYZFile
from qchem.Molecule import Molecule
from qchem.Parser import OrcaOutput
from qchem.Data.Enums import OrcaInputTemplate, OrcaCalculationType
from qchem.Calculation.OrcaCalculation import runOrcaCalculation
from qchem.Calculation.BaseOrcaCalculation import BaseOrcaCalculation


class GOAT(BaseOrcaCalculation):
    """Performs a GOAT (Global Optimizer Algorithm) Calculation on a Molecule. Exposes the Conformers and their Contributions to the Ensemble"""

    #
    # Need to be Set
    #
    calculationType: str = OrcaCalculationType.GOAT_XTB.value
    """The Keyword for the Calculation to run on the Molecule (For Pipelines replace with name)"""

    conformers: list[Molecule]
    """List of the Conformer Molecules Found by GOAT"""

    conformerContribution: list[float]  # Switch to DataFrame?
    """List of the Percentage Contributions each Molecule has to the Ensemble"""

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
        # Set Basis and Functional to Empty
        variables["basis"] = ""
        variables["functional"] = ""

        # Make a Super Call (Use the Base Class Init for some Boilerplate Setup)
        super().__init__(
            name, molecule, template, index, cores, isLocal, stdout, **variables
        )

        # Set the Values
        self.conformers = []

    def runCalculation(self):
        """Runs the GOAT Calculation and Saves the list of Conformer Molecules and their Individual Contributions to the Ensemble
        
        ## Parameters : \n
            self - Default Parameter for the Class Instance
        
        ## Returns : \n
            None - No Return Value
        """

        # Start the Clock
        startTime = time.time()

        # Add a Print Statement to say we are running
        print(f"Running GOAT on {self.name}...")

        # Create the Calculation Object
        calculation = runOrcaCalculation(
            self.name, self.inputFile, self.index, self.isLocal, STDOut=False
        )

        # Get the Calculation Time
        self.calculationTime = time.time() - startTime

        # Save the Output File Path
        self.outputFilePath = calculation.outputFilePath
        self.orcaCachePath = calculation.orcaCachePath

        # Extract the Conformer Molecules from the Resulting XYZ File
        self.extractConformers()

        # Extract the Contributions
        conformerOutput = OrcaOutput(calculation.outputFilePath).conformers
        self.conformerContribution = conformerOutput[conformerOutput.columns[3]].values

        # Display a Print Statement for the GOAT Completion
        print(f"Finished GOAT on {self.name}! ({self.clockTime(self.calculationTime)})")

    def extractConformers(self):
        """Extracts all the Conformer Molecules from the .finalensemble.xyz File
        
        ## Parameters : \n
            self - Default Parameter for the Class Instance
            
        ## Returns : \n
            None - No Return Value
        """
        # Get the Number of Atoms we should Expect
        if self.isFileReference():
            atomNum = XYZFile(os.path.join(self.orcaCachePath, self.molecule)).atomCount
        else:
            atomNum = self.molecule.atomCount

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