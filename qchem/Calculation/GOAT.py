import time
import os
from qchem.Molecule import Molecule
from qchem.XYZFile import XYZFile
from qchem.Parser import OrcaOutput
from qchem.Calculation.OrcaCalculation import OrcaCalculation
from qchem.Data.Enums import OrcaInputTemplate
from qchem.Calculation.BaseOrcaCalculation import BaseOrcaCalculation
from .OrcaCalcs import RunOrcaCalculation

class GOAT(BaseOrcaCalculation):

    #
    # Need to be Set
    #
    calculationType: str = "GOAT XTB"

    # Most of this is Boilerplate, Should make a Abstract class to inherit this boilerplate logic
    #molecule: Molecule | str
    """The Molecule to be Optimized"""

    #cores: int
    """The Number of Cores the Geometry Optimization will Utilize"""

    #isLocal: bool
    """Boolean Flag to determine if the Optimization is Local or not (Local = Runs on Device, Non-Local = Runs in Docker Container)"""

    #name: str
    """Name of the Molecule being GeoOptimized"""

    #calculationTime: float
    """Time Elapsed for the GOAT Optimization to Complete"""

    #orcaCachePath: str
    """The Path to the Orca Cache folder for the Calculation"""

    #outputFilePath: str
    """The Path to the Output File"""

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
        variables["basis"] = ""
        variables["functional"] = ""

        # Make a Super Call (Use the Base Class Init for some Boilerplate Setup)
        super().__init__(name, molecule, template, index, cores, isLocal, stdout, **variables)

        # Set the Values
        self.conformers = []

    def RunCalculation(self):
        """Runs through the GOAT Calculation"""

        # GOATTemplate: OrcaInputTemplate | str
        #calculationType = f"GOAT XTB PAL{self.cores}"

        # Start the Clock
        startTime = time.time()

        # Check if we are using a File reference, create the appropriate Input File associated with it
        #if self.IsFileReference():
        #    inputFile = OrcaInputFile(
        #        OrcaInputTemplate.BASIC,
        #        calculation=calculationType,
        #        basis="",
        #        functional="",
        #        xyzfile=self.molecule,
        #    )
#
        #else:
        #    inputFile = OrcaInputFile(
        #        OrcaInputTemplate.BASICXYZ,
        #        calculation=calculationType,
        #        basis="",
        #        functional="",
        #        xyz=self.molecule.XYZBody(),
        #    )

        # Add a Print Statement to say we are running
        print(f"Running GOAT (Global Optimizer Algorithm) on {self.name}...")

        # Create the Calculation Object
        calculation = RunOrcaCalculation(self.name, self.inputFile, self.index, self.isLocal, STDOut=False)
        
        #calculation = OrcaCalculation(
        #    self.name, self.inputFile, self.index, isLocal=self.isLocal, stdout=False, **self.variables
        #)

        # Run the Calculation
        #calculation.RunCalculation()

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