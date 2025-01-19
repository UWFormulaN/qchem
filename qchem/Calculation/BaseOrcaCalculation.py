import os
from typing import Any
from qchem.Data.Enums import OrcaInputTemplate
from .OrcaInputFile import OrcaInputFile
from ..Molecule import Molecule
from abc import ABC, abstractmethod

class BaseOrcaCalculation(ABC):

    molecule: Molecule | str
    """String Path to the XYZ file or Molecule Object of the Molecule to run a Calculation on"""

    template: str | OrcaInputTemplate
    """Template for the Input File that gets filled out to run the Calculation"""

    calculationType: str = ""
    """The Keyword for the Calculation to run on the Molecule (For Pipelines replace with name)"""

    inputFile: OrcaInputFile
    """The auto generated Input File that is used for the Calculation"""

    inputFilePath: str
    """The path the Input File will be saved to, often orcaCachePath + inputFileName.inp"""

    outputFilePath: str
    """The path the Output File will be saved to, often orcaCachePath + outputFileName.out"""

    variables: dict[str, Any]
    """Additional Variables that will be pasted into the Input File Template"""

    index: int
    """Number to identify individual Docker Orca Calculations running in parallel"""

    name: str
    """Friendly name for the calculation, used for Input and Output file names"""

    orcaCachePath: str
    """The path to the folder that stores temporary and resulting Calculation Files"""

    calculationTime: float
    """The resulting total time for the calculation in seconds."""

    STDOut: bool
    """A Boolean flag to indicate if Standard Output logs should be printed"""

    isLocal: bool
    """A Boolean flag to indicate if the calculation runs locally or in Docker (True = Local, False = Docker)"""

    cores: int
    """Number of CPU Cores allocated to the calculation"""

    defaultName: str = "Molecule"
    """Default Calculation Name to use if unspecified. Will check if Molecule Object already has a name first."""

    def __init__(
        self,
        name: str,
        molecule: str | Molecule,
        template: str | OrcaInputTemplate,
        index: int,
        cores: int,
        isLocal: bool,
        stdout: bool,
        **variables,
    ):

        # Value Type Checking
        if not (self.calculationType and isinstance(self.calculationType, (str))):
            raise ValueError("Calculation Type must be defined and a String")

        if not isinstance(name, (str)) or name == "":
            raise ValueError("Name of the Calculation must be specified ")

        if not (molecule and isinstance(molecule, (str, Molecule))):
            raise ValueError("Molecule must be a Molecule Object or String")

        if not isinstance(template, (str, OrcaInputTemplate)):
            raise ValueError("Template must be a String or OrcaInputTemplate")

        if not isinstance(index, (int)):
            raise ValueError("Index must be an integer")

        if not isinstance(cores, (int)):
            raise ValueError("Cores must be an integer")

        if not isinstance(isLocal, (bool)):
            raise ValueError("IsLocal must be a boolean")

        if not isinstance(stdout, (bool)):
            raise ValueError("STDOut must be a boolean")

        # Set the Values
        self.name = name
        self.molecule = molecule
        self.template = template
        self.index = index
        self.cores = cores
        self.isLocal = isLocal
        self.STDOut = stdout
        self.variables = variables

        # Set Name to the Name of the Molecule
        if self.name == self.defaultName:
            if isinstance(molecule, Molecule):
                self.name = molecule.Name
            else:
                self.name = (self.calculationType + " Molecule").replace(" ", "_")

        # Set Appropriate XYZ Format and Default Templates (Only if unspecified)
        if self.isFileReference():
            self.variables["xyzfile"] = self.molecule
            if self.template == "":
                if self.cores == 1:
                    self.template = OrcaInputTemplate.BASIC
                else:
                    self.template = OrcaInputTemplate.BASICPARALLEL
        else:
            self.variables["xyz"] = self.molecule.XYZBody()
            if self.template == "":
                if self.cores == 1:
                    self.template = OrcaInputTemplate.BASICXYZ
                else:
                    self.template = OrcaInputTemplate.BASICXYZPARALLEL

        # Define the Calculation Type and the Number of Cores to Use
        self.variables["calculation"] = self.calculationType
        self.variables["cores"] = self.cores

        # Generate Cache Paths
        orcaCache = "OrcaCache"  # Change this to QChemCache?
        self.orcaCachePath = os.path.join(os.getcwd(), orcaCache, self.name)

        # Convert to Functions?
        self.outputFilePath = os.path.join(self.orcaCachePath, self.getOutputFileName())
        self.inputFilePath = os.path.join(self.orcaCachePath, self.getInputFileName())

        # Create the Input File
        self.inputFile = OrcaInputFile(self.template, **self.variables)

    @abstractmethod
    def runCalculation(self):
        """Runs the Calculation Algorithm, is an Abstract Method/Function that needs to be overridden.

        Parameters:\n
            self - Default Parameter for the Class Instance

        Returns: None - No Return Value
        """
        pass

    def isFileReference(self):
        """Checks if the Molecule is defined as a File Reference as a str or a Molecule Object

        Parameters : \n
            self - Default Parameter for the Class Instance

        Returns : bool - True if the Molecule is a File Reference, False if it is a Molecule Object
        """
        if isinstance(self.molecule, (str)):
            return True
        else:
            return False

    def getInputFileName(self):
        """Returns the Input File Name with it's extension"""
        return f"{self.name}.inp"

    def getOutputFileName(self):
        """Returns the Output File Name with it's extension"""
        return f"{self.name}.out"

    def getOutput(self) -> str:
        """Returns the Entire Content of the Output File as a Single String"""
        # Open the Output File and Grab the Content
        with open(self.outputFilePath, "r") as file:
            self.CalculationOutput = file.read()

    def createDirectories(self):
        """Creates the necessary Cache Folders to store the Calculation"""
        # Make Cache Folder if it doesn't Exist
        if not os.path.exists(self.orcaCachePath):
            os.makedirs(self.orcaCachePath)

        # Make a folder for the Specific Calculation
        if not os.path.exists(self.orcaCachePath):
            os.makedirs(self.orcaCachePath)

    def clockTime(self, seconds):
        """Converts Seconds to a Human Readable Time String"""
        # Convert Seconds to Hours, Minutes, and Seconds
        days = seconds // 86400
        hours = (seconds % 86400) // 3600
        minutes = (seconds % 3600) // 60
        remainingSeconds = seconds % 60

        # Generate the Time String
        parts = []
        if days > 0:
            parts.append(f"{int(hours)} day{'s' if days > 1 else ''}")
        if hours > 0:
            parts.append(f"{int(hours)} hour{'s' if hours > 1 else ''}")
        if minutes > 0:
            parts.append(f"{int(minutes)} minute{'s' if minutes > 1 else ''}")
        if remainingSeconds > 0:
            parts.append(
                f"{int(remainingSeconds)} second{'s' if remainingSeconds > 1 else ''}"
            )

        # Return the Time String
        return ", ".join(parts) if parts else "0 seconds"

    def basisSetFunctionalCompliant(self):
        """Checks if the Basis Set and Functional Values have been defined. Used to simply call and confirm in Calculation Classes that require a Basis Set and Functional to be defined"""
        if not ("basis" in self.variables):
            raise ValueError("BasisSet not defined! Provide Basis Set Name as a String")

        if not ("functional" in self.variables):
            raise ValueError(
                "Functional not defined! Provide Functional Name as a String"
            )

    def setCalculationType(self, calculationType):
        """Sets the Calculation Type in the Variables Dictionary"""
        self.variables["calculation"] = calculationType
