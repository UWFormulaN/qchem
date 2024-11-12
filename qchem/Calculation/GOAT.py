import time
import os
from qchem.Molecule import Molecule
from qchem.XYZFile import XYZFile
from qchem.Parser import OrcaOutput
from qchem.Calculation.OrcaInputFile import OrcaInputFile
from qchem.Calculation.OrcaCalculation import OrcaCalculation
from qchem.Data.Enums import OrcaInputTemplate, OrcaCalculationType

class GOAT:
    
    # Most of this is Boilerplate, Should make a Abstract class to inherit this boilerplate logic
    
    molecule: Molecule | str
    """The Molecule to be Optimized"""
    
    cores: int
    """The Number of Cores the Geometry Optimization will Utilize"""
    
    isLocal : bool
    """Boolean Flag to determine if the Optimization is Local or not (Local = Runs on Device, Non-Local = Runs in Docker Container)"""
    
    name: str
    """Name of the Molecule being GeoOptimized"""
    
    calculationTime: float
    """Time Elapsed for the GOAT Optimization to Complete"""
    
    orcaCachePath : str
    """The Path to the Orca Cache folder for the Calculation"""
    
    outputFilePath : str
    """The Path to the Output File"""
    
    conformers: list[Molecule]
    """List of Conformer Molecules Calculated from GOAT"""
    
    conformerContribution : list[float]
    """A List of the Contributions each Conformer provides to the Ensemble"""
    
    def __init__(self, molecule: str | Molecule, cores:int = 1, isLocal:bool = False, name:str = ""):
        
        # Check if Values are empty or of Wrong Type
        if (not (molecule and isinstance(molecule, (str, Molecule)))):
            raise ValueError("Molecule is not defined! Provide a Path to the XYZ file or a Molecule Object")
        
        if (not isinstance(cores, (int))):
            raise ValueError("Cores is not an Integer! Provide an Integer Value")
            
        if (not isinstance(isLocal, (bool))):
            raise ValueError("IsLocal is not a Boolean! Provide a Boolean Value")
        
        if (not isinstance(name, (str))):
            raise ValueError("Name is not a String! Provide a String Value")
        
        # Set the Name of the Molecule
        if (name == ""):
            if (isinstance(molecule, (Molecule))):
                self.name = molecule.Name
            else:
                print("No Name Provided for the Molecule, using Default Name: GOATMolecule")
                self.name = "GOATMolecule"
        else:
            self.name = name
        
        # Set the Values
        self.molecule = molecule
        self.cores = cores
        self.isLocal = isLocal
        self.conformers = []
        
    def IsFileReference(self):
        """Determines if the Molecule is stored as a File Reference, or is a direct Molecule Object"""
        if (isinstance(self.molecule, (str))):
            return True
        else:
            return False
        
    def RunCalculation (self):
        """Runs through the GOAT Calculation"""
        
        #GOATTemplate: OrcaInputTemplate | str
        calculationType = f"GOAT XTB PAL{self.cores}"
        
        # Start the Clock
        startTime = time.time()
        
        # Check if we are using a File reference, create the appropriate Input File associated with it
        if (self.IsFileReference()):
            inputFile = OrcaInputFile(OrcaInputTemplate.BASIC,
                                      calculation = calculationType,
                                      basis = "",
                                      functional = "",
                                      xyzfile = self.molecule)
            
        else:
            inputFile = OrcaInputFile(OrcaInputTemplate.BASICXYZ,
                                      calculation = calculationType,
                                      basis = "",
                                      functional = "",
                                      xyz = self.molecule.XYZBody())
        
        # Create the Calculation Object
        calculation = OrcaCalculation(self.name, inputFile, isLocal=self.isLocal, stdout=False)
        
        # Add a Print Statement to say we are running
        print(f"Running GOAT (Global Optimizer Algorithm) on {self.name}...")
        
        # Run the Calculation
        calculation.RunCalculation()
        
        # Get the Calculation Time
        self.calculationTime = time.time() - startTime
    
        # Save the Output File Path
        self.outputFilePath = calculation.OutputFilePath
        self.orcaCachePath = calculation.OrcaCachePath
        
        # Extract the Conformer Molecules from the Resulting XYZ File
        self.ExtractConformers()
        
        # Extract the Contributions
        conformerOutput = OrcaOutput(calculation.OutputFilePath).conformers
        self.conformerContribution = conformerOutput[conformerOutput.columns[3]].values
        
        # Display a Print Statement for the GOAT Completion
        print(f"Finished GOAT (Global Optimizer Algorithm) on {self.name}! ({self.ClockTime(self.calculationTime)})")
    
    def ExtractConformers (self):
        """Extracts all the Conformer Molecules"""
        # Get the Number of Atoms we should Expect
        if (self.IsFileReference()):
            atomNum = XYZFile(self.molecule).AtomCount
        else:
            atomNum = self.molecule.AtomCount
        
        # Open the File
        ensembleXYZFile = open(os.path.join(self.orcaCachePath, f"{self.name}.finalensemble.xyz"))
        
        # Get all the Lines from the File
        allLines = ensembleXYZFile.readlines()
        
        # Get the Expected Length of a XYZ File
        XYZLength = atomNum + 2
        
        # Calculate the Number of 
        moleculeCount = int((len(allLines))/(XYZLength))
        
        for i in range(moleculeCount):
            # Get the Lines for a Single XYZ File
            molLines = allLines[i*XYZLength:(i+1)*XYZLength:]
            
            # Make the Name
            moleculeName = f"{self.name}_Conf_{i}"
            
            # Load as a XYZ File
            xyz = XYZFile(molecule=molLines, name=moleculeName)
            
            # Convert to a Molecule
            molecule = Molecule(moleculeName, xyz)
            
            # Add to the Conformer List
            self.conformers.append(molecule)
            
    def ClockTime(self, seconds):
        """Converts Seconds to a Human Readable Time String"""
        # Convert Seconds to Hours, Minutes, and Seconds
        hours = seconds // 3600
        minutes = (seconds % 3600) // 60
        remainingSeconds = seconds % 60
        
        # Generate the Time String
        parts = []
        if hours > 0:
            parts.append(f"{int(hours)} hour{'s' if hours > 1 else ''}")
        if minutes > 0:
            parts.append(f"{int(minutes)} minute{'s' if minutes > 1 else ''}")
        if remainingSeconds > 0:
            parts.append(f"{int(remainingSeconds)} second{'s' if remainingSeconds > 1 else ''}")
        
        # Return the Time String
        return ", ".join(parts) if parts else "0 seconds"