import time
import pandas as pd
from qchem.Molecule import Molecule
from qchem.Parser import OrcaOutput
from qchem.Calculation.OrcaInputFile import OrcaInputFile
from qchem.Calculation.OrcaCalculation import OrcaCalculation
from qchem.Data.Enums import OrcaInputTemplate

class Frequency:
    
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
    
    basisSet: str
    """The Basis Set to be used for the Optimization"""
    
    functional: str
    """The Density Functional to be used for the Optimization"""
    
    vibrationalFrequencies : pd.DataFrame
    """Vibrational Frequencies of the Molecule"""
    
    IRFrequencies : pd.DataFrame
    """IR Frequencies of the Molecule"""
    
    def __init__(self, molecule: str | Molecule, basisSet: str, functional: str, cores:int = 1, isLocal:bool = False, name:str = ""):
        
        # Check if Values are empty or of Wrong Type
        if (not (molecule and isinstance(molecule, (str, Molecule)))):
            raise ValueError("Molecule is not defined! Provide a Path to the XYZ file or a Molecule Object")
        
        if (not (basisSet and isinstance(basisSet, (str)))):
            raise ValueError("BasisSet is not defined! Provide the Name of the Basis Set as a String")
        
        if (not functional and isinstance(functional, (str))):
            raise ValueError("Functional is not defined! Provide the Name of the Functional as a String")
        
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
        self.basisSet = basisSet
        self.functional = functional
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
        
        # Define the Calculation Type
        calculationType = f"FREQ"
        
        # Start the Clock
        startTime = time.time()
        
        # Check if we are using a File reference, create the appropriate Input File associated with it
        if (self.IsFileReference()):
            inputFile = OrcaInputFile(OrcaInputTemplate.BASIC,
                                      calculation = calculationType,
                                      basis = self.basisSet,
                                      functional = self.functional,
                                      xyzfile = self.molecule)
        else:
            inputFile = OrcaInputFile(OrcaInputTemplate.BASICXYZ,
                                      calculation = calculationType,
                                      basis = self.basisSet,
                                      functional = self.functional,
                                      xyz = self.molecule.XYZBody())
        
        # Create the Calculation Object
        calculation = OrcaCalculation(self.name, inputFile, isLocal=self.isLocal, stdout=False)
        
        # Add a Print Statement to say we are running
        print(f"Running Frequency Analysis on {self.name}...")
        
        # Run the Calculation
        calculation.RunCalculation()
        
        # Get the Calculation Time
        self.calculationTime = time.time() - startTime
        
        # Save the Output File Path
        self.outputFilePath = calculation.OutputFilePath
        self.orcaCachePath = calculation.OrcaCachePath
        
        # Load the Output File
        outputFile = OrcaOutput(calculation.OutputFilePath)
        
        # Extract the Vibrational Frequency from the Output File
        self.vibrationalFrequencies = outputFile.get_vibrational_frequencies()
        
        # Load the IR Frequency from the 
        self.IRFrequencies = outputFile.get_ir_frequencies()
        
        # Display a Print Statement for the Frequency Completion
        print(f"Finished Running Freqency Analysis on {self.name}! ({self.ClockTime(self.calculationTime)})")
    
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