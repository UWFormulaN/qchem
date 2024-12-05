import time
import pandas as pd
from qchem.Molecule import Molecule
from qchem.Parser import OrcaOutput
from qchem.Calculation.OrcaInputFile import OrcaInputFile
from qchem.Calculation.OrcaCalculation import OrcaCalculation
from qchem.Data.Enums import OrcaInputTemplate
from .BaseOrcaCalculation import BaseOrcaCalculation
from .OrcaCalcs import RunOrcaCalculation

class Frequency(BaseOrcaCalculation):
    
    # Most of this is Boilerplate, Should make a Abstract class to inherit this boilerplate logic
    defaultName: str = "FREQMolecule"
    
    basisSet: str
    """The Basis Set to be used for the Optimization"""
    
    functional: str
    """The Density Functional to be used for the Optimization"""
    
    vibrationalFrequencies : pd.DataFrame
    """Vibrational Frequencies of the Molecule"""
    
    IRFrequencies : pd.DataFrame
    """IR Frequencies of the Molecule"""
    
    def __init__(self, molecule: str | Molecule, basisSet: str, functional: str, index: int = 1, cores:int = 1, isLocal:bool = False, name:str = defaultName):
        
        # Make a Super Call (Use the Base Class Init for some Boilerplate Setup)
        super().__init__(name, molecule, index, cores, isLocal, True)
        
        if (not (basisSet and isinstance(basisSet, (str)))):
            raise ValueError("BasisSet is not defined! Provide the Name of the Basis Set as a String")
        
        if (not functional and isinstance(functional, (str))):
            raise ValueError("Functional is not defined! Provide the Name of the Functional as a String")
        
        if (isinstance(molecule, Molecule) and name == self.defaultName):
            self.name = molecule.Name
        
        # Set the Values
        self.basisSet = basisSet
        self.functional = functional
        self.conformers = []
        
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
        
        # Add a Print Statement to say we are running
        print(f"Running Frequency Analysis on {self.name}...")
        
        # Run the Orca Calculation
        calculation = RunOrcaCalculation(self.name, inputFile, self.index, self.isLocal, False)
        
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
        print(f"Finished Running Freqency Analysis on {self.name}! ({self.ClockTime(self.calculationTime)})")