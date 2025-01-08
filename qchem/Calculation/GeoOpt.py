import time
import os
import pandas as pd
from qchem.Molecule import Molecule
from qchem.Parser import OrcaOutput
from qchem.Calculation.OrcaInputFile import OrcaInputFile
from qchem.Calculation.OrcaCalculation import OrcaCalculation
from qchem.Data.Enums import OrcaInputTemplate, OrcaCalculationType
from .BaseOrcaCalculation import BaseOrcaCalculation
from .OrcaCalcs import RunOrcaCalculation

class GeoOpt(BaseOrcaCalculation):
    
    optMolecule: Molecule
    """The Final Optimized Molecule"""
    
    basisSet: str
    """The Basis Set to be used for the Optimization"""
    
    functional: str
    """The Density Functional to be used for the Optimization"""
    
    optimizedMoleculePath: str
    """Path to the Optimized Molecule File"""
    
    calculation: OrcaCalculation
    """Reference to the Last Orca Calculation Object for the GeoOpt"""

    def __init__(self, molecule: str | Molecule, basisSet: str, functional: str, index: int = 1, cores:int = 1, isLocal:bool = False, name:str = "OPTMolecule"):
        
        # Make a Super Call (Use the Base Class Init for some Boilerplate Setup)
        super().__init__(name, molecule, index, cores, isLocal, True)
        
        # Check if Values are empty or of Wrong Type
        if (not (basisSet and isinstance(basisSet, (str)))):
            raise ValueError("BasisSet is not defined! Provide the Name of the Basis Set as a String")
        
        if (not functional and isinstance(functional, (str))):
            raise ValueError("Functional is not defined! Provide the Name of the Functional as a String")
            
        # Set the Values
        self.basisSet = basisSet
        self.functional = functional

    def RunCalculation (self):
        """Runs through a Geometry Optimization on the Molecule and repeats until properly converged"""
        # Get Default Values
        OPTtemplate = ""
        xyzMol = ""
        optIndex = 1
        isOptimized = False
        freqFailCount = 0
        calculationType = f"{OrcaCalculationType.OPTIMIZATION.value} {OrcaCalculationType.FREQUENCY.value}"
        
        startTime = time.time()
        
        # Determine the Proper Template to use based on the Molecule Type
        if (self.IsFileReference()):
            OPTtemplate = OrcaInputTemplate.BASICPARALLEL
            xyzMol = self.molecule
            
            inputFile = OrcaInputFile(OPTtemplate,
                                    calculation = calculationType,
                                    basis=self.basisSet,
                                    functional=self.functional,
                                    cores = self.cores,
                                    xyzfile=xyzMol)
        else:
            OPTtemplate = OrcaInputTemplate.BASICXYZPARALLEL
            xyzMol = self.molecule.XYZBody()
            
            inputFile = OrcaInputFile(OPTtemplate, 
                                    calculation = calculationType,
                                    basis=self.basisSet,
                                    functional=self.functional,
                                    cores = self.cores,
                                    xyz=xyzMol)
        
        # Start the Optimization Loop
        while (not isOptimized):
            
            iterStartTime = time.time()
            
            # Generate Print Statement for User on the Optimization Attempt
            print(f"Running Optimization Attempt {optIndex} on {self.name}...")
            
            # Generate an Indexed Name for each Calculation
            calcName = self.name if (optIndex == 1) else self.name + f"_{optIndex}"
            
            # Run the Calculation
            calculation = RunOrcaCalculation(calcName, inputFile, isLocal=self.isLocal, STDOut=False)
            
            # Rune the Calculation
            #calculation.RunCalculation()
            
            # Get the Output File
            outputFile = OrcaOutput(calculation.outputFilePath)
            
            # Check if Frequencies are Empty and cannot be found
            if outputFile.vibrational_frequencies is None or (isinstance(outputFile.vibrational_frequencies, pd.DataFrame) and outputFile.vibrational_frequencies.empty):
                print("No Frequencies Found! Optimization Failed!")
                freqFailCount += 1
                if (freqFailCount >= 3):
                    print("Failed to Optimize Molecule after 3 Attempts! Aborting Optimization!")
                    return
            else:
                # Check if the Molecule is Fully Optimized
                if (self.IsOptimized(outputFile.vibrational_frequencies["frequency"])):
                    self.calculationTime = time.time() - startTime
                    print(f"Molecule {self.name} is Optimized! ({self.ClockTime(self.calculationTime)})")
                    isOptimized = True
                    self.calculation = calculation
                    self.optimizedMoleculePath = os.path.join(calculation.orcaCachePath,calculation.name + ".xyz")
                    self.optMolecule = Molecule(self.name, self.optimizedMoleculePath)
                    break
                
            calcTime = time.time() - iterStartTime
            print(f"Finished Optimization Attempt {optIndex} on {self.name} ({self.ClockTime(calcTime)})")
            
            # Update the Molecule and Optimization Template for the Next Iteration
            self.optimizedMoleculePath = os.path.join(calculation.orcaCachePath, calculation.name + ".xyz")
            OPTtemplate = OrcaInputTemplate.BASICXYZPARALLEL
            xyzMol = Molecule(self.name, self.optimizedMoleculePath).XYZBody()
            optIndex += 1
            
            # Generate the Input File
            inputFile = OrcaInputFile(OPTtemplate, 
                                    calculation = calculationType,
                                    basis=self.basisSet,
                                    functional=self.functional,
                                    cores = self.cores,
                                    xyz=xyzMol)
            
            # Create the Calculation Object
            #calculation = OrcaCalculation(self.name + f"_{optIndex}", inputFile, isLocal=self.isLocal, stdout=False)
            
    def IsOptimized(self, frequencies):
        """Determines if the Molecule is Optimized to the Valid Minimum"""
        for freq in frequencies:
            if (freq < 0):
                return False
        
        return True
    