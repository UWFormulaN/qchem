import time
import pandas as pd
from qchem.Molecule import Molecule
from qchem.OrcaInputFile import OrcaInputFile
from qchem.Data.Enums import OrcaInputTemplate, OrcaCalculationType
from qchem.OrcaCalculation import OrcaCalculation
from qchem.Parser import OrcaOutput

class GeoOpt:
    
    molecule: Molecule | str
    """The Molecule to be Optimized"""
    
    optMolecule: Molecule
    """The Final Optimized Molecule"""
    
    basisSet: str
    """The Basis Set to be used for the Optimization"""
    
    functional: str
    """The Density Functional to be used for the Optimization"""
    
    isLocal : bool
    """Boolean Flag to determine if the Optimization is Local or not (Local = Runs on Device, Non-Local = Runs in Docker Container)"""
    
    name: str
    """Name of the Molecule being GeoOptimized"""
    
    optimizationTime : float
    """Time Elapsed for the Optimization to Complete"""
    
    optimizedMoleculePath: str
    """Path to the Optimized Molecule File"""
    
    calculation: OrcaCalculation
    """Reference to the Last Orca Calculation Object for the GeoOpt"""
    
    def __init__(self, molecule, basisSet, functional, cores:int = 1, isLocal:bool = False, name:str = ""):
        
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
                print("No Name Provided for the Molecule, using Default Name: GeoOptMolecule")
                self.name = "GeoOptMolecule"
        else:
            self.name = name
            
        # Set the Values
        self.molecule = molecule
        self.basisSet = basisSet
        self.functional = functional
        self.cores = cores
        self.isLocal = isLocal
        
    def IsFileReference(self):
        """Determines if the Molecule is stored as a File Reference, or is a direct Molecule Object"""
        if (isinstance(self.molecule, (str))):
            return True
        else:
            return False
        
    def Optimize (self):
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
        
        # Create the Calculation Object
        calculation = OrcaCalculation(self.name, inputFile, isLocal=self.isLocal, stdout=False)    
        
        # Start the Optimization Loop
        while (not isOptimized):
            
            iterStartTime = time.time()
            
            # Generate Print Statement for User on the Optimization Attempt
            print(f"Running Optimization Attempt {optIndex} on {self.name}")
            
            # Rune the Calculation
            calculation.RunCalculation()
            
            # Get the Output File
            outputFile = OrcaOutput(calculation.OutputFilePath)
            
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
                    calcTime = time.time() - startTime
                    print(f"Molecule {self.name} is Optimized! ({self.ClockTime(calcTime)})")
                    isOptimized = True
                    self.calculation = calculation
                    break
                
            calcTime = time.time() - iterStartTime
            print(f"Finished Optimization Attempt {optIndex} on {self.name} ({self.ClockTime(calcTime)})")
            
            # Update the Molecule and Optimization Template for the Next Iteration
            self.optimizedMoleculePath = calculation.OrcaCachePath + f"\\{calculation.CalculationName}.xyz"
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
            calculation = OrcaCalculation(self.name + f"_{optIndex}", inputFile, isLocal=self.isLocal, stdout=False)
            
    def IsOptimized(self, frequencies):
        """Determines if the Molecule is Optimized to the Valid Minimum"""
        for freq in frequencies:
            if (freq < 0):
                return False
        
        return True
    
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
        