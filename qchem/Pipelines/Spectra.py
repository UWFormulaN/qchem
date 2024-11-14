import time
import pandas as pd
import matplotlib.pyplot as plt
from qchem.Molecule import Molecule
from qchem.Parser import OrcaOutput
from qchem.Calculation.GeoOpt import GeoOpt
from qchem.Calculation.GOAT import GOAT
from qchem.Calculation.Frequency import Frequency
from qchem.Calculation.OrcaInputFile import OrcaInputFile
from qchem.Calculation.OrcaCalculation import OrcaCalculation
from qchem.Data.Enums import OrcaInputTemplate

class Spectra:
    
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
    
    IRSpectra : pd.DataFrame
    """The Resulting IR Spectra plotted"""
    
    def __init__ (self, molecule: str | Molecule, basisSet: str, functional: str, cores:int = 1, isLocal:bool = False, name:str = ""):
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
    
    def IsFileReference(self):
        """Determines if the Molecule is stored as a File Reference, or is a direct Molecule Object"""
        if (isinstance(self.molecule, (str))):
            return True
        else:
            return False
        
    def RunCalculation (self):
        """Runs through all Calculations required to produce a Spectra Graph"""
    
        print("\nRunning GeoOpt\n")
    
        # Create a Geo Opt Calculation Object
        geoOptCalc = GeoOpt(self.molecule, self.basisSet, self.functional, self.cores, self.isLocal, f"{self.name}_GEOOPT")
    
        # Run the GeoOptimization on the Molecule
        geoOptCalc.RunCalculation()
        
        print("\nFinished GeoOpt\n")
        
        print("\nRunning GOAT\n")
        
        # Create GOAT Calculation Object 
        goatCalc = GOAT(geoOptCalc.optMolecule, self.cores, self.isLocal, f"{self.name}_GOAT")
        
        # Run the GOAT Calculation
        goatCalc.RunCalculation()
        
        print("\nFinished GeoOpt\n")
        
        print("\nRunning Frequency Analysis\n")
        
        # Get the Number of Conformers Created
        conformersNum = len(goatCalc.conformers)
        
        # Create a Dictionary for the Frequency and IR Intensity Key Pairs
        IRFreqDict = { }
        
        # Loop through all the Conformers and Run a Frequency Calculation
        for i in range(conformersNum):
            
            # Create the Frequency Calculation
            freqCalc = Frequency(self.molecule, self.basisSet, self.functional, self.cores, self.isLocal, f"{self.name}_FREQ_{i}")
            
            # Run the Frequency Calculation
            freqCalc.RunCalculation()
            
            # Extract the Frequency Values
            frequencies = freqCalc.IRFrequencies["frequency"].values 
            
            # Extract the IR Intensity Values and Multiply by the conformer Contribution
            IRIntensity = freqCalc.IRFrequencies["IR_intensity"].values * goatCalc.conformerContribution[i]
            
            # Add the Values to the Dictionary
            for frequency, intensity in zip(frequencies, IRIntensity):
                if frequency in IRFreqDict:
                    IRFreqDict[frequency] += intensity
                else:
                    IRFreqDict[frequency] = intensity
        
        print("\nFinished Frequency Analysis\n")
        
        print("\nMaking Final Touches\n")
        
        # Extract the Sorted Frequencies
        Frequencies = sorted(IRFreqDict.keys(), reverse=True)
        
        # Extract the Resulting Sorted IR Intensity
        IRIntensities = [IRFreqDict[freq] for freq in Frequencies]
        
        # Normalize and Reverse the Intensity
        IRIntensities =  1 - (IRIntensities / max(IRIntensities))
        
        # Add the Values to a Data Frame
        self.IRSpectra = pd.DataFrame({
            "Frequency" : Frequencies,
            "IRIntensity" : IRIntensities
        })
        
        print("\nFinished Making Spectra\n")

    def PlotSpectra (self):
        # Plots the Spectra
        plt.figure()
        plt.plot(self.IRSpectra["Frequency"], self.IRSpectra["IRIntensity"])
        plt.xlabel("Frequency (1/cm)")
        plt.ylabel("IR Intensity")
        plt.show()
        
        
        
        
        
        