import time
import numpy as np
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
        
        
        self.IRSpectra = pd.DataFrame({
            "Wavenumber" : [],
            "IRIntensity" : []
        })
        
        
        # Create a Dictionary for the Frequency and IR Intensity Key Pairs
        IRFreqDict = { }
        
        # Add a bunch of Dummy Data
        #for i in range(0, 4000, 10):
        #    IRFreqDict[i] = 0
        
        # Convert this to a cluster calculation
        
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
            for wavenumber, intensity in zip(frequencies, IRIntensity):
                
                if wavenumber in self.IRSpectra["Wavenumber"].values:
                    index = self.IRSpectra[self.IRSpectra["Wavenumber"] == wavenumber].index
                    self.IRSpectra["IRIntensity"][index] += intensity
                else:
                    self.IRSpectra.loc[len(self.IRSpectra)] = [wavenumber, intensity]
                
                if wavenumber in IRFreqDict:
                    IRFreqDict[wavenumber] += intensity
                else:
                    IRFreqDict[wavenumber] = intensity
        
        print("\nFinished Frequency Analysis\n")
        
        print("\nMaking Final Touches\n")
        
        # Extract the Sorted Frequencies
        Frequencies = sorted(IRFreqDict.keys(), reverse=True)
        
        # Extract the Resulting Sorted IR Intensity
        IRIntensities = [IRFreqDict[freq] for freq in Frequencies]
        
        
        self.IRSpectra.sort_values(by="Wavenumber")
        
        # Add the Values to a Data Frame
       # self.IRSpectra = pd.DataFrame({
       #     "Wavenumber" : Frequencies,
       #     "IRIntensity" : IRIntensities
       # })
        
        self.IRSpectra.to_csv(f"{self.name}_Spectra.csv", index=False, header=None)
        
        print("Final Spectra")
        
        print(self.IRSpectra)
        
        print("\nFinished Making Spectra\n")
        
        # I think these are good settings
        kernelSize = int(len(IRIntensities)/8) + 1 if (int(len(IRIntensities)/8) % 2 == 0) else int(len(IRIntensities)/8)
        #print(kernelSize)
        sigma = 5 # 10
        
        kernel = self.gaussianKernel(kernelSize, sigma)

        IRIntensities = self.gaussianBlur(IRIntensities, kernel)
        
        # Normalize and Reverse the Intensity
        IRIntensities =  1 - (IRIntensities / max(IRIntensities))
        
        self.PlotSpectra()
        
        # Store the Raw DataFrame
        # Add a function to save Spectrum locally
        # Add a function to parse the spectrum
        

    def gaussianKernel(self, size, sigma):
        if size % 2 == 0:
            size -= 1
        kernel = np.exp(-np.linspace(-size // 2, size // 2, size)**2 / (2 * sigma**2))
        return kernel / kernel.sum()

    @staticmethod
    def GaussianBlur ( data, sigma):
        
        # Get Size
        size = len(data)
        
        # Make sure Size is Inpair
        if size % 2 == 0:
            size -= 1
        
        # Generate Kernel
        kernel = np.exp(-np.linspace(-size // 2, size // 2, size)**2 / (2 * sigma**2))
        kernel =  kernel / kernel.sum()
        
        # Pad the Data
        padded_data = np.pad(data, len(kernel) // 2, mode='reflect')  # Handle edges
        
        # Blurr the Data by Applying Gaussian Blurring
        blurred_data = np.convolve(padded_data, kernel, mode='valid')
        
        return blurred_data
        
    def PlotSpectra (self):
        # Plots the Spectra
        plt.figure()
        plt.plot(self.IRSpectra["Wavenumber"], self.IRSpectra["IRIntensity"])
        plt.xlabel("Wavenumber (1/cm)")
        plt.ylabel("IR Intensity")
        plt.gca().invert_xaxis()
        plt.show()
    
    @staticmethod
    def LoadSpectra (spectra : str | pd.DataFrame):
        
        # Check if the file is a DataFrame or Path
        if not isinstance(spectra, (str, pd.DataFrame)):
            raise ValueError("Spectra is not a Path to a DataFrame file or a DataFrame Object")
        
        # Determine the dataFrame Type
        if (isinstance(spectra, str)):
            # Read the Dataframe with no Header
            spectraNoHeader = pd.read_csv(spectra, header=None)
            
            # Extract the first Rows Values
            firstRow = spectraNoHeader.iloc[0].values
            
            # Make Bool Array
            isHeader = [isinstance(i, str) for i in firstRow]
            
            # Determine if the First Row is a Header or Not
            if isHeader: # First Line is the Header, we can Load it
                return pd.read_csv(spectra)
            else: # No Header, manually add them
                return pd.read_csv(spectra, names=["Wavenumber", "IRIntensity"])
            
        elif (isinstance(spectra, pd.DataFrame)):
            # Is a Dataframe, do nothing
            
            # Determine the Columns to check for
            columns = ["Wavenumber", "IRIntensity"]
            
            # Make sure all Columns are present and correctly Named
            if not all(i in spectra for i in columns):
                raise ValueError("The DataFrame provided doesn't have columns named \"Wavenumber\" and \"IRIntensity\"")
    
            return spectra.copy()
                

    @staticmethod
    def PlotSpectra (spectra: str | pd.DataFrame, sigma: int = 5, maxWaveNum = 4000):
        
        ColumnNum = 2
        
        # Load the Spectra
        IRSpectra = Spectra.LoadSpectra(spectra)
        
        # Add a bunch of Dummy Data
        for i in range(0, maxWaveNum, 1):
            IRSpectra.loc[ColumnNum] = [i, 0]
            
            #IRFreqDict[i] = 0
        
        # Sort the Columns
        IRSpectra.sort_values(by="Wavenumber", ascending=False) #Reversing?
        
        # Apply a Gaussian Blurring to the Intensities
        IRIntensity = Spectra.GaussianBlur(IRSpectra["Wavenumber"], sigma)
        
        # Normalize and Reverse Intensities
        IRIntensity = 1 - (IRIntensity / max(IRIntensity))
            
        # Add the Spectra DF Data
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        