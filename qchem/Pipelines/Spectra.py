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

    isLocal: bool
    """Boolean Flag to determine if the Optimization is Local or not (Local = Runs on Device, Non-Local = Runs in Docker Container)"""

    name: str
    """Name of the Molecule being GeoOptimized"""

    calculationTime: float
    """Time Elapsed for the GOAT Optimization to Complete"""

    orcaCachePath: str
    """The Path to the Orca Cache folder for the Calculation"""

    outputFilePath: str
    """The Path to the Output File"""

    basisSet: str
    """The Basis Set to be used for the Optimization"""

    IRSpectra: pd.DataFrame
    """The Resulting IR Spectra plotted"""

    def __init__(
        self,
        molecule: str | Molecule,
        basisSet: str,
        functional: str,
        cores: int = 1,
        isLocal: bool = False,
        name: str = "",
    ):
        # Check if Values are empty or of Wrong Type
        if not (molecule and isinstance(molecule, (str, Molecule))):
            raise ValueError(
                "Molecule is not defined! Provide a Path to the XYZ file or a Molecule Object"
            )

        if not (basisSet and isinstance(basisSet, (str))):
            raise ValueError(
                "BasisSet is not defined! Provide the Name of the Basis Set as a String"
            )

        if not functional and isinstance(functional, (str)):
            raise ValueError(
                "Functional is not defined! Provide the Name of the Functional as a String"
            )

        if not isinstance(cores, (int)):
            raise ValueError("Cores is not an Integer! Provide an Integer Value")

        if not isinstance(isLocal, (bool)):
            raise ValueError("IsLocal is not a Boolean! Provide a Boolean Value")

        if not isinstance(name, (str)):
            raise ValueError("Name is not a String! Provide a String Value")

        # Set the Name of the Molecule
        if name == "":
            if isinstance(molecule, (Molecule)):
                self.name = molecule.Name
            else:
                print(
                    "No Name Provided for the Molecule, using Default Name: GOATMolecule"
                )
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
        if isinstance(self.molecule, (str)):
            return True
        else:
            return False

    def RunCalculation(self):
        """Runs through all Calculations required to produce a Spectra Graph"""

        print("\nRunning GeoOpt\n")

        # Create a Geo Opt Calculation Object
        geoOptCalc = GeoOpt(
            self.molecule,
            self.basisSet,
            self.functional,
            self.cores,
            self.isLocal,
            f"{self.name}_GEOOPT",
        )

        # Run the GeoOptimization on the Molecule
        geoOptCalc.RunCalculation()

        print("\nFinished GeoOpt\n")

        print("\nRunning GOAT\n")

        # Create GOAT Calculation Object
        goatCalc = GOAT(
            geoOptCalc.optMolecule, self.cores, self.isLocal, f"{self.name}_GOAT"
        )

        # Run the GOAT Calculation
        goatCalc.RunCalculation()

        print("\nFinished GeoOpt\n")

        print("\nRunning Frequency Analysis\n")

        # Get the Number of Conformers Created
        conformersNum = len(goatCalc.conformers)

        self.IRSpectra = pd.DataFrame({"Wavenumber": [], "IRIntensity": []})

        # IRContribution = pd.DataFrame({
        #    "Name" : [],
        #    "Contribution" : []
        # })
        #
        # for i in range(conformersNum):
        #    IRContribution.loc[len(IRContribution)] = [f"{self.name} {i}", goatCalc.conformerContribution[i]]

        # Load the IR Contributions
        IRContribution = pd.DataFrame(
            {
                "Name": [f"{self.name} {i}" for i in range(conformersNum)],
                "Contribution": goatCalc.conformerContribution[:conformersNum],
            }
        )

        # Save the Contributions to Excel File
        IRContribution.to_csv(f"{self.name}_Contributions.csv", index=False)

        # Convert this to a cluster calculation

        # Loop through all the Conformers and Run a Frequency Calculation
        for i in range(conformersNum):

            # Create the Frequency Calculation
            freqCalc = Frequency(
                goatCalc.conformers[i],
                self.basisSet,
                self.functional,
                self.cores,
                self.isLocal,
                f"{self.name}_FREQ_{i}",
            )

            # Run the Frequency Calculation
            freqCalc.RunCalculation()

            FreqSpectra = pd.DataFrame(
                {
                    "Wavenumber": freqCalc.IRFrequencies["frequency"].values,
                    "IRIntensity": freqCalc.IRFrequencies["IR_intensity"].values,
                }
            )

            FreqSpectra.to_csv(
                f"{self.name}_IRIntensity_{i}.csv",
                index=False,
            )

            FreqSpectra["IRIntensity"] = (
                FreqSpectra["IRIntensity"].values * goatCalc.conformerContribution[i]
            )

            # This Method Works, no need for a Dictionary and all that
            self.IRSpectra = pd.concat([self.IRSpectra, FreqSpectra], ignore_index=True)

        print("\nFinished Frequency Analysis\n")

        print("\nMaking Final Touches\n")

        self.IRSpectra.groupby("Wavenumber", as_index=False).agg({"IRIntensity": "sum"})

        self.IRSpectra = self.IRSpectra.sort_values(by="Wavenumber", ascending=False)

        self.IRSpectra.to_csv(f"{self.name}_Spectra.csv", index=False)

        print("Final Spectra")

        print(self.IRSpectra)

        print("\nFinished Making Spectra\n")

        print("Plotting")
        Spectra.PlotSpectra(self.IRSpectra, "First")

        # Spectra.PlotSpectra(IRSpec, "SpecTest")

        # Store the Raw DataFrame
        # Add a function to save Spectrum locally
        # Add a function to parse the spectrum

        # Ok so Old Method is Working Properly

    def gaussianKernel(self, size, sigma):
        if size % 2 == 0:
            size -= 1
        kernel = np.exp(-np.linspace(-size // 2, size // 2, size) ** 2 / (2 * sigma**2))
        return kernel / kernel.sum()

    @staticmethod
    def GaussianBlur(data, sigma):

        # Get Size
        size = int(len(data))

        # Make sure Size is Inpair
        if size % 2 == 0:
            size -= 1

        # Generate Kernel
        kernel = np.exp(-np.linspace(-size // 2, size // 2, size) ** 2 / (2 * sigma**2))
        kernel = kernel / kernel.sum()

        # Pad the Data
        padded_data = np.pad(data, len(kernel) // 2, mode="reflect")  # Handle edges

        # Blurr the Data by Applying Gaussian Blurring
        blurred_data = np.convolve(padded_data, kernel, mode="valid")

        return blurred_data

    def PlotSpectra2(self):
        # Plots the Spectra
        plt.figure()
        plt.plot(self.IRSpectra["Wavenumber"], self.IRSpectra["IRIntensity"])
        plt.xlabel("Wavenumber (1/cm)")
        plt.ylabel("IR Intensity")
        plt.gca().invert_xaxis()
        plt.show()

    @staticmethod
    def LoadSpectra(spectra: str | pd.DataFrame):

        # Check if the file is a DataFrame or Path
        if not isinstance(spectra, (str, pd.DataFrame)):
            raise ValueError(
                "Spectra is not a Path to a DataFrame file or a DataFrame Object"
            )

        # Determine the dataFrame Type
        if isinstance(spectra, str):
            # Read the Dataframe with no Header
            spectraNoHeader = pd.read_csv(spectra, header=None)

            # Extract the first Rows Values
            firstRow = spectraNoHeader.iloc[0].values

            # Make Bool Array
            isHeader = [isinstance(i, str) for i in firstRow]

            # Determine if the First Row is a Header or Not
            if isHeader:  # First Line is the Header, we can Load it
                return pd.read_csv(spectra)
            else:  # No Header, manually add them
                return pd.read_csv(spectra, names=["Wavenumber", "IRIntensity"])

        elif isinstance(spectra, pd.DataFrame):
            # Is a Dataframe, do nothing

            # Determine the Columns to check for
            columns = ["Wavenumber", "IRIntensity"]

            # Make sure all Columns are present and correctly Named
            if not all(i in spectra for i in columns):
                raise ValueError(
                    'The DataFrame provided doesn\'t have columns named "Wavenumber" and "IRIntensity"'
                )

            print("Exists")

            return spectra.copy()

    @staticmethod
    def PlotSpectra(
        spectra: str | pd.DataFrame,
        plotName: str = "Spectra",
        sigma: int = 5,
        maxWaveNum=4000,
    ):

        # Load the Spectra
        IRSpectra = Spectra.LoadSpectra(spectra)

        # Add a bunch of Dummy Data
        for i in range(0, maxWaveNum, 10):
            IRSpectra.loc[len(IRSpectra)] = [i, 0]

        print("Printing Spectra")
        print(IRSpectra)

        # Sort the Columns
        IRSpectra = IRSpectra.sort_values(
            by="Wavenumber", ascending=False
        )  # Reversing?

        # Apply a Gaussian Blurring to the Intensities
        IRSpectra["IRIntensity"] = Spectra.GaussianBlur(IRSpectra["IRIntensity"], sigma)

        # Normalize and Reverse Intensities
        IRSpectra["IRIntensity"] = 1 - (
            IRSpectra["IRIntensity"].values / max(IRSpectra["IRIntensity"].values)
        )

        # Plots the Spectra
        plt.figure()
        plt.plot(IRSpectra["Wavenumber"], IRSpectra["IRIntensity"])
        plt.xlabel("Wavenumber (1/cm)")
        plt.ylabel("IR Intensity")
        plt.gca().invert_xaxis()
        plt.savefig(f"{plotName}.png")
        plt.show()
