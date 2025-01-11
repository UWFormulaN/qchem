import os
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from qchem.Molecule import Molecule
from qchem.Calculation.GOAT import GOAT
from qchem.Calculation.GeoOpt import GeoOpt
from qchem.Data.Enums import OrcaInputTemplate
from qchem.Calculation.Frequency import Frequency
from qchem.Calculation.BaseOrcaCalculation import BaseOrcaCalculation

class Spectra(BaseOrcaCalculation):

    calculationType: str = "Spectra"
    
    IRSpectra: pd.DataFrame
    """The Resulting IR Spectra plotted"""

    def __init__(
        self,
        molecule: str | Molecule,
        template: str | OrcaInputTemplate = "",
        index: int = 1,
        cores: int = 1,
        isLocal: bool = False,
        name: str = "FREQMolecule",
        stdout: bool = True,
        **variables
    ):
        # Make a Super Call (Use the Base Class Init for some Boilerplate Setup)
        super().__init__(name, molecule, template, index, cores, isLocal, stdout, **variables)

        # Check if the Calculation has a Basis Set and a Functional Defined (Specific to Certain Calculations)
        self.BasisSetFunctionalCompliant()
    
        # Delete Cores from Variables cause it causes Issues        
        self.variables.pop("cores")

    def RunCalculation(self):
        """Runs through all Calculations required to produce a Spectra Graph"""

        # Start the Timer
        startTime = time.time()

        # Make Cache Folder if it doesn't Exist
        if not os.path.exists(self.orcaCachePath):
            os.makedirs(self.orcaCachePath)

        print("\nRunning GeoOpt!\n")

        # Create a Geo Opt Calculation Object
        geoOptCalc = GeoOpt(self.molecule, True, self.template, self.index, self.cores, self.isLocal, f"{self.name}_GEOOPT", False, **self.variables)
    
        # Run the GeoOptimization on the Molecule
        geoOptCalc.RunCalculation()

        print("\nFinished GeoOpt!\n")
        print("\nRunning GOAT!\n")

        # Create GOAT Calculation Object
        goatCalc = GOAT(geoOptCalc.optMolecule, self.template, self.index, self.cores, self.isLocal, f"{self.name}_GOAT", False, **self.variables)
        
        # Run the GOAT Calculation
        goatCalc.RunCalculation()

        print("\nFinished GOAT!\n")
        print("\nRunning Frequency Analysis!\n")

        # Get the Number of Conformers Created
        conformersNum = len(goatCalc.conformers)

        # Create Blank DataFrame
        self.IRSpectra = pd.DataFrame({"Wavenumber": [], "IRIntensity": []})

        # Load the IR Contributions
        IRContribution = pd.DataFrame(
            {
                "Name": [f"{self.name} {i}" for i in range(conformersNum)],
                "Contribution": goatCalc.conformerContribution[:conformersNum],
            }
        )

        # Save the Contributions to Excel File
        IRContribution.to_csv(
            os.path.join(self.orcaCachePath, f"{self.name}_Contributions.csv"), index=False
        )

        # Loop through all the Conformers and Run a Frequency Calculation
        for i in range(conformersNum):

            # Create the Frequency Calculation
            freqCalc = Frequency(goatCalc.conformers[i], self.template, self.index, self.cores, self.isLocal, f"{self.name}_FREQ_{i}", False, **self.variables)
            
            # Run the Frequency Calculation
            freqCalc.RunCalculation()

            # Load Individual Spectra Results
            FreqSpectra = pd.DataFrame(
                {
                    "Wavenumber": freqCalc.IRFrequencies["frequency"].values,
                    "IRIntensity": freqCalc.IRFrequencies["IR_intensity"].values,
                }
            )

            # Save Individual Spectra
            FreqSpectra.to_csv(
                os.path.join(self.orcaCachePath, f"{self.name}_IRIntensity_{i}.csv"),
                index=False,
            )

            # Multiply Spectra by Contribution
            FreqSpectra["IRIntensity"] = (
                FreqSpectra["IRIntensity"].values * goatCalc.conformerContribution[i]
            )

            # This Method Works, no need for a Dictionary and all that
            self.IRSpectra = pd.concat([self.IRSpectra, FreqSpectra], ignore_index=True)

        print("\nFinished Frequency Analysis!\n")
        print("\nMaking Final Touches\n")

        # Group Common Wavenumbers and Sum their Values
        self.IRSpectra.groupby("Wavenumber", as_index=False).agg({"IRIntensity": "sum"})

        # Sort Wavenumbers from Largest Wavenumber -> Lowest Wavenumber
        self.IRSpectra = self.IRSpectra.sort_values(by="Wavenumber", ascending=False)

        # Save Full Spectra
        self.IRSpectra.to_csv(
            os.path.join(self.orcaCachePath, f"{self.name}_Spectra.csv"), index=False
        )

        # Get Total Time for Spectra
        calcTime = time.time() - startTime

        print(f"\nFinished Making {self.name} Spectra! ({self.ClockTime(calcTime)})\n")

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
        padded_data = np.pad(data, len(kernel) // 2, mode="reflect")

        # Blurr the Data by Applying Gaussian Blurring
        blurred_data = np.convolve(padded_data, kernel, mode="valid")

        return blurred_data

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
            if all(isHeader):  # First Line is the Header, we can Load it
                return pd.read_csv(spectra)
            else:  # No Header, manually add them
                return pd.read_csv(spectra, names=["Wavenumber", "IRIntensity"])

        elif isinstance(spectra, pd.DataFrame):
            
            # Determine the Columns to check for
            columns = ["Wavenumber", "IRIntensity"]

            # Make sure all Columns are present and correctly Named
            if not all(i in spectra for i in columns):
                raise ValueError(
                    'The DataFrame provided doesn\'t have columns named "Wavenumber" and "IRIntensity"'
                )

            # Return a Copy so that Original isn't destroyed
            return spectra.copy()

    @staticmethod
    def PlotSpectra(
        spectra: str | pd.DataFrame,
        plotName: str = "Spectra",
        sigma: int = 5,
        maxWaveNum=4000,
        spacing=10,
        showPlot = True
    ):
        """Create and Plots an IR Spectra, takes a path to a CSV or Pandas Dataframe with the Calculated"""

        # Load the Spectra
        IRSpectra = Spectra.LoadSpectra(spectra)

        # Add a bunch of Dummy Data
        for i in range(0, maxWaveNum, spacing):
            IRSpectra.loc[len(IRSpectra)] = [i, 0]

        # Group Duplicate Wavenumbers if they Appear
        IRSpectra.groupby("Wavenumber", as_index=False).agg({"IRIntensity": "sum"})

        # Sort the Columns from Largest Wavenumber -> Smallest Wavenumber
        IRSpectra = IRSpectra.sort_values(by="Wavenumber", ascending=False)

        # Apply a Gaussian Blurring to the Intensities
        IRSpectra["IRIntensity"] = Spectra.GaussianBlur(IRSpectra["IRIntensity"], sigma)

        # Normalize and Inverse Intensity
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
        
        if showPlot:
            plt.show()