import os
import time
import pandas as pd
from qchem.Molecule import Molecule
from qchem.Parser import OrcaOutput
from .BaseOrcaCalculation import BaseOrcaCalculation
from qchem.Calculation.OrcaInputFile import OrcaInputFile
from qchem.Data.Enums import OrcaInputTemplate, OrcaCalculationType
from qchem.Calculation.OrcaCalculation import runOrcaCalculation, OrcaCalcResult


class GeoOpt(BaseOrcaCalculation):
    """Performs a Frequency Calculation on a Molecule. Can be set to Optimize until convergence or only run a single Calculation"""

    #
    # Need to be Set
    #
    calculationType: str = (
        f"{OrcaCalculationType.OPTIMIZATION.value} {OrcaCalculationType.FREQUENCY.value}"
    )
    """The Keyword for the Calculation to run on the Molecule (For Pipelines replace with name)"""

    fullOptimization: bool
    """A Boolean Flag determining if the Molecule will be Completely Optimized or will run a Single Optimization Iteration (True = Runs until Convergence, False = Runs for a Single Iteration)"""

    optMolecule: Molecule
    """The Resulting Optimized Molecule once the Calculation is Complete"""

    optimizedMoleculePath: str
    """The Path to the XYZ File of the Resulting Optimized Molecule"""

    calculation: OrcaCalcResult
    """The reference to the latest Calculation Results from the GeoOpt """

    def __init__(
        self,
        molecule: str | Molecule,
        fullOptimization: bool = True,
        template: str | OrcaInputTemplate = "",
        index: int = 1,
        cores: int = 1,
        isLocal: bool = False,
        name: str = "Molecule",
        stdout: bool = True,
        **variables,
    ):

        # Make a Super Call (Use the Base Class Init for some Boilerplate Setup)
        super().__init__(
            name, molecule, template, index, cores, isLocal, stdout, **variables
        )

        # Check if the Calculation has a Basis Set and a Functional Defined (Specific to Certain Calculations)
        self.basisSetFunctionalCompliant()

        self.fullOptimization = fullOptimization

    def runCalculation(self):
        """Runs the GeoOpt Calculation and Saves the Optimized Molecule and it's Path
        
        Parameters: \n
            self - Default Parameter for the Class Instance
            
        Returns : \n
            None - No Return Value
        """

        # Single Optimization
        if not self.fullOptimization:
            self.singleOptimization()
        else:
            self.completeOptimization()

    def singleOptimization(self):
        """Runs a Single Optimization Attempt on the Molecule. If the Molecule has not Converged a warning will be sent to the Terminal
        
        Parameters : \n
            self - Default Parameter for the Class Instance
            
        Returns : \n
            None - No Return Value
        """

        print(f"Running OPT on {self.name}...")

        # Start the Timer
        startTime = time.time()

        # Run the Calculation
        calculation = runOrcaCalculation(
            self.name, self.inputFile, isLocal=self.isLocal, STDOut=False
        )

        # Get the Output File
        outputFile = OrcaOutput(calculation.outputFilePath)

        # Check if the Vibrational Frequencies Exist and/or is not empty
        if outputFile.vibrationalFrequencies is None or (
            isinstance(outputFile.vibrationalFrequencies, pd.DataFrame)
            and outputFile.vibrationalFrequencies.empty
        ):
            print("No Frequencies Found! Optimization Failed!")

        else:
            # Check if the Molecule is Fully Optimized
            if self.isOptimized(outputFile.vibrationalFrequencies["frequency"]):
                self.calculationTime = time.time() - startTime
                print(
                    f"Molecule {self.name} is Optimized! ({self.clockTime(self.calculationTime)})"
                )
                self.calculation = calculation
                self.optimizedMoleculePath = os.path.join(
                    calculation.orcaCachePath, calculation.name + ".xyz"
                )
                self.optMolecule = Molecule(self.name, self.optimizedMoleculePath)
            else:
                print(f"Molecule {self.name} is not Optimized!")

        # Stop Timer and Get Total time in Seconds
        calcTime = time.time() - startTime

        # Print a Finishing Statement with the Time
        print(f"Finished OPT on {self.name} ({self.clockTime(calcTime)})")

    def completeOptimization(self):
        """Runs an Algorithm that will continuously run Geometry Optimizations until the Molecule has been completely Optimized and Converged
        
        Parameters : \n
            self - Default Parameter for the Class Instance
            
        Returns : \n
            None - No Return Value"""
        # Start the Timer
        startTime = time.time()

        # Set the Optimization Index, and Initialize the Optimization Flag
        optIndex = 1
        isOptimized = False
        freqFailCount = 0

        # Full Optimization Loop
        while not isOptimized:

            # Start the Individual Iteration Timer
            iterStartTime = time.time()

            # Generate Print Statement for User on the Optimization Attempt
            print(f"Running OPT {optIndex} on {self.name}...")

            # Generate an Indexed Name for each Calculation
            calcName = self.name if (optIndex == 1) else self.name + f"_{optIndex}"

            # Run the Calculation
            calculation = runOrcaCalculation(
                calcName, self.inputFile, isLocal=self.isLocal, STDOut=False
            )

            # Get the Output File
            outputFile = OrcaOutput(calculation.outputFilePath)

            # Check if the Vibrational Frequencies Exist and/or is not empty
            if outputFile.vibrationalFrequencies is None or (
                isinstance(outputFile.vibrationalFrequencies, pd.DataFrame)
                and outputFile.vibrationalFrequencies.empty
            ):
                print("No Frequencies Found! Optimization Failed!")
                freqFailCount += 1
                if freqFailCount >= 3:
                    print(
                        "Failed to Optimize Molecule after 3 Attempts! Aborting Optimization!"
                    )
                    return
            else:
                # Check if the Molecule is Fully Optimized
                if self.isOptimized(outputFile.vibrationalFrequencies["frequency"]):
                    self.calculationTime = time.time() - startTime
                    print(
                        f"Molecule {self.name} is Optimized! ({self.clockTime(self.calculationTime)})"
                    )
                    isOptimized = True
                    self.calculation = calculation
                    self.optimizedMoleculePath = os.path.join(
                        calculation.orcaCachePath, calculation.name + ".xyz"
                    )
                    self.optMolecule = Molecule(self.name, self.optimizedMoleculePath)
                    break

            calcTime = time.time() - iterStartTime
            print(
                f"Finished OPT {optIndex} on {self.name} ({self.clockTime(calcTime)})"
            )

            # Update the Molecule and Optimization Template for the Next Iteration
            self.optimizedMoleculePath = os.path.join(
                calculation.orcaCachePath, calculation.name + ".xyz"
            )

            # Update the Molecule and Optimization Template for the Next Iteration
            self.template = OrcaInputTemplate.BASICXYZPARALLEL
            self.variables["xyz"] = Molecule(
                self.name, self.optimizedMoleculePath
            ).XYZBody()
            optIndex += 1

            # Generate the Input File
            self.inputFile = OrcaInputFile(self.template, **self.variables)

    def isOptimized(self, frequencies: list[float]):
        """Checks if the Molecule has been Optimized using the Vibrational Frequencies. Returns a boolean indicating if Optimized (Fully Optimized = All Frequencies > 0)
        
        Parameters : \n
            self - Default Parameter for the Class Instance \n
            frequencies - A list of Vibrational Frequencies
        
        Returns : \n
            bool - True if Molecule is Optimized, False if not Optimized
        """
        for freq in frequencies:
            if freq < 0:
                return False

        return True
