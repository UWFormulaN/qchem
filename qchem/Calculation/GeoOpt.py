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

    #
    # Need to be Set
    #
    calculationType: str = (
        f"{OrcaCalculationType.OPTIMIZATION.value} {OrcaCalculationType.FREQUENCY.value}"
    )

    fullOptimization: bool
    """Boolean Flag to determine if the Optimization is Full or not (True = Runs until Convergence, False = Runs for a Single Iteration)"""

    optMolecule: Molecule
    """The Final Optimized Molecule"""

    optimizedMoleculePath: str
    """Path to the Optimized Molecule File"""

    calculation: OrcaCalcResult
    """Reference to the Last Orca Calculation Object for the GeoOpt"""

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
        """Runs through a Geometry Optimization on the Molecule and repeats until properly converged"""

        # Single Optimization
        if not self.fullOptimization:
            self.singleOptimization()
        else:
            self.fullOptimization()

    def singleOptimization(self):
        """Runs a Single Optimization Attempt on the Molecule (Will Output a message if the Molecule is not fully optimized)"""

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

    def fullOptimization(self):
        """Runs a Full Geo Optimization on the Molecule. Will repeat the Optimization loop until the Molecule has fully converged."""
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

    def isOptimized(self, frequencies):
        """Returns a Boolean Flag as to if the Molecule is fully Optimized (Fully Optimized = All Frequencies > 0)"""
        for freq in frequencies:
            if freq < 0:
                return False

        return True
