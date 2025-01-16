import os
import time
import pandas as pd
from qchem.Molecule import Molecule
from qchem.Parser import OrcaOutput
from .BaseOrcaCalculation import BaseOrcaCalculation
from qchem.Calculation.OrcaInputFile import OrcaInputFile
from qchem.Data.Enums import OrcaInputTemplate, OrcaCalculationType
from qchem.Calculation.OrcaCalculation import RunOrcaCalculation, OrcaCalcResult


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
        self.BasisSetFunctionalCompliant()

        self.fullOptimization = fullOptimization

    def RunCalculation(self):
        """Runs through a Geometry Optimization on the Molecule and repeats until properly converged"""

        # Single Optimization
        if not self.fullOptimization:
            self.SingleOptimization()
        else:
            self.FullOptimization()

    def SingleOptimization(self):

        print(f"Running OPT on {self.name}...")

        # Start the Timer
        startTime = time.time()

        # Run the Calculation
        calculation = RunOrcaCalculation(
            self.name, self.inputFile, isLocal=self.isLocal, STDOut=False
        )

        # Get the Output File
        outputFile = OrcaOutput(calculation.outputFilePath)

        if outputFile.vibrational_frequencies is None or (
            isinstance(outputFile.vibrational_frequencies, pd.DataFrame)
            and outputFile.vibrational_frequencies.empty
        ):
            print("No Frequencies Found! Optimization Failed!")

        else:
            # Check if the Molecule is Fully Optimized
            if self.IsOptimized(outputFile.vibrational_frequencies["frequency"]):
                self.calculationTime = time.time() - startTime
                print(
                    f"Molecule {self.name} is Optimized! ({self.ClockTime(self.calculationTime)})"
                )
                self.calculation = calculation
                self.optimizedMoleculePath = os.path.join(
                    calculation.orcaCachePath, calculation.name + ".xyz"
                )
                self.optMolecule = Molecule(self.name, self.optimizedMoleculePath)
            else:
                print(f"Molecule {self.name} is not Optimized!")

        calcTime = time.time() - startTime
        print(f"Finished OPT on {self.name} ({self.ClockTime(calcTime)})")

    def FullOptimization(self):
        # Start the Timer
        startTime = time.time()

        # Set the Optimization Index, and Initialize the Optimization Flag
        optIndex = 1
        isOptimized = False
        freqFailCount = 0

        # Full Optimization Loop
        while not isOptimized:

            iterStartTime = time.time()

            # Generate Print Statement for User on the Optimization Attempt
            print(f"Running OPT {optIndex} on {self.name}...")

            # Generate an Indexed Name for each Calculation
            calcName = self.name if (optIndex == 1) else self.name + f"_{optIndex}"

            # Run the Calculation
            calculation = RunOrcaCalculation(
                calcName, self.inputFile, isLocal=self.isLocal, STDOut=False
            )

            # Get the Output File
            outputFile = OrcaOutput(calculation.outputFilePath)

            # Check if Frequencies are Empty and cannot be found
            if outputFile.vibrational_frequencies is None or (
                isinstance(outputFile.vibrational_frequencies, pd.DataFrame)
                and outputFile.vibrational_frequencies.empty
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
                if self.IsOptimized(outputFile.vibrational_frequencies["frequency"]):
                    self.calculationTime = time.time() - startTime
                    print(
                        f"Molecule {self.name} is Optimized! ({self.ClockTime(self.calculationTime)})"
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
                f"Finished OPT {optIndex} on {self.name} ({self.ClockTime(calcTime)})"
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

    def IsOptimized(self, frequencies):
        """Determines if the Molecule is Optimized to the Valid Minimum"""
        for freq in frequencies:
            if freq < 0:
                return False

        return True
