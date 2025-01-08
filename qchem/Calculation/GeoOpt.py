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

    #
    # Need to be Set
    #
    calculationType: str = f"{OrcaCalculationType.OPTIMIZATION.value} {OrcaCalculationType.FREQUENCY.value}"

    optMolecule: Molecule
    """The Final Optimized Molecule"""

    optimizedMoleculePath: str
    """Path to the Optimized Molecule File"""

    calculation: OrcaCalculation
    """Reference to the Last Orca Calculation Object for the GeoOpt"""

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
        
        #self.variables["calculation"] = f"{OrcaCalculationType.OPTIMIZATION.value} {OrcaCalculationType.FREQUENCY.value}"
        
        #self.inputFile = OrcaInputFile(self.template, **self.variables)

    def RunCalculation(self):
        """Runs through a Geometry Optimization on the Molecule and repeats until properly converged"""
        
        startTime = time.time()
        optIndex = 1
        isOptimized = False
        freqFailCount = 0
        
        # Get Default Values
        #OPTtemplate = ""
        #xyzMol = ""
        
        #calculationType = f"{OrcaCalculationType.OPTIMIZATION.value} {OrcaCalculationType.FREQUENCY.value}"

        

        ## Determine the Proper Template to use based on the Molecule Type
        #if self.IsFileReference():
        #    OPTtemplate = OrcaInputTemplate.BASICPARALLEL
        #    xyzMol = self.molecule
#
        #    inputFile = OrcaInputFile(
        #        OPTtemplate,
        #        calculation=calculationType,
        #        basis=self.basisSet,
        #        functional=self.functional,
        #        cores=self.cores,
        #        xyzfile=xyzMol,
        #    )
        #else:
        #    OPTtemplate = OrcaInputTemplate.BASICXYZPARALLEL
        #    xyzMol = self.molecule.XYZBody()
#
        #    inputFile = OrcaInputFile(
        #        OPTtemplate,
        #        calculation=calculationType,
        #        basis=self.basisSet,
        #        functional=self.functional,
        #        cores=self.cores,
        #        xyz=xyzMol,
        #    )

        # Add a toggle for Single Or Full Optimization

        # Start the Optimization Loop
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

            # Rune the Calculation
            # calculation.RunCalculation()

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
            
            self.template = OrcaInputTemplate.BASICXYZPARALLEL
            self.variables["xyz"] = Molecule(self.name, self.optimizedMoleculePath).XYZBody()
            optIndex += 1
            
            #OPTtemplate = OrcaInputTemplate.BASICXYZPARALLEL
            #xyzMol = Molecule(self.name, self.optimizedMoleculePath).XYZBody()
            #optIndex += 1

            # Generate the Input File
            self.inputFile = OrcaInputFile(self.template, **self.variables)
            
            #inputFile = OrcaInputFile(
            #    OPTtemplate,
            #    calculation=calculationType,
            #    basis=self.basisSet,
            #    functional=self.functional,
            #    cores=self.cores,
            #    xyz=xyzMol,
            #)

            # Create the Calculation Object
            # calculation = OrcaCalculation(self.name + f"_{optIndex}", inputFile, isLocal=self.isLocal, stdout=False)

    def IsOptimized(self, frequencies):
        """Determines if the Molecule is Optimized to the Valid Minimum"""
        for freq in frequencies:
            if freq < 0:
                return False

        return True
