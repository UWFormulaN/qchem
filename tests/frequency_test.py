# Adding the project root to sys.path
import shutil
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Importing the required modules
from qchem.Molecule import Molecule
from qchem.Calculation.GOAT import GOAT
from qchem.Calculation.Frequency import Frequency
from qchem.Data import OrcaBasisSet, OrcaDensityFunctional
from qchem.Parser import OrcaOutput
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__))))

#
# Frequency Tests
#

#
# SETTINGS : Modify these as needed for Testing
#
LocalTest = False # I Suggest using False here as XTB may be required and that is for sure installed on Docker Container
Cores = 1

#
# Test 1 : Load Molecule and Use Molecule as Input
#
def Test1 ():
    # Load the Propane Molecule
    mol = Molecule("PropaneFREQ", os.path.join("tests", "test_files", "propane.xyz"))

    # Define the Frequency Object
    frequency = Frequency(mol, cores=Cores, isLocal=LocalTest, basis=OrcaBasisSet.DEF2_SVP.value, functional=OrcaDensityFunctional.B3LYP.value)
    
    # Run the Frequency Calculation
    frequency.RunCalculation()

#
# Test 2 : Use File Reference as Input for Molecular Geometry
#
def Test2 ():
    # Create the Orca Cache Directory
    if os.path.exists(os.path.join("OrcaCache", "Propane_FREQ_Ref")):
        shutil.rmtree(os.path.join("OrcaCache", "Propane_FREQ_Ref"))
    os.mkdir(os.path.join("OrcaCache", "Propane_FREQ_Ref"))
    
    # Copy the Proper XYZ File to be used
    shutil.copy(os.path.join("tests", "test_files", "propane.xyz"), os.path.join("OrcaCache", "Propane_FREQ_Ref", "propane.xyz"))
    
    # Define the Frequency Object
    frequency = Frequency("propane.xyz", cores=Cores, isLocal=LocalTest, basis=OrcaBasisSet.DEF2_SVP.value, functional=OrcaDensityFunctional.B3LYP.value, name="Propane_FREQ_Ref")
    
    # Run the Frequency Calculation
    frequency.RunCalculation()
    
#
# Running Tests
#
Test1()
Test2()