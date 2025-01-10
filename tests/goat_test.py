# Adding the project root to sys.path
import shutil
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Importing the required modules
from qchem.Molecule import Molecule
from qchem.Calculation.GOAT import GOAT
from qchem.Parser import OrcaOutput
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__))))

#
# GOAT Tests
#

#
# SETTINGS : Modify these as needed for Testing
#
LocalTest = False # I Suggest using False here as XTB may be required and that is for sure installed on Docker Container
Cores = 8

#
# Test 1 : Load Molecule and Use Molecule as Input
#
def Test1 ():
    # Load the Propane Molecule
    mol = Molecule("PropaneGOAT", os.path.join("tests", "test_files", "propane.xyz"))

    # Create a GOAT Object
    goat = GOAT(mol, cores=Cores, isLocal=LocalTest, index=1)

    # Run the GOAT Calculation
    goat.RunCalculation()

#
# Test 2 : Use File Reference as Input for Molecular Geometry
#
def Test2 ():
    # Create the Orca Cache Directory
    if os.path.exists(os.path.join("OrcaCache", "Propane_GOAT_Ref")):
        shutil.rmtree(os.path.join("OrcaCache", "Propane_GOAT_Ref"))
    os.mkdir(os.path.join("OrcaCache", "Propane_GOAT_Ref"))
    
    # Copy the Proper XYZ File to be used
    shutil.copy(os.path.join("tests", "test_files", "propane.xyz"), os.path.join("OrcaCache", "Propane_GOAT_Ref", "propane.xyz"))
    
    # Create a GOAT Object
    goat = GOAT("propane.xyz", cores=Cores, name="Propane_GOAT_Ref", isLocal=LocalTest, index=2)

    # Run the GOAT Calculation
    goat.RunCalculation()
    
#
# Running Tests
#
if __name__ == "__main__":
    Test1()
    Test2()