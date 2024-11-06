# Adding the project root to sys.path
import shutil
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Importing the required modules
from qchem.Molecule import Molecule
from qchem.Calculation.GOAT import GOAT
from qchem.Data import OrcaBasisSet, OrcaDensityFunctional
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__))))

#
# File for Testing GOAT
#

# GOAT Test

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
    goat = GOAT(mol, Cores, isLocal=LocalTest)

    # Run the GOAT Calculation
    goat.RunCalculation()

#
# Test 2 : Use File Reference as Input for Molecular Geometry
#
def Test2 ():
    # Create a GeoOpt Object (NOTE : You need to Create the OrcaCache Folder and the Propane_Ref Folder and paste the propane.xyz file into the folder for this test to work)\
    # Create the Orca Cache Directory
    if os.path.exists(os.path.join("OrcaCache", "Propane_GOAT_Ref")):
        os.rmdir(os.path.join("OrcaCache", "Propane_GOAT_Ref"))
    os.mkdir(os.path.join("OrcaCache", "Propane_GOAT_Ref"))
    
    # Copy the Proper XYZ File to be used
    shutil.copy(os.path.join("tests", "test_files", "propane.xyz"), os.path.join("OrcaCache", "Propane_GOAT_Ref", "propane.xyz"))
    
    # Create a GOAT Object
    goat = GOAT("propane.xyz", Cores, name="Propane_GOAT_Ref", isLocal=LocalTest)

    # Run the GOAT Calculation
    goat.RunCalculation()

#
# Test 3 : Use a Molecule That doesn't Converge on the first Attempt
#
#def Test3 ():
#    # Load the Propane Molecule
#    mol = Molecule("Caffeine", os.path.join("tests", "test_files", "caffeine.xyz"))
#
#    # Create a GeoOpt Object
#    geoOpt = GeoOpt(mol, OrcaBasisSet.MINI.value, OrcaDensityFunctional.B3LYP.value, Cores, isLocal=LocalTest)
#
#    # Run the Optimization
#    geoOpt.Optimize()

#
# Running Tests
#
#Test1()
Test2()
#Test3()