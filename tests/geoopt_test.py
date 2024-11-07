# Adding the project root to sys.path
import shutil
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Importing the required modules
from qchem.Molecule import Molecule
from qchem.Calculation.GeoOpt import GeoOpt
from qchem.Data import OrcaBasisSet, OrcaDensityFunctional
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__))))

#
# GeoOpt Tests
#

#
# SETTINGS : Modify these as needed for Testing
#
LocalTest = False
Cores = 8

#
# Test 1 : Load Molecule and Use Molecule as Input
#
def Test1 ():
    # Load the Propane Molecule
    mol = Molecule("PropaneGEOOPT", os.path.join("tests", "test_files", "propane.xyz"))

    # Create a GeoOpt Object
    geoOpt = GeoOpt(mol, OrcaBasisSet.DEF2_SVP.value, OrcaDensityFunctional.B3LYP.value, Cores, isLocal=LocalTest)

    # Run the Optimization
    geoOpt.Optimize()

#
# Test 2 : Use File Reference as Input for Molecular Geometry
#
def Test2 ():
    # Create the Orca Cache Directory
    if os.path.exists(os.path.join("OrcaCache", "Propane_GEOOPT_Ref")):
        shutil.rmtree(os.path.join("OrcaCache", "Propane_GEOOPT_Ref"))
    os.mkdir(os.path.join("OrcaCache", "Propane_GEOOPT_Ref"))
    
    # Copy the Proper XYZ File to be used
    shutil.copy(os.path.join("tests", "test_files", "propane.xyz"), os.path.join("OrcaCache", "Propane_GEOOPT_Ref", "propane.xyz"))
    
    # Create a GeoOpt Object
    geoOpt = GeoOpt("propane.xyz", OrcaBasisSet.DEF2_SVP.value, OrcaDensityFunctional.B3LYP.value, Cores, name="Propane_GEOOPT_Ref", isLocal=LocalTest)

    # Run the Optimization
    geoOpt.Optimize()

#
# Test 3 : Use a Molecule That doesn't Converge on the first Attempt
#
def Test3 ():
    # Load the Propane Molecule
    mol = Molecule("CaffeineGEOOPT", os.path.join("tests", "test_files", "caffeine.xyz"))

    # Create a GeoOpt Object
    geoOpt = GeoOpt(mol, OrcaBasisSet.MINI.value, OrcaDensityFunctional.B3LYP.value, Cores, isLocal=LocalTest)

    # Run the Optimization
    geoOpt.Optimize()

#
# Running Tests
#
Test1()
Test2()
Test3()