# Adding the project root to sys.path
import shutil
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Importing the required modules
from qchem.Molecule import Molecule
from qchem.Pipelines.Spectra import Spectra
from qchem.Calculation.GOAT import GOAT
from qchem.Calculation.Frequency import Frequency
from qchem.Data import OrcaBasisSet, OrcaDensityFunctional
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
    
    # Load the Ethane Molecule
    mol = Molecule("Ethane", os.path.join("tests", "test_files", "ethane.xyz"))
    
    # Create a Spectra Object
    spectra = Spectra(mol, OrcaBasisSet.DEF2_SVP.value, OrcaDensityFunctional.B3LYP.value, Cores, LocalTest)
    
    # Run the Spectra Calculation
    spectra.RunCalculation()
    
    # Display the Plot in a Tab
    Spectra.PlotSpectra(spectra.IRSpectra, "Ethane Test", showPlot=False)

#
# Test 2 : Use File Reference as Input for Molecular Geometry
#
def Test2 ():
    # Create the Orca Cache Directory
    if os.path.exists(os.path.join("OrcaCache", "Propane_Spectra_Ref_GEOOPT")):
        shutil.rmtree(os.path.join("OrcaCache", "Propane_Spectra_Ref_GEOOPT"))
    os.mkdir(os.path.join("OrcaCache", "Propane_Spectra_Ref_GEOOPT"))
    
    # Copy the Proper XYZ File to be used
    shutil.copy(os.path.join("tests", "test_files", "propane.xyz"), os.path.join("OrcaCache", "Propane_Spectra_Ref_GEOOPT", "propane.xyz"))
    
     # Create a Spectra Object
    spectra = Spectra("propane.xyz", OrcaBasisSet.DEF2_SVP.value, OrcaDensityFunctional.B3LYP.value, Cores, LocalTest, name="Propane_Spectra_Ref")
    
    # Run the Spectra Calculation
    spectra.RunCalculation()
    
    # Display the Plot in a Tab
    Spectra.PlotSpectra(spectra.IRSpectra, "Propane", showPlot=False)
    
    
#
# Test 1 : Load Molecule and Use Molecule as Input
#
def Test3 ():
    
    # Load the Aspirin Molecule
    mol = Molecule("Aspirin", os.path.join("tests", "test_files", "aspirin_raw.xyz"))
    
    # Create a Spectra Object
    spectra = Spectra(mol, OrcaBasisSet.DEF2_SVP.value, OrcaDensityFunctional.B3LYP.value, Cores, LocalTest)
    
    # Run the Spectra Calculation
    spectra.RunCalculation()
    
    # Display the Plot in a Tab
    spectra.PlotSpectra(spectra.IRSpectra, "Aspirin", showPlot=False)

#
# Test 1 : Load Molecule and Use Molecule as Input
#
def Test4 ():
    
    # Load the Aspirin Molecule
    mol = Molecule("Dypyrro", os.path.join("tests", "test_files", "dypyrroMethane.xyz"))
    
    # Create a Spectra Object
    spectra = Spectra(mol, OrcaBasisSet.DEF2_SVP.value, OrcaDensityFunctional.B3LYP.value, Cores, LocalTest)
    
    # Run the Spectra Calculation
    spectra.RunCalculation()
    
    # Display the Plot in a Tab
    spectra.PlotSpectra(spectra.IRSpectra, "Dypyrro", showPlot=False)

#
# Running Tests
#
Test1()
Test2()
Test3()
Test4()