import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__))))

# Importing the required modules
from qchem.Molecule import Molecule
from qchem.Calculation.GeoOpt import GeoOpt
from qchem.Data import OrcaBasisSet, OrcaDensityFunctional
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__))))
from qchem.Parser import OrcaOutput
from qchem.Pipelines.Spectra import Spectra


if __name__ == "__main__":
    # We're looking to define our molecule from the xyz file we generated as a qchem object
    # The Molecule class takes in the name and the path to the xyz file


    mol = Molecule("exampleMolecule", "example.xyz") # Replace this with your actual molecule file and name

    """
    Now that we have the molecule defined, we can create a GeoOpt object to then run our calculation.

    GeoOpting the molecule will find the geometry with the lowest energy by using gradient descent. To avoid 
    getting stuck at a saddlepoint (gradient is zero but not a minimum), we will also run a 'frequency' calculation afterwards.
    """


    # We'll input our molecule as well as the theory level we want to optimize at.

    # We can also parallelize the calculation by setting cores. Orca uses openmpi for parallel CPU processing.

    # isLocal=True means that the calculation will use orca that is installed locally on your machine.
    # False means that it will call a docker container to run the calculation.
    
    # If we send the file through Discorca, Discorca will run it locally, so we will still put 'isLocal=True'

    # We're going to set fullOptimization to False, meaning that if the calculation gets stuck on a saddlepoint ## Ryan will write more here

    geoOpt = GeoOpt(mol, basis="def2-svp", functional="B3LYP", cores=1, isLocal=True, fullOptimization=False)
    
    # Now we can run the calculation by calling the runCalculation method from the geoOpt object.
    # This will create an input file, run the calculation on the machine, and create an output file with the results.
    geoOpt.runCalculation()

    #We can then go and find the output file to do more work with it.
    output = OrcaOutput(geoOpt.outputFilePath)

    # We'll first look at vibrational frequencies. Let's output that here

    print(output.vibrationalFrequencies)

    #Spectra.plotSpectra(output.IRFrequencies)


# ALSO Install.txt is needed, also xyz file




