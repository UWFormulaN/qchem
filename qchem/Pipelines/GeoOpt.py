
from qchem.Molecule import Molecule
from qchem.OrcaInputFile import OrcaInputFile
from qchem.Data.Enums import OrcaInputTemplate
from qchem.OrcaCalculation import OrcaCalculation

class GeoOpt:
    
    molecule: Molecule | str
    """The Molecule to be Optimized"""
    
    optMolecule: Molecule
    """The Final Optimized Molecule"""
    
    basisSet: str
    """The Basis Set to be used for the Optimization"""
    
    functional: str
    """The Density Functional to be used for the Optimization"""
    
    isLocal : bool
    """Boolean Flag to determine if the Optimization is Local or not (Local = Runs on Device, Non-Local = Runs in Docker Container)"""
    
    name: str
    """Name of the Molecule being GeoOptimized"""
    
    def __init__(self, molecule, basisSet, functional, cores:int = 1, isLocal:bool = False, name:str = ""):
        
        # Check if Values are empty or of Wrong Type
        if (not (molecule and isinstance(molecule, (str, Molecule)))):
            raise ValueError("Molecule is not defined! Provide a Path to the XYZ file or a Molecule Object")
        
        if (not (basisSet and isinstance(basisSet, (str)))):
            raise ValueError("BasisSet is not defined! Provide the Name of the Basis Set as a String")
        
        if (not functional and isinstance(functional, (str))):
            raise ValueError("Functional is not defined! Provide the Name of the Functional as a String")
            
        if (not isinstance(cores, (int))):
            raise ValueError("Cores is not an Integer! Provide an Integer Value")
            
        if (not isinstance(isLocal, (bool))):
            raise ValueError("IsLocal is not a Boolean! Provide a Boolean Value")
        
        if (not isinstance(name, (str))):
            raise ValueError("Name is not a String! Provide a String Value")
        
        # Set the Name of the Molecule
        if (name == ""):
            if (isinstance(molecule, (Molecule))):
                self.name = molecule.Name
            else:
                print("No Name Provided for the Molecule, using Default Name: GeoOptMolecule")
                self.name = "GeoOptMolecule"
        else:
            self.name = name
            
        self.molecule = molecule
        self.basisSet = basisSet
        self.functional = functional
        self.cores = cores
        self.isLocal = isLocal
        
        
    def IsFileReference(self):
        """Determines if the Molecule is stored as a File Reference, or is a direct Molecule Object"""
        if (isinstance(self.molecule, (str))):
            return True
        else:
            return False
        
    def Optimize (self):
        # Create proper input file for the Molecule
        # Run the Optimization
        # Perform FREQ Analysis, then check if all the Vibrational Frequencies are positive, if so we have a truly optimized molecule
        
        # Get Default Values
        OPTtemplate = ""
        xyzMol = ""
        
        # Determine the Proper Template to use based on the Molecule Type
        if (self.IsFileReference()):
            OPTtemplate = OrcaInputTemplate.GEOOPT
            xyzMol = self.molecule
        else:
            OPTtemplate = OrcaInputTemplate.GEOOPTXYZ
            xyzMol = self.molecule.XYZBody()
            
        # Generate the Input File
        inputFile = OrcaInputFile(OPTtemplate, 
                                basis=self.basisSet,
                                 functional=self.functional,
                                 cores = self.cores,
                                 xyz=xyzMol)
        
        # Create the Calculation Object
        calculation = OrcaCalculation(self.name, inputFile, isLocal=self.isLocal)
        
        # Rune the Calculation
        calculation.RunCalculation()
        
        
        
        
        
        
        
        
        
        
        
        
    
    
    
    