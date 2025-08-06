import pandas as pd
import numpy as np
from openbabel import pybel

class Molecule:
    """Standard molecule class for handling coordinates and conversions between data formats."""

    #
    # Need to be Set
    #
    input: str
    """The input string representing the molecule. Can be any format supported 
    by Open Babel, but formats other than SMILES should be specified with in_format.
    """

    in_format: str
    """The format of the input string (e.g., "smi", "xyz", "pdb"). Used to determine 
    how to read the input.
    """

    def __init__(self, input, in_format="smi"):
        self.input = input
        self.in_format = in_format
        self.coords = self._get_clean_xyz()

    def _get_clean_xyz(self):
        xyz = self.convert(self.input, self.in_format, "xyz", localopt=True)
        lines = xyz.splitlines()
        return "\n".join(lines[2:])
    
    @staticmethod
    def convert(input: str, in_format: str, out_format: str, localopt=False) -> str:
        """Converts a molecule from one format to another using Open Babel. 
        Used to set initial coordinates, and can be called outside of the class as a static method.
        
        ## Parameters : \n
            input - The input string representing the molecule. Can be any format supported by Open Babel.
            in_format - The format of the input string (e.g., "smi", "xyz", "pdb").
            out_format - The desired output format (e.g., "smi", "xyz", "pdb").
            localopt - If True, performs a quick geometry optimization on the molecule

        ## Returns : \n
            None - No Return Value
        """
        mol = pybel.readstring(in_format, input)
        mol.addh()
        mol.make3D()
        if localopt or in_format == "smi":
            mol.localopt()

        return mol.write(out_format)

    def save(self, filename):
        pass


if __name__ == "__main__":
    tester = Molecule("CCO")
    print(tester.coords)
    print("==============")
    print(tester.convert("CCO", "smi", "xyz"))