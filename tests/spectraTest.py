# Importing the required modules
from qchem.Pipelines.Spectra import Spectra
import os
import pytest

#
# Spectra Tests
#

@pytest.mark.parametrize(
    "file",
    [
        "Ethane_Spectra.csv",
        "Ethane_Spectra_NoHeader.csv"
    ],
)
def LoadTest (file):
    assert(not Spectra.LoadSpectra(os.path.join("test_files", "Spectra", file)) == None)