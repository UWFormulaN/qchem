import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import pytest
from qchem.Parser import OrcaOutput

# Test file paths
TEST_OUTPUT_DIR = os.path.join("tests", "test_files", "output_files")
ASPIRIN_FTIR = os.path.join(TEST_OUTPUT_DIR, "aspirin_ftir.out")


def testFileLoading():
    """Test basic file loading and initialization"""
    output = OrcaOutput(ASPIRIN_FTIR)
    assert output.filePath == ASPIRIN_FTIR
    assert output.name == "aspirin_ftir"
    assert len(output.lines) > 0


def testCalculationTypeDetection():
    """Test detection of calculation types"""
    output = OrcaOutput(ASPIRIN_FTIR)
    assert "FREQ" in output.calculationTypes
    assert len(output.calculationTypes) > 0


def testEnergyExtraction():
    """Test extraction of energy values"""
    output = OrcaOutput(ASPIRIN_FTIR)
    assert isinstance(output.energy, float)
    assert isinstance(output.SCFEnergies, list)
    assert len(output.SCFEnergies) > 0


def testTimingExtraction():
    """Test extraction of timing information"""
    output = OrcaOutput(ASPIRIN_FTIR)
    assert output.finalTimings is not None
    assert "Time" in output.finalTimings.columns
    assert "Timing" in output.finalTimings.columns


@pytest.mark.parametrize("attribute", ["mayerPopulation", "loedwin", "dipole", "absolutedipole", "vibrationalFrequencies", "IRFrequencies"])
def testPropertyExtraction(attribute):
    """Test extraction of various molecular properties"""
    output = OrcaOutput(ASPIRIN_FTIR)
    assert hasattr(output, attribute)
    assert getattr(output, attribute) is not None


def testGibbsEnergy():
    """Test extraction of Gibbs free energy"""
    output = OrcaOutput(ASPIRIN_FTIR)
    gibbs = output.getGibbsEnergy()
    assert isinstance(gibbs, tuple)
    assert len(gibbs) == 2
    assert isinstance(gibbs[0], float)
    assert isinstance(gibbs[1], str)


def testSaveToTxt(tmp_path):
    """Test saving parsed data to text file"""
    output = OrcaOutput(ASPIRIN_FTIR)
    save_path = os.path.join(tmp_path, "test_output.txt")
    output.saveToTxt(save_path)
    assert os.path.exists(save_path)
    with open(save_path, "r") as f:
        content = f.read()
        assert "Final Single Point Energy" in content
        assert "SCF Energies" in content
        assert "Final Timings" in content


def testVibrationalFrequencies():
    """Test extraction of vibrational frequencies"""
    output = OrcaOutput(ASPIRIN_FTIR)
    freqs = output.getVibrationalFrequencies()
    assert not freqs.empty
    assert "mode" in freqs.columns
    assert "frequency" in freqs.columns
    assert len(freqs) > 0


def testDipoleMoments():
    """Test extraction of dipole moments"""
    output = OrcaOutput(ASPIRIN_FTIR)
    vector = output.getDipoleVector()
    magnitude = output.getDipoleMagnitude()
    assert isinstance(vector, tuple)
    assert len(vector) == 3
    assert all(isinstance(x, float) for x in vector)
    assert isinstance(magnitude, float)
