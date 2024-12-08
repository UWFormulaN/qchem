import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import pytest
from qchem.Parser import OrcaOutput

# Test file paths
TEST_OUTPUT_DIR = os.path.join("tests", "test_files", "output_files")
ASPIRIN_FTIR = os.path.join(TEST_OUTPUT_DIR, "aspirin_ftir.out")


def test_file_loading():
    """Test basic file loading and initialization"""
    output = OrcaOutput(ASPIRIN_FTIR)
    assert output.file_path == ASPIRIN_FTIR
    assert output.name == "aspirin_ftir"
    assert len(output.lines) > 0


def test_calculation_type_detection():
    """Test detection of calculation types"""
    output = OrcaOutput(ASPIRIN_FTIR)
    assert "FREQ" in output.calc_types
    assert len(output.calc_types) > 0


def test_energy_extraction():
    """Test extraction of energy values"""
    output = OrcaOutput(ASPIRIN_FTIR)
    assert isinstance(output.energy, float)
    assert isinstance(output.scf_energies, list)
    assert len(output.scf_energies) > 0


def test_timing_extraction():
    """Test extraction of timing information"""
    output = OrcaOutput(ASPIRIN_FTIR)
    assert output.final_timings is not None
    assert "Time" in output.final_timings.columns
    assert "Timing" in output.final_timings.columns


@pytest.mark.parametrize("attribute", ["mayer_population", "loedwin", "dipole", "absolutedipole", "vibrational_frequencies", "IR_frequencies"])
def test_property_extraction(attribute):
    """Test extraction of various molecular properties"""
    output = OrcaOutput(ASPIRIN_FTIR)
    assert hasattr(output, attribute)
    assert getattr(output, attribute) is not None


def test_gibbs_energy():
    """Test extraction of Gibbs free energy"""
    output = OrcaOutput(ASPIRIN_FTIR)
    gibbs = output.get_gibbs_energy()
    assert isinstance(gibbs, tuple)
    assert len(gibbs) == 2
    assert isinstance(gibbs[0], float)
    assert isinstance(gibbs[1], str)


def test_save_to_txt(tmp_path):
    """Test saving parsed data to text file"""
    output = OrcaOutput(ASPIRIN_FTIR)
    save_path = os.path.join(tmp_path, "test_output.txt")
    output.save_to_txt(save_path)
    assert os.path.exists(save_path)
    with open(save_path, "r") as f:
        content = f.read()
        assert "Final Single Point Energy" in content
        assert "SCF Energies" in content
        assert "Final Timings" in content


def test_vibrational_frequencies():
    """Test extraction of vibrational frequencies"""
    output = OrcaOutput(ASPIRIN_FTIR)
    freqs = output.get_vibrational_frequencies()
    assert not freqs.empty
    assert "mode" in freqs.columns
    assert "frequency" in freqs.columns
    assert len(freqs) > 0


def test_dipole_moments():
    """Test extraction of dipole moments"""
    output = OrcaOutput(ASPIRIN_FTIR)
    vector = output.get_dipole_vector()
    magnitude = output.get_dipole_magnitude()
    assert isinstance(vector, tuple)
    assert len(vector) == 3
    assert all(isinstance(x, float) for x in vector)
    assert isinstance(magnitude, float)
