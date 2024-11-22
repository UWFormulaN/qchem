import pytest
from qchem import OrcaInputFile
from qchem.Data.Enums import OrcaInputTemplate


@pytest.mark.parametrize(
    "test_input",
    [
        OrcaInputTemplate.BASIC,
        "tests/test.inp",
        "!&{calculation} &{basis} &{functional}\n*xyzfile 0 1 &{xyzfile}\n",
    ],
)
def test_input_file(test_input):
    assert (
        OrcaInputFile(
            test_input,
            calculation="OPT",
            basis="def2-SVP",
            functional="PBE",
            xyzfile="aspirin.xyz",
        ).GenerateInputFile()
        == "!OPT def2-SVP PBE\n*xyzfile 0 1 aspirin.xyz\n"
    )
