import pandas as pd
import pytest
from tools import utils


@pytest.fixture
def sample_xyz_path():
    return "./test_files/nano.xyz"


def test_read_xyz(sample_xyz_path):
    xyz_df = utils.read_xyz(sample_xyz_path)
    expected_df = pd.DataFrame(
        {
            "Atom": ["C", "C", "C", "C", "N"],
            "X": [-24.745720, -24.572582, -25.889140, -26.172895, -26.846360],
            "Y": [-11.253937, -12.601896, -13.190549, -11.013601, -12.212589],
            "Z": [8.115931, 8.024661, 7.996101, 8.139223, 8.080417],
        }
    )
    pd.testing.assert_frame_equal(xyz_df, expected_df)
