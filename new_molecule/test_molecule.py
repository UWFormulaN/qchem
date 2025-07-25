import pytest
from molecule import Molecule

@pytest.mark.parametrize(
    "input_str, in_format, expected_atoms",
    [
        ("CCO", "smi", ["C", "C", "O", "H", "H", "H", "H", "H", "H"]),
        (
            "9\n\n"
            "C 1.05821 -0.02481 -0.07213\n"
            "C 0.56076 -0.07917 1.35742\n"
            "O 1.02443 1.05468 2.07596\n"
            "H 0.71915 0.89318 -0.56329\n"
            "H 0.69813 -0.88509 -0.64322\n"
            "H 2.15275 -0.01577 -0.10113\n"
            "H -0.53281 -0.07527 1.38263\n"
            "H 0.91438 -0.98339 1.86201\n"
            "H 1.99649 1.03585 2.05767",
            "xyz",
            ["C", "C", "O", "H", "H", "H", "H", "H", "H"],
        ),
    ]
)
def test_molecule_coords(input_str, in_format, expected_atoms):
    mol = Molecule(input_str, in_format)
    lines = mol.coords.strip().splitlines()

    # Check number of atoms
    assert len(lines) == len(expected_atoms), "Atom count mismatch"

    # Check that each line starts with the correct atom
    atoms_in_output = [line.split()[0] for line in lines]
    assert atoms_in_output == expected_atoms, "Atom types or order mismatch"

