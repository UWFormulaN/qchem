import pytest
import numpy as np
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

@pytest.mark.parametrize(
    "input_str, in_format",
    [
        ("CCO", "smi"),
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
        ),
    ]
)
def test_no_duplicate_coordinates(input_str, in_format):
    mol = Molecule(input_str, in_format)
    coords = []

    for line in mol.coords.strip().splitlines():
        parts = line.split()
        if len(parts) != 4:
            continue  # skip malformed lines
        xyz = tuple(map(float, parts[1:]))  # ignore atom symbol
        coords.append(xyz)

    # Convert to numpy array for vectorized distance checks
    coords_np = np.array(coords)

    # Check for exact duplicates
    unique_coords = np.unique(coords_np, axis=0)
    assert len(unique_coords) == len(coords), "Duplicate coordinates found!"

    # Optional: Check that no pair of atoms are closer than 0.1 Å (excluding bonded cases)
    for i in range(len(coords_np)):
        for j in range(i + 1, len(coords_np)):
            dist = np.linalg.norm(coords_np[i] - coords_np[j])
            assert dist > 0.1, f"Atoms {i} and {j} are too close: {dist:.3f} Å"
