from enum import Enum

# Enum for Orca Basis Sets
class OrcaBasisSet(Enum):
    DEF2_SVP = "DEF2-SVP"

# Enum for Orca Calculation Types
class OrcaCalculationType(Enum):
    OPTIMIZATION = "OPT"
    INFRARED_SPECTRUM = "FREQ"
    HARTREE_FOCK = "HF"

# Enum for Orca Density Functionals
class OrcaDensityFunctional(Enum):
    B3LYP = "B3LYP",
    PBE = "PBE"

# Enum for Orca Input File Templates
class OrcaInputTemplate(Enum):
    BASIC = "!&{calculation} &{basis} &{functional}\n*xyzfile 0 1 &{xyzfile}\n"

#     Karlsruhe basis sets
# Some of the various valence adaptations of Karlsruhe basis sets[9] are briefly described below.

# def2-SV(P) – Split valence with polarization functions on heavy atoms (not hydrogen)
# def2-SVP – Split valence polarization
# def2-SVPD – Split valence polarization with diffuse functions
# def2-TZVP – Valence triple-zeta polarization
# def2-TZVPD – Valence triple-zeta polarization with diffuse functions
# def2-TZVPP – Valence triple-zeta with two sets of polarization functions
# def2-TZVPPD – Valence triple-zeta with two sets of polarization functions and a set of diffuse functions
# def2-QZVP – Valence quadruple-zeta polarization
# def2-QZVPD – Valence quadruple-zeta polarization with diffuse functions
# def2-QZVPP – Valence quadruple-zeta with two sets of polarization functions
# def2-QZVPPD – Valence quadruple-zeta with two sets of polarization functions and a set of diffuse functions