from enum import Enum

# Enum for Orca Basis Sets
class OrcaBasisSet(Enum):
    """Stores Common Basis Sets for Orca Calculations"""
    DEF2_SVP = "DEF2-SVP"
    """DEF2_SVP Basis Set"""
    MINI = "MINI"
    """MINI Basis Set"""

# Enum for Orca Calculation Types
class OrcaCalculationType(Enum):
    """Stores Common and Automated Calculation Types"""
    OPTIMIZATION = "OPT"
    """Used for a GeoOpt (Geometry Optimization) Calculation"""
    FREQUENCY = "FREQ"
    """Used for a Freq (Frequency) Calculation"""
    HARTREE_FOCK = "HF"
    """Used for a HF (Single Point Energy Hartree Fock) Calculation"""
    GOAT = "GOAT"
    """Used for a GOAT XTB (Global Optimizer Algorithm with XTB) Calculation"""
    GOAT_XTB = "GOAT XTB"
    """Used for a GOAT XTB (Global Optimizer Algorithm with XTB) Calculation"""

# Enum for Orca Density Functionals
class OrcaDensityFunctional(Enum):
    """Stores Common Density Functionals used for Orca Calculations"""
    B3LYP = "B3LYP"
    """B3LYP Density Functional"""
    PBE = "PBE"
    """PBE Density Functional"""

# Enum for Orca Input File Templates
class OrcaInputTemplate(Enum):
    """Stores Common Input File Templates used for Orca Calculations"""
    BASIC = "!&{calculation} &{basis} &{functional}\n*xyzfile 0 1 &{xyzfile}\n"
    """Basic Input File using File Reference"""
    BASICXYZ = "!&{calculation} &{basis} &{functional}\n* xyz 0 1 \n&{xyz}\n*"
    """Basic Input File using pasted XYZ info"""
    BASICPARALLEL = "!&{calculation} &{basis} &{functional}\n%pal nprocs &{cores} end\n*xyzfile 0 1 &{xyzfile}\n"
    """Basic Input file for Multicore Calculations using File Reference"""
    BASICXYZPARALLEL = "!&{calculation} &{basis} &{functional}\n%pal nprocs &{cores} end\n* xyz 0 1 \n&{xyz}\n*"
    """Basic Input file for Multicore Calculations using pasted XYZ Info"""

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