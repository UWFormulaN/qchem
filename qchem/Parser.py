# This is top priority for now, get dataframes out of Orca output files for info.

import pandas as pd
import re

# This section is for parsing geometry files
# It should include xyz, smiles, mol, etc.
 
def read_xyz(path:str) -> pd.core.frame.DataFrame:
    return pd.read_csv(
        path, sep=r"\s+", skiprows=1, names=["Atom", "X", "Y", "Z"], engine="python"
    )

# This section is for parsing specific results from output files
# So far, it contains end timings, SCF Energy, loedwin and Mayer Analysis

# Add class "Project" or "OrcaProject" that takes a folder, extracts all the .out files and makes them into OrcaOutput objects
# Will also need a class variable ".name" for OrcaOutput that will be the file name

class OrcaOutput:
    def __init__(
        self,
        file_path
    ):
        self.file_path = file_path
        with open(self.file_path, 'r') as file:
            self.lines = file.readlines()
        self.final_timings = self.read_final_timings()
        self.mayer_population = self.read_mayer()
        self.loedwin = self.read_Loedwin()
        self.dipole = self.read_dipole()
        self.absolutedipole = self.read_absolute_dipole()
        #self.quadrupole = self.read_quadrupole()
        self.name = re.search(r"[^\\]+$",self.file_path).group()[:-4]
        self.energy = self.finalenergy()

    def read_final_timings(self) -> pd.core.frame.DataFrame:   
        for i, line in enumerate(self.lines[-20:]):
            if line.strip() == 'Timings for individual modules:':
                start_index = i + 2
                time_lines = self.lines[-20:][start_index:-2]

                df = pd.DataFrame([line.split('...') for line in time_lines], columns=['Timing', 'Time'])
                #df = df.apply(pd.to_numeric, errors='coerce').dropna()
                df[['Time', 'B']] = df['Time'].str.split('sec', n=1, expand=True)
                return df.drop(['B'],axis=1)
    def finalenergy(self) -> float:
        for line in self.lines:
            if line.strip()[0:25] == "FINAL SINGLE POINT ENERGY":
                return float(line.split()[-1])

    def read_mayer(self) -> list:
        # Reads out Mayer Population Analysis data from an Orca ouput file

        #  NA   - Mulliken gross atomic population
        #  ZA   - Total nuclear charge
        #  QA   - Mulliken gross atomic charge
        #  VA   - Mayer's total valence
        #  BVA  - Mayer's bonded valence
        #  FA   - Mayer's free valence

        # Returns data of the form [df1, df2, ...] where df1, df2, ... are dataframes for each Mayer Population Analysis section
        # (typically start and end of geometry optimization)

        for line in self.lines:
            if line.strip()[0:15] == "Number of atoms":
                atom_count = int(line.split()[-1])

        mayer_populations = []

        for i, line in enumerate(self.lines):
            if (
                line.strip()
                == "ATOM       NA         ZA         QA         VA         BVA        FA"
            ):
                start_index = i + 1
                end_index = i + atom_count + 1

                mayer_lines = self.lines[start_index:end_index]

                df = pd.DataFrame(
                    [line.split()[1:8] for line in mayer_lines],
                    columns=["ATOM", "NA", "ZA", "QA", "VA", "BVA", "FA"],
                )

                df[["NA", "ZA", "QA", "VA", "BVA", "FA"]] = df[
                    ["NA", "ZA", "QA", "VA", "BVA", "FA"]
                ].astype(float)

                mayer_populations.append(df)

        return mayer_populations

    def read_Loedwin(self) -> list:
    # Reads out Loedwin Atomic Charges data
        for line in self.lines:
            if line.strip()[0:15] == 'Number of atoms':
                atom_count = int(line.split()[-1])

        loedwin = []
        start_switch = False
        end_switch = False

        for i, line in enumerate(self.lines):
            if line.strip() == 'LOEWDIN ATOMIC CHARGES':
                start_index = i + 2
                start_switch = True
            if line.strip() == 'LOEWDIN REDUCED ORBITAL CHARGES':
                end_index = i-2
                end_switch = True
            if (start_switch & end_switch) == True:
                loedwin_lines = self.lines[start_index:end_index]

                df = pd.DataFrame([line.split()[0:5] for line in loedwin_lines], columns=['Count', 'Atom', ':', 'Atomic Charge'])
                df = df.drop([':'],axis=1)
                loedwin.append(df)
                start_switch = False
                end_switch = False

        return loedwin
    
    def read_SCF_Energy(self) -> list:
    # Reads out Energy steps for SCF Iterations

    # ITER    - Step number
    # Energy  - Energy
    # Delta-E - Change in energy from previous step

    # Returns data of the form [df1, df2, ...] where df1, df2, ... are dataframes for each SCF Iterations section
    # (which occur at every single point energy calculation)
        scf_energies = []
        start_switch = False
        end_switch = False
    
        for i, line in enumerate(self.lines):
            if line.strip() == 'SCF ITERATIONS':
                start_index = i + 1
                start_switch = True
            if line.strip() == '***Rediagonalizing the Fockian in SOSCF/NRSCF***':
                end_index = i + 1
                end_switch = True
            if (start_switch & end_switch) == True:
                SCF_lines = self.lines[start_index:end_index]

                df = pd.DataFrame([line.split()[0:3] for line in SCF_lines], columns=['ITER', 'Energy', 'Delta-E'])
                # There is more data in here, but I don't understand the difference of before and after "Removing any level shift"
                df = df.apply(pd.to_numeric, errors='coerce').dropna()
                scf_energies.append(df)

                # This is a little silly, but it works lol
                start_switch = False
                end_switch = False

        return(scf_energies)
    
    def read_dipole(self) -> tuple:
        for line in self.lines:
            if line.strip()[0:19] == "Total Dipole Moment":
                return tuple(map(float, line.split()[4:]))
    def read_absolute_dipole(self) -> float:
        for line in self.lines:
            if line.strip()[0:15] == "Magnitude (a.u.)":
                return line.split()[2]
            
    def read_quadrupole(self) -> tuple:
        for i, line in enumerate(self.lines):
            if line.strip()[0:24] == "QUADRUPOLE MOMENT (A.U.)":
                df = pd.read_csv(self.lines[i+3:i+7], delim_whitespace=True, names=["Type", "XX", "YY", "ZZ", "XY", "XZ", "YZ"])
                return df
            

def read_final_timings(output_path: str) -> pd.core.frame.DataFrame:
    with open(output_path, 'r') as file:
        lines = file.readlines()
    
    for i, line in enumerate(lines[-20:]):
        if line.strip() == 'Timings for individual modules:':
            start_index = i + 2
            time_lines = lines[-20:][start_index:-2]

            df = pd.DataFrame([line.split('...') for line in time_lines], columns=['Timing', 'Time'])
            #df = df.apply(pd.to_numeric, errors='coerce').dropna()
            df[['Time', 'B']] = df['Time'].str.split('sec', n=1, expand=True)
            df = df.drop(['B'],axis=1)
            return df

def read_SCF_Energy(output_path: str) -> list:
# Reads out Energy steps for SCF Iterations

# ITER    - Step number
# Energy  - Energy
# Delta-E - Change in energy from previous step

# Returns data of the form [df1, df2, ...] where df1, df2, ... are dataframes for each SCF Iterations section
# (which occur at every single point energy calculation)
    with open(output_path, 'r') as file:
        lines = file.readlines()

    scf_energies = []
    start_switch = False
    end_switch = False
    
    for i, line in enumerate(lines):
        if line.strip() == 'SCF ITERATIONS':
            start_index = i + 1
            start_switch = True
        if line.strip() == '***Rediagonalizing the Fockian in SOSCF/NRSCF***':
            end_index = i + 1
            end_switch = True
        if (start_switch & end_switch) == True:
            SCF_lines = lines[start_index:end_index]

            df = pd.DataFrame([line.split()[0:3] for line in SCF_lines], columns=['ITER', 'Energy', 'Delta-E'])
            # There is more data in here, but I don't understand the difference of before and after "Removing any level shift"
            df = df.apply(pd.to_numeric, errors='coerce').dropna()
            scf_energies.append(df)

            # This is a little silly, but it works lol
            start_switch = False
            end_switch = False

    return(scf_energies)

def read_mayer(output_path: str) -> list:
    # Reads out Mayer Population Analysis data from an Orca ouput file

    #  NA   - Mulliken gross atomic population
    #  ZA   - Total nuclear charge
    #  QA   - Mulliken gross atomic charge
    #  VA   - Mayer's total valence
    #  BVA  - Mayer's bonded valence
    #  FA   - Mayer's free valence

    # Returns data of the form [df1, df2, ...] where df1, df2, ... are dataframes for each Mayer Population Analysis section
    # (typically start and end of geometry optimization)
    with open(output_path, "r") as file:
        lines = file.readlines()

    for line in lines:
        if line.strip()[0:15] == "Number of atoms":
            atom_count = int(line.split()[-1])

    mayer_populations = []

    for i, line in enumerate(lines):
        if (
            line.strip()
            == "ATOM       NA         ZA         QA         VA         BVA        FA"
        ):
            start_index = i + 1
            end_index = i + atom_count + 1

            mayer_lines = lines[start_index:end_index]

            df = pd.DataFrame(
                [line.split()[1:8] for line in mayer_lines],
                columns=["ATOM", "NA", "ZA", "QA", "VA", "BVA", "FA"],
            )

            df[["NA", "ZA", "QA", "VA", "BVA", "FA"]] = df[
                ["NA", "ZA", "QA", "VA", "BVA", "FA"]
            ].astype(float)

            mayer_populations.append(df)

    return mayer_populations

def read_Loedwin(output_path: str) -> list:
# Reads out Energy steps for SCF Iterations

# ITER    - Step number
# Energy  - Energy
# Delta-E - Change in energy from previous step

# Returns data of the form [df1, df2, ...] where df1, df2, ... are dataframes for each SCF Iterations section
# (which occur at every single point energy calculation)
    with open(output_path, 'r') as file:
        lines = file.readlines()

    for line in lines:
        if line.strip()[0:15] == 'Number of atoms':
            atom_count = int(line.split()[-1])

    loedwin = []
    start_switch = False
    end_switch = False
    
    for i, line in enumerate(lines):
        if line.strip() == 'LOEWDIN ATOMIC CHARGES':
            start_index = i + 2
            start_switch = True
        if line.strip() == 'LOEWDIN REDUCED ORBITAL CHARGES':
            end_index = i-2
            end_switch = True
        if (start_switch & end_switch) == True:
            loedwin_lines = lines[start_index:end_index]

            df = pd.DataFrame([line.split()[0:5] for line in loedwin_lines], columns=['Count', 'Atom', ':', 'Atomic Charge'])
            df = df.drop([':'],axis=1)
            loedwin.append(df)
            start_switch = False
            end_switch = False

    return(loedwin)


# This Subsection is for examining gibbs free energy data and calculating solvent partition coefficients

def read_gibbs(output_path: str) -> tuple:
    with open(output_path, 'r') as file:
        lines = file.readlines()
    
    for i, line in enumerate(lines):
        if line.strip()[0:23] == 'Final Gibbs free energy':
            gibbs = float(re.search(r'-?\d+\.\d+', line).group())
            unit = re.search(r'\b\w+\b$',line).group()
            return gibbs, unit
        
def find_partition(gibbs_1:tuple, gibbs_2:tuple, temp=293.15) -> float:
    delta_gibbs = gibbs_1[0]*convert_energy(gibbs_1[1],'Eh') - gibbs_2[0]*convert_energy(gibbs_2[1],'Eh')
    delta_gibbs = delta_gibbs*convert_energy('Eh','J/mol')
    return -delta_gibbs/(8.314*temp)
    
def convert_energy(unit1: str, unit2: str) -> float:
    energy = {'Eh': 1, 'eV': 27.2107, 'kJ/mol': 2625.5, 'kcal/mol': 627.503, 'J/mol':2625500}
    return energy[unit2]/energy[unit1]

def relative_polarity():
    return
