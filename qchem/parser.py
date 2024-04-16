# This is top priority for now, get dataframes out of Orca output files for info.

import pandas as pd
import seaborn as sns

# This section is for parsing geometry files
# It should include xyz, smiles, mol, etc.
 
def read_xyz(path):
    return pd.read_csv(
        path, sep=r"\s+", skiprows=1, names=["Atom", "X", "Y", "Z"], engine="python"
    )

# This section is for parsing specific results from output files
# So far, it contains SCF Energy and Mayer Analysis

def read_SCF_Energy(output_path):
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

def read_mayer(output_path):
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

