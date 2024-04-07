import pandas as pd
import seaborn as sns


def read_xyz(path):
    return pd.read_csv(
        path, sep=r"\s+", skiprows=1, names=["Atom", "X", "Y", "Z"], engine="python"
    )


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
