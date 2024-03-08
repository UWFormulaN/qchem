import pandas as pd
import seaborn as sns


def read_xyz(path):
    return pd.read_csv(
        path, sep=r"\s+", skiprows=1, names=["Atom", "X", "Y", "Z"], engine="python"
    )

def hello():
    return
