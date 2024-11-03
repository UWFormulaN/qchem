import pandas as pd
import re


class OrcaOutput:
    def __init__(self, file_path):
        """Initialize the OrcaOutput class with the path to the ORCA output file.
        Reads the file content and processes different properties such as final timings,
        Mayer population analysis, Loewdin charges, dipole moments, etc."""
        self.file_path = file_path
        with open(self.file_path, "r") as file:
            self.lines = file.readlines()
        # Call methods to extract relevant data from the file
        self.final_timings = self.read_final_timings()
        self.mayer_population = self.read_mayer()
        self.loedwin = self.read_Loedwin()
        self.dipole = self.read_dipole()
        self.absolutedipole = self.read_absolute_dipole()
        self.energy = self.finalenergy()
        # Extract the name of the file without the path and extension
        self.name = re.search(r"[^\\]+$", self.file_path).group()[:-4]
        self.vibrational_frequencies = self.get_vibrational_frequencies() if "FREQ" in self.calc_type() else None
        self.IR_frequencies = self.get_IR_frequencies() if "FREQ" in self.calc_type() else None
        self.chemical_shifts = self.get_chemical_shifts() if "NMR" in self.calc_type() else None

    def read_xyz(self, path: str) -> pd.DataFrame:
        """Reads XYZ format file to extract atomic coordinates and returns it as a DataFrame."""
        return pd.read_csv(path, sep=r"\s+", skiprows=1, names=["Atom", "X", "Y", "Z"], engine="python")

    def read_final_timings(self) -> pd.DataFrame:
        """Extracts the final timings from the ORCA output file.
        These timings give the time spent on various computational modules."""
        for i, line in enumerate(self.lines[-20:]):
            if line.strip() == "Timings for individual modules:":
                start_index = i + 2
                time_lines = self.lines[-20:][start_index:-2]

                # Process timing data into a DataFrame
                df = pd.DataFrame([line.split("...") for line in time_lines], columns=["Timing", "Time"])
                df[["Time", "B"]] = df["Time"].str.split("sec", n=1, expand=True)
                return df.drop(["B"], axis=1)

    def finalenergy(self) -> float:
        """Extracts the final single-point energy from the ORCA output file."""
        for line in self.lines:
            if line.strip()[0:25] == "FINAL SINGLE POINT ENERGY":
                return float(line.split()[-1])

    def read_mayer(self) -> list:
        """Extracts Mayer Population Analysis data, which provides details on atomic populations,
        charges, and valence. It returns a list of DataFrames."""
        for line in self.lines:
            if line.strip()[0:15] == "Number of atoms":
                atom_count = int(line.split()[-1])

        mayer_populations = []
        for i, line in enumerate(self.lines):
            if line.strip() == "ATOM       NA         ZA         QA         VA         BVA        FA":
                start_index = i + 1
                end_index = i + atom_count + 1

                mayer_lines = self.lines[start_index:end_index]

                # Process Mayer population data into a DataFrame
                df = pd.DataFrame(
                    [line.split()[1:8] for line in mayer_lines],
                    columns=["ATOM", "NA", "ZA", "QA", "VA", "BVA", "FA"],
                )
                df[["NA", "ZA", "QA", "VA", "BVA", "FA"]] = df[["NA", "ZA", "QA", "VA", "BVA", "FA"]].astype(float)

                mayer_populations.append(df)

        return mayer_populations

    def read_Loedwin(self) -> list:
        """Extracts Loewdin atomic charges from the ORCA output file.
        Returns a list of DataFrames containing the atomic charges."""
        for line in self.lines:
            if line.strip()[0:15] == "Number of atoms":
                atom_count = int(line.split()[-1])

        loedwin = []
        start_switch = False
        end_switch = False

        for i, line in enumerate(self.lines):
            if line.strip() == "LOEWDIN ATOMIC CHARGES":
                start_index = i + 2
                start_switch = True
            if line.strip() == "LOEWDIN REDUCED ORBITAL CHARGES":
                end_index = i - 2
                end_switch = True
            if (start_switch & end_switch) == True:
                loedwin_lines = self.lines[start_index:end_index]

                # Process Loewdin charges into a DataFrame
                df = pd.DataFrame(
                    [line.split()[0:5] for line in loedwin_lines], columns=["Count", "Atom", ":", "Atomic Charge"]
                )
                df = df.drop([":"], axis=1)
                loedwin.append(df)
                start_switch = False
                end_switch = False

        return loedwin

    def read_SCF_Energy(self) -> list:
        """Extracts the energy values from SCF iterations during the self-consistent field process.
        Returns a list of DataFrames containing iteration data such as step, energy, and energy change."""
        scf_energies = []
        start_switch = False
        end_switch = False

        for i, line in enumerate(self.lines):
            if line.strip() == "SCF ITERATIONS":
                start_index = i + 1
                start_switch = True
            if line.strip() == "***Rediagonalizing the Fockian in SOSCF/NRSCF***":
                end_index = i + 1
                end_switch = True
            if (start_switch & end_switch) == True:
                SCF_lines = self.lines[start_index:end_index]

                # Process SCF energy data into a DataFrame
                df = pd.DataFrame([line.split()[0:3] for line in SCF_lines], columns=["ITER", "Energy", "Delta-E"])
                df = df.apply(pd.to_numeric, errors="coerce").dropna()
                scf_energies.append(df)

                start_switch = False
                end_switch = False

        return scf_energies

    def read_dipole(self) -> tuple:
        """Extracts the total dipole moment (vector) from the ORCA output file."""
        for line in self.lines:
            if line.strip()[0:19] == "Total Dipole Moment":
                return tuple(map(float, line.split()[4:]))

    def read_absolute_dipole(self) -> float:
        """Extracts the magnitude of the dipole moment from the ORCA output file."""
        for line in self.lines:
            if line.strip()[0:15] == "Magnitude (a.u.)":
                return float(line.split()[2])

    def read_gibbs(self) -> tuple:
        """Extracts the Gibbs free energy from the ORCA output file, returning a tuple of the energy and its unit."""
        for i, line in enumerate(self.lines):
            if line.strip()[0:23] == "Final Gibbs free energy":
                gibbs = float(re.search(r"-?\d+\.\d+", line).group())
                unit = re.search(r"\b\w+\b$", line).group()
                return gibbs, unit

    def find_partition(self, gibbs_1: tuple, gibbs_2: tuple, temp=293.15) -> float:
        """Calculates the partition coefficient between two Gibbs free energies at a given temperature."""
        delta_gibbs = gibbs_1[0] * self.convert_energy(gibbs_1[1], "Eh") - gibbs_2[0] * self.convert_energy(gibbs_2[1], "Eh")
        delta_gibbs = delta_gibbs * self.convert_energy("Eh", "J/mol")
        return -delta_gibbs / (8.314 * temp)

    def convert_energy(self, unit1: str, unit2: str) -> float:
        """Converts energy values between different units such as Eh, eV, kJ/mol, kcal/mol, and J/mol."""
        energy = {"Eh": 1, "eV": 27.2107, "kJ/mol": 2625.5, "kcal/mol": 627.503, "J/mol": 2625500}
        return energy[unit2] / energy[unit1]

    def relative_polarity(self):
        """Placeholder for a function to calculate relative polarity."""
        return

    def save_to_txt(self, output_path: str):
        """Writes all the extracted information to a text file in a human-readable format."""
        with open(output_path, "w") as file:
            file.write(f"File: {self.file_path}\n")
            file.write(f"Final Single Point Energy: {self.energy}\n\n")

            file.write("SCF Energies:\n")
            for df in self.read_SCF_Energy():
                file.write(df.to_string(index=False))
                file.write("\n\n")

            file.write("Final Timings:\n")
            if isinstance(self.final_timings, pd.DataFrame):
                file.write(self.final_timings.to_string(index=False))
                file.write("\n\n")

            file.write("Mayer Population Analysis:\n")
            for i, df in enumerate(self.mayer_population):
                file.write(f"Step {i + 1}:\n")
                file.write(df.to_string(index=False))
                file.write("\n\n")

            file.write("Loewdin Charges:\n")
            for i, df in enumerate(self.loedwin):
                file.write(f"Step {i + 1}:\n")
                file.write(df.to_string(index=False))
                file.write("\n\n")

            file.write(f"Total Dipole Moment: {self.dipole}\n")
            file.write(f"Magnitude of Dipole Moment: {self.absolutedipole}\n\n")

            file.write("Gibbs Free Energy:\n")
            gibbs = self.read_gibbs()
            if gibbs:
                file.write(f"Final Gibbs Free Energy: {gibbs[0]} {gibbs[1]}\n")

    def calc_type(self) -> list[str]:
        """Determine the type of ORCA calculation from the output file"""
        calc_types = []
        
        for line in self.lines:
            if "VIBRATIONAL FREQUENCIES" in line:
                calc_types.append("FREQ") if not "FREQ" in calc_types else None
            elif "NMR CHEMICAL SHIFTS" in line:
                calc_types.append("NMR") if not "NMR" in calc_types else None
            elif "GEOMETRY OPTIMIZATION" in line:
                calc_types.append("OPT") if not "OPT" in calc_types else None
        return calc_types

    def get_vibrational_frequencies(self) -> pd.DataFrame:
        """Extract vibrational frequencies"""
        freqs = []
        for i, line in enumerate(self.lines):
            if "VIBRATIONAL FREQUENCIES" in line:
                start_idx = i + 5  # Skip header lines
                while self.lines[start_idx].strip():
                    parts = self.lines[start_idx].split()
                    if len(parts) >= 2:  # Mode, Frequency
                        freqs.append(
                            {
                                "mode": int(parts[0].strip(":")),
                                "frequency": float(parts[1]),
                            }
                        )
                    start_idx += 1

        if (len(freqs) == 0):
            print("Returning empty DataFrame")
            return pd.DataFrame(columns=["mode", "frequency"])
        return pd.DataFrame(freqs)

    def get_IR_frequencies(self) -> pd.DataFrame:
        """Extract IR vibrational frequencies and IR intensities"""
        freqs = []
        for i, line in enumerate(self.lines):
            if "IR SPECTRUM" in line:
                start_idx = i + 4  # Skip header lines
                while self.lines[start_idx].strip():
                    parts = self.lines[start_idx].split()
                    if len(parts) >= 7:  # Mode, freq, eps, Int, T**2, TX, TY, TZ
                        freqs.append(
                            {
                                "mode": int(parts[0].strip(":")),
                                "frequency": float(parts[1]),
                                "IR_intensity": float(parts[3]),  # km/mol
                            }
                        )
                    start_idx += 1
        return pd.DataFrame(freqs)

    def get_chemical_shifts(self) -> pd.DataFrame:
        """Extract NMR chemical shifts"""
        shifts = []
        for i, line in enumerate(self.lines):
            if "CHEMICAL SHIFTS" in line:
                start_idx = i + 5
                while self.lines[start_idx].strip():
                    parts = self.lines[start_idx].split()
                    if len(parts) >= 4:
                        shifts.append(
                            {
                                "atom": parts[0],
                                "nucleus": parts[1],
                                "isotropic": float(parts[2]),
                                "anisotropic": float(parts[3]),
                            }
                        )
                    start_idx += 1
        return pd.DataFrame(shifts)

    def save_freq_data(self, output_path: str):
        """Save frequency calculation results"""
        with open(output_path, "w") as file:
            file.write(f"Frequency Calculation Results for {self.file_path}\n\n")

            if isinstance(self.frequencies, pd.DataFrame):
                file.write("Vibrational Frequencies and IR Intensities:\n")
                file.write(self.frequencies.to_string(index=False))
                file.write("\n\n")

            # Add thermochemistry data if available
            gibbs = self.read_gibbs()
            if gibbs:
                file.write(f"Gibbs Free Energy: {gibbs[0]} {gibbs[1]}\n")

    def save_nmr_data(self, output_path: str):
        """Save NMR calculation results"""
        with open(output_path, "w") as file:
            file.write(f"NMR Calculation Results for {self.file_path}\n\n")

            if isinstance(self.chemical_shifts, pd.DataFrame):
                file.write("Chemical Shifts:\n")
                file.write(self.chemical_shifts.to_string(index=False))
                file.write("\n\n")

