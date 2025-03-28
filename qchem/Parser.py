import pandas as pd
import re
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from qchem.Molecule import Molecule
from qchem.XYZFile import XYZFile

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__))))

class OrcaOutput:
    def __init__(self, filePath: str):
        """Initialize OrcaOutput with ORCA output file path and extract all data."""
        # Store file path and read contents
        self.filePath = filePath
        with open(self.filePath, "r") as file:
            self.lines = file.readlines()

        # Extract filename without path/extension using regex
        self.name = re.search(r"[^\\]+$", self.filePath).group()[:-4]
        
        self.removalList = []

        # Extract all data from the Output
        self.extractRegexIndices()
        self.globalExtractor()

    def saveToTxt(self, outputPath: str):
        """Save all extracted data to formatted text file."""
        with open(outputPath, "w") as file:
            # Write basic information
            file.write(f"File: {self.filePath}\n")
            file.write(f"Final Single Point Energy: {self.energy}\n\n")

            # Write SCF convergence data
            file.write("SCF Energies:\n")
            for energy in self.SCFEnergies:
                file.write(f"{energy}\n\n")

            # Write timing information
            file.write("Final Timings:\n")
            if isinstance(self.finalTimings, pd.DataFrame):
                file.write(self.finalTimings.to_string(index=False))
                file.write("\n\n")

            # Write population analysis
            file.write("Mayer Population Analysis:\n")
            for i, df in enumerate(self.mayerPopulation):
                file.write(f"Step {i + 1}:\n")
                file.write(df.to_string(index=False))
                file.write("\n\n")

            # Write charge analysis
            file.write("Loewdin Charges:\n")
            for i, df in enumerate(self.loedwin):
                file.write(f"Step {i + 1}:\n")
                file.write(df.to_string(index=False))
                file.write("\n\n")

            # Write dipole information
            file.write(f"Total Dipole Moment: {self.dipole}\n")
            file.write(f"Magnitude of Dipole Moment: {self.absoluteDipole}\n\n")

            # Write Gibbs energy if available
            file.write("Gibbs Free Energy:\n")
            gibbs = self.gibbsEnergy
            if gibbs:
                file.write(f"Final Gibbs Free Energy: {gibbs[0]} {gibbs[1]}\n\n")

            # Write GOAT data if available
            if "GOAT" in self.calculationTypes:
                file.write("GOAT Calculation Results:\n\n")

                file.write("Conformer Analysis:\n")
                if not self.conformers.empty:
                    file.write(self.conformers.to_string(index=False))
                    file.write("\n\n")

                if self.GOATSummary:
                    file.write("Summary:\n")
                    file.write(f"Conformers below 3 kcal/mol: {self.GOATSummary['conformersBelow3KCal']}\n")
                    file.write(f"Lowest energy conformer: {self.GOATSummary['lowestEnergy']} Eh\n")
                    file.write(f"Sconf at 298.15 K: {self.GOATSummary['sconf']} cal/(molK)\n")
                    file.write(f"Gconf at 298.15 K: {self.GOATSummary['gconf']} kcal/mol\n\n")

    def readXYZFile(self, path: str) -> pd.DataFrame:
        """Read XYZ format atomic coordinates."""
        return pd.read_csv(path, sep=r"\s+", skiprows=1, names=["Atom", "X", "Y", "Z"], engine="python")

    def determineCalculationType(self, line):
        """Determine types of calculations in output file."""
        
        if "VIBRATIONAL FREQUENCIES" in line:
            self.calculationTypes.append("FREQ") if not "FREQ" in self.calculationTypes else None
        elif "CHEMICAL SHIELDINGS (ppm)" in line:
            self.calculationTypes.append("NMR") if not "NMR" in self.calculationTypes else None
        elif "GEOMETRY OPTIMIZATION" in line:
            self.calculationTypes.append("OPT") if not "OPT" in self.calculationTypes else None
        elif "GOAT Global Iter" in line:
            self.calculationTypes.append("GOAT") if not "GOAT" in self.calculationTypes else None
            
    def getFinalTimings(self, index) :
        """Extract computational timing information."""
        
        startIndex = index + 2
        timeLines = self.lines[startIndex:-2]

        # Process timing data into a DataFrame
        df = pd.DataFrame([line.split("...") for line in timeLines], columns=["Timing", "Time"])
        df[["Time", "B"]] = df["Time"].str.split("sec", n=1, expand=True)
        self.finalTimings = df.drop(["B"], axis=1)
        
    def getFinalEnergy(self, line):
        """Extract final single-point energy from output."""
        #if line.strip().startswith("FINAL SINGLE POINT ENERGY"):
        self.energy = float(line.split()[-1])
                    
    def getSCFEnergies(self, index):
        """Extract SCF iteration energies and convergence data."""
        
        # Skip header lines to get to actual energy value
        for j in range(index + 1, min(index + 5, len(self.lines))):
            if "Total Energy" in self.lines[j]:
                # Extract energy value
                energy = float(self.lines[j].split()[-2])
                self.SCFEnergies.append(energy)
                break
            elif any(char.isdigit() for char in self.lines[j]):
                # Direct energy value line
                energy = float(self.lines[j].split()[-2])
                self.SCFEnergies.append(energy)
                break
        

    def getMayerPopulation(self, index):
        """Extract Mayer population analysis data for atomic properties."""
        # Get total number of atoms
        
        mayerPopulations = []
        startIndex = index + 1
        endIndex = index + self.atomCount + 1
        mayerLines = self.lines[startIndex:endIndex]

        # Convert to DataFrame with atomic properties
        df = pd.DataFrame(
            [line.split()[1:8] for line in mayerLines],
            columns=["ATOM", "NA", "ZA", "QA", "VA", "BVA", "FA"],
        )
        # Convert numeric columns
        df[["NA", "ZA", "QA", "VA", "BVA", "FA"]] = df[["NA", "ZA", "QA", "VA", "BVA", "FA"]].astype(float)
        mayerPopulations.append(df)
        
        self.mayerPopulation = mayerPopulations

    def getDipoleVector(self, line):
        """Extract x, y, z components of dipole moment vector."""
        
        self.dipole = tuple(map(float, line.split()[4:]))  # Convert to floats

    def getDipoleMagnitude(self, line, index):
        """Extract magnitude of total dipole moment."""
        self.absoluteDipole =  float(line.split()[3])

    def getVibrationalFrequencies(self, index):
        """Extract vibrational frequencies from frequency calculation."""

        freqs = []
        
        startIdx = index + 5  # Skip header lines
        # Process each frequency line
        while self.lines[startIdx].strip():
            parts = self.lines[startIdx].split()
            if len(parts) >= 2:
                freqs.append(
                    {
                        "mode": int(parts[0].strip(":")),
                        "frequency": float(parts[1]),
                    }
                )
            startIdx += 1

        # Return empty DataFrame if no frequen cies found
        if len(freqs) == 0:
            self.vibrationalFrequencies = pd.DataFrame(columns=["mode", "frequency"])
            return
        
        self.vibrationalFrequencies = pd.DataFrame(freqs)

    def getGibbsEnergy(self, line):
        """Extract Gibbs free energy and units."""
        
        #if line.strip()[0:23] == "Final Gibbs free energy":
        gibbs = float(re.search(r"-?\d+\.\d+", line).group())
        unit = re.search(r"\b\w+\b$", line).group()
        self.gibbsEnergy = gibbs, unit
            
    def getSolvationEnergy(self, line, index):
        """Extract solvation energy (Eh) from output

        ## Parameters : \n
            self : OrcaOutput - Default Parameter for the Class Instance

        ## Returns : \n
            float - Solvation energy of the solute in the specified solvent (Eh)
        """
        
        #if "Gsolv" in line: 
        self.solvationEnergy = float(line.strip()[29:-9])

    def getConformerInfo(self, index):
        """Extract conformer energies and populations."""
        
        conformers = []
        startIdx = index + 2  # Skip header
        while "------" not in self.lines[startIdx]:
            startIdx += 1
        startIdx += 1  # Skip separator line

        while self.lines[startIdx].strip() and not "Conformers below" in self.lines[startIdx]:
            parts = self.lines[startIdx].split()
            if len(parts) >= 5:
                conformers.append(
                    {
                        "conformer": int(parts[0]),
                        "energy": float(parts[1]),
                        "degeneracy": int(parts[2]),
                        "totalPercent": float(parts[3]),
                        "cumulativePercent": float(parts[4]),
                    }
                )
            startIdx += 1

        self.conformerInfo = pd.DataFrame(conformers)

    def getGoatSummary(self, line):
        """Extract GOAT calculation summary."""
        
        if "Conformers below" in line:
            self.GOATSummary["conformersBelow3KCal"] = int(line.split(":")[1])
        elif "Lowest energy conformer" in line:
            self.GOATSummary["lowestEnergy"] = float(line.split(":")[1].split()[0])
        elif "Sconf at" in line:
            self.GOATSummary["sconf"] = float(line.split(":")[1].split()[0])
        elif "Gconf at" in line:
            self.GOATSummary["gconf"] = float(line.split(":")[1].split()[0])
            
    def getIRFrequencies(self, index):
        """Extract IR frequencies and intensities."""
        
        freqs = []
       
        startIdx = index + 4  # Skip header lines
        while self.lines[startIdx].strip():
            parts = self.lines[startIdx].split()
            if len(parts) >= 7:  # Mode, freq, eps, Int, T**2, TX, TY, TZ
                freqs.append(
                    {
                        "mode": int(parts[0].strip(":")),
                        "frequency": float(parts[1]),
                        "IRIntensity": float(parts[3]),  # km/mol
                    }
                )
            startIdx += 1
            
        self.IRFrequencies = pd.DataFrame(freqs)

    def getChemicalShifts(self, index):
        """Extract NMR chemical shifts."""
        
        shifts = []
        startIdx = index + 6
        while self.lines[startIdx].strip():
            parts = self.lines[startIdx].split()
            if len(parts) >= 4:
                shifts.append(
                    {
                        "atom": parts[0],
                        "nucleus": parts[1],
                        "isotropic": float(parts[2]),
                        "anisotropic": float(parts[3]),
                    }
                )
            startIdx += 1
        self.chemicalShifts = pd.DataFrame(shifts)

    def getLoewdinCharges(self, index):
        """Extract Loewdin atomic charges."""

        startIdx = index + 2
        currentCharges = []
        # Read charges until blank line
        while self.lines[startIdx].strip():
            parts = self.lines[startIdx].split()
            for i, part in enumerate(parts):
                if part == ":":
                    parts = parts[:i] + parts[i + 1 :]
            currentCharges.append({"atomNum": int(parts[0]), "atom": parts[1], "charge": float(parts[2])})
            startIdx += 1
            
        self.loedwin.append(pd.DataFrame(currentCharges))

    def extractConformers(self):
        """Extract conformer molecular structures."""
        # Get the Number of Atoms we should Expect
        if self.isFileReference():
            atomNum = XYZFile(self.molecule).atomCount
        else:
            atomNum = self.molecule.atomCount

        # Open the File
        ensembleXYZFile = open(os.path.join(self.orcaCachePath, f"{self.name}.finalensemble.xyz"))

        # Get all the Lines from the File
        allLines = ensembleXYZFile.readlines()

        # Get the Expected Length of a XYZ File
        XYZLength = atomNum + 2

        # Calculate the Number of
        moleculeCount = int((len(allLines)) / (XYZLength))

        for i in range(moleculeCount):
            # Get the Lines for a Single XYZ File
            molLines = allLines[i * XYZLength : (i + 1) * XYZLength :]

            # Make the Name
            moleculeName = f"{self.name}_Conf_{i}"

            # Load as a XYZ File
            xyz = XYZFile(molecule=molLines, name=moleculeName)

            # Convert to a Molecule
            molecule = Molecule(moleculeName, xyz)

            # Add to the Conformer List
            self.conformers.append(molecule)
            
    def extractRegexIndices(self):
        patterns = {
            "FinalEnergy": re.compile(r"^FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", re.MULTILINE),
            "SCF": re.compile(r"^TOTAL SCF ENERGY", re.MULTILINE),
            "FinalTiming": re.compile(r"^Timings for individual modules:", re.MULTILINE),
            "MayerPopStart": re.compile(r"ATOM\s+NA\s+ZA\s+QA\s+VA\s+BVA\s+FA"),
            "NumberOfAtoms": re.compile(r"^Number of atoms\s+(\d+)", re.MULTILINE),
            "Dipole": re.compile(r"^Total Dipole Moment", re.MULTILINE),
            "DipoleMagnitude": re.compile(r"^Magnitude \(a\.u\.\)", re.MULTILINE),
            "GibbsEnergy": re.compile(r"^Final Gibbs free energy", re.MULTILINE),
            "SolvationEnergy": re.compile(r"\s+Gsolv\s+=\s+([-+]?\d+\.\d+)", re.MULTILINE),
            "Loedwin": re.compile(r"^LOEWDIN ATOMIC CHARGES", re.MULTILINE),
            "VibrationalFrequencies": re.compile(r"^VIBRATIONAL FREQUENCIES", re.MULTILINE),
            "IRFrequencies": re.compile(r"^IR SPECTRUM", re.MULTILINE),
            "ChemicalShifts": re.compile(r"^CHEMICAL SHIELDING SUMMARY \(ppm\)", re.MULTILINE),
            "ConformerInfo": re.compile(r"^# Final ensemble info #", re.MULTILINE),
            "GoatConformerCount": re.compile(r"^Conformers below", re.MULTILINE),
            "GoatLowestEnergy": re.compile(r"^Lowest energy conformer", re.MULTILINE),
            "GoatSconf": re.compile(r"^Sconf at", re.MULTILINE),
            "GoatGconf": re.compile(r"^Gconf at", re.MULTILINE),
        }
        
        self.matchIndices = {key: [] for key in patterns.keys()}
        
        for key, pattern in patterns.items():
            for lineNum, line in enumerate(self.lines):
                if pattern.match(line):
                    self.matchIndices[key].append(lineNum)
        
            
    def globalExtractor(self):
        
        self.calculationTypes = []
        self.SCFEnergies = []
        self.solvationEnergy = None
        self.GOATSummary = {}
        self.loedwin = []
        self.vibrationalFrequencies = pd.DataFrame(columns=["mode", "frequency"])
        self.IRFrequencies = pd.DataFrame(columns=["mode", "frequency", "IRIntensity"])
        self.conformerInfo = pd.DataFrame(columns=["conformer", "energy", "degeneracy", "totalPercent", "cumulativePercent"])
        self.GOATSummary = None
        self.chemicalShifts = None
     
        for i, line in enumerate(self.lines):
            self.determineCalculationType(line)
            if line.strip().startswith("Number of atoms"):
                self.atomCount = int(line.split()[-1])
      
        for i in self.matchIndices["FinalEnergy"]:
            self.getFinalEnergy(self.lines[i], i)
            
        for i in self.matchIndices["SCF"]:
            self.getSCFEnergies(self.lines[i], i)
        
        for i in self.matchIndices["FinalTiming"]:
            self.getFinalTimings(self.lines[i], i)
        
        for i in self.matchIndices["MayerPopStart"]:
            self.getMayerPopulation(self.lines[i], i)
            
        for i in self.matchIndices["Dipole"]:
            self.getDipoleVector(self.lines[i], i)
        
        for i in self.matchIndices["DipoleMagnitude"]:
            self.getDipoleMagnitude(self.lines[i], i)
        
        for i in self.matchIndices["GibbsEnergy"]:
            self.getGibbsEnergy(self.lines[i], i)
        
        for i in self.matchIndices["SolvationEnergy"]:
            self.getSolvationEnergy(self.lines[i], i)
        
        for i in self.matchIndices["Loedwin"]:
            self.getLoewdinCharges(self.lines[i], i)
        
        if "FREQ" in self.calculationTypes:
            for i in self.matchIndices["VibrationalFrequencies"]:
                self.getVibrationalFrequencies(self.lines[i], i)
                
            for i in self.matchIndices["IRFrequencies"]:
                self.getIRFrequencies(self.lines[i], i)
                
        if "NMR" in self.calculationTypes:
            for i in self.matchIndices["ChemicalShifts"]:
                self.getChemicalShifts(self.lines[i], i)
        
        if "GOAT" in self.calculationTypes:
            for i in self.matchIndices["ConformerInfo"]:
                self.getConformerInfo(self.lines[i], i)
            
            for i in self.matchIndices["GoatConformerCount"]:
                self.getGoatSummary(self.lines[i], i)
        
            for i in self.matchIndices["GoatLowestEnergy"]:
                self.getGoatSummary(self.lines[i], i)
            
            for i in self.matchIndices["GoatSconf"]:
                self.getGoatSummary(self.lines[i], i)
            
            for i in self.matchIndices["GoatGconf"]:
                self.getGoatSummary(self.lines[i], i)


if __name__ == "__main__":
    for file in os.listdir(os.path.join(os.getcwd(), "tests", "test_files", "output_files")):
        if file.endswith(".out"):
            loadPath = os.path.join(os.getcwd(), "tests", "test_files", "output_files", file)
            savePath = os.path.join(
                os.getcwd(), "tests", "test_files", "output_files", "extracted", f"{file[:-4]}_parsed.txt"
            )
            output = OrcaOutput(loadPath)
            output.saveToTxt(savePath)
            
            
