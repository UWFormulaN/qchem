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
            
    def getFinalTimings(self, line, index) :
        """Extract computational timing information."""
        
        if index < len(self.lines) - 20 or not line.strip() == "Timings for individual modules:":
            return None
            
        startIndex = index + 2
        timeLines = self.lines[startIndex:-2]

        # Process timing data into a DataFrame
        df = pd.DataFrame([line.split("...") for line in timeLines], columns=["Timing", "Time"])
        df[["Time", "B"]] = df["Time"].str.split("sec", n=1, expand=True)
        self.finalTimings = df.drop(["B"], axis=1)
        self.removalList.append("FinalTiming")
        #self.extractors.pop("FinalTiming")
        
    def getFinalEnergy(self, line, index):
        """Extract final single-point energy from output."""
        if line.strip().startswith("FINAL SINGLE POINT ENERGY"):
           self.energy = float(line.split()[-1])
           self.removalList.append("FinalEnergy")
           #self.extractors.pop("FinalEnergy")
                    
    def getSCFEnergies(self, line, index):
        """Extract SCF iteration energies and convergence data."""
        
        if not ("TOTAL SCF ENERGY" in line or "FINAL SINGLE POINT ENERGY" in line):
            return
            
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
        self.removalList.append("SCF")
        #self.extractors.pop("SCF")
        

    def getMayerPopulation(self, line, index):
        """Extract Mayer population analysis data for atomic properties."""
        # Get total number of atoms
        if line.strip().startswith("Number of atoms"):
            self.atomCount = int(line.split()[-1])

        if not line.strip() == "ATOM       NA         ZA         QA         VA         BVA        FA":
            return
        
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
        self.removalList.append("MayerPop")
        #self.extractors.pop("MayerPop")

    def getDipoleVector(self, line, index):
        """Extract x, y, z components of dipole moment vector."""
        
        if line.strip().startswith("Total Dipole Moment"):
            self.dipole = tuple(map(float, line.split()[4:]))  # Convert to floats
            #self.extractors.pop("Dipole")
            self.removalList.append("Dipole")

    def getDipoleMagnitude(self, line, index):
        """Extract magnitude of total dipole moment."""
        if line.strip().startswith("Magnitude (a.u.)"):
            self.absoluteDipole =  float(line.split()[3])
            #self.extractors.pop("DipoleMagnitude")
            self.removalList.append("DipoleMagnitude")

    def getVibrationalFrequencies(self, line, index):
        """Extract vibrational frequencies from frequency calculation."""

        if not "FREQ" in self.calculationTypes:
            self.vibrationalFrequencies = pd.DataFrame(columns=["mode", "frequency"])
            return
        
        if not "VIBRATIONAL FREQUENCIES" in line:
            return
        
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
        self.removalList.append("VibrationalFrequencies")
        #self.extractors.pop("VibrationalFrequencies")

    def getGibbsEnergy(self, line, index):
        """Extract Gibbs free energy and units."""
        
        if line.strip()[0:23] == "Final Gibbs free energy":
            gibbs = float(re.search(r"-?\d+\.\d+", line).group())
            unit = re.search(r"\b\w+\b$", line).group()
            self.gibbsEnergy = gibbs, unit
            #self.extractors.pop("GibbsEnergy")
            self.removalList.append("GibbsEnergy")
            
    def getSolvationEnergy(self, line, index):
        """Extract solvation energy (Eh) from output

        ## Parameters : \n
            self : OrcaOutput - Default Parameter for the Class Instance

        ## Returns : \n
            float - Solvation energy of the solute in the specified solvent (Eh)
        """
        
        if "Gsolv" in line: 
            self.solvationEnergy = float(line.strip()[29:-9])
            #self.extractors.pop("SolvationEnergy")
            self.removalList.append("SolvationEnergy")

    def getConformerInfo(self, line, index):
        """Extract conformer energies and populations."""
        
        if not "GOAT" in self.calculationTypes:
            return
        
        if not "# Final ensemble info #" in line:
            return
        
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
        self.removalList.append("ConformerInfo")
        #self.extractors.pop("ConformerInfo")

    def getGoatSummary(self, line, index):
        """Extract GOAT calculation summary."""
        
        if not "GOAT" in self.calculationTypes:
            return
        
        if "Conformers below" in line:
            self.GOATSummary["conformersBelow3KCal"] = int(line.split(":")[1])
        elif "Lowest energy conformer" in line:
            self.GOATSummary["lowestEnergy"] = float(line.split(":")[1].split()[0])
        elif "Sconf at" in line:
            self.GOATSummary["sconf"] = float(line.split(":")[1].split()[0])
        elif "Gconf at" in line:
            self.GOATSummary["gconf"] = float(line.split(":")[1].split()[0])
            
    def getIRFrequencies(self, line, index):
        """Extract IR frequencies and intensities."""
        
        if not "FREQ" in self.calculationTypes:
            return
        
        if not "IR SPECTRUM" in line:
            return
        
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
        #self.extractors.pop("IRFrequencies")
        self.removalList.append("IRFrequencies")

    def getChemicalShifts(self, line, index):
        """Extract NMR chemical shifts."""
        
        if not "NMR" in self.calculationTypes:
            return
        
        if not "CHEMICAL SHIELDING SUMMARY (ppm)" in line:
            return
        
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
        #self.extractors.pop("ChemicalShifts")
        self.removalList.append("ChemicalShifts")

    def getLoewdinCharges(self, line, index):
        """Extract Loewdin atomic charges."""

        if not "LOEWDIN ATOMIC CHARGES" in line:
            return 
        
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
        #self.extractors.pop("Loedwin")
        self.removalList.append("Loedwin")

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
        
        self.extractors = {
            "SCF" : self.getSCFEnergies,
            "FinalTiming" : self.getFinalTimings,
            "FinalEnergy" : self.getFinalEnergy,
            "MayerPop" : self.getMayerPopulation,
            "Dipole" : self.getDipoleVector,
            "DipoleMagnitude" : self.getDipoleMagnitude,
            "GibbsEnergy" : self.getGibbsEnergy,
            "SolvationEnergy" : self.getSolvationEnergy,
            "Loedwin" : self.getLoewdinCharges,
        }
        
        for i, line in enumerate(self.lines):
            self.determineCalculationType(line)
            
        if "FREQ" in self.calculationTypes:
            
            self.extractors["VibrationalFrequencies"] = self.getVibrationalFrequencies
            self.extractors["IRFrequencies"] = self.getIRFrequencies
            
        if "GOAT" in self.calculationTypes:
            self.extractors["ConformerInfo"] = self.getConformerInfo
            self.extractors["GOATSummary"] = self.getGoatSummary
            
        if "NMR" in self.calculationTypes:
            self.extractors["ChemicalShifts"] = self.getChemicalShifts
            
        # Loop through all the lines in the file and grab it's line index alongside it
        for i, line in enumerate(self.lines):
            for extractor in self.extractors.keys():
                self.extractors[extractor](line, i)
            
            for removal in self.removalList:
                self.extractors.pop(removal)
                #extractor(line, i)
            self.removalList = []
            
            #self.getSCFEnergies(line, i)
            #self.getFinalTimings(line, i)
            #self.getFinalEnergy(line, i)
            #self.getMayerPopulation(line, i)
            #self.getDipoleVector(line, i)
            #self.getDipoleMagnitude(line, i)
            #self.getVibrationalFrequencies(line, i)
            #self.getGibbsEnergy(line, i)
            #self.getSolvationEnergy(line, i)
            #self.getConformerInfo(line, i)
            #self.getGoatSummary(line, i)
            #self.getIRFrequencies(line, i)
            #self.getChemicalShifts(line, i)
            #self.getLoewdinCharges(line, i)

if __name__ == "__main__":
    for file in os.listdir(os.path.join(os.getcwd(), "tests", "test_files", "output_files")):
        if file.endswith(".out"):
            loadPath = os.path.join(os.getcwd(), "tests", "test_files", "output_files", file)
            savePath = os.path.join(
                os.getcwd(), "tests", "test_files", "output_files", "extracted", f"{file[:-4]}_parsed.txt"
            )
            output = OrcaOutput(loadPath)
            output.saveToTxt(savePath)
            
            
