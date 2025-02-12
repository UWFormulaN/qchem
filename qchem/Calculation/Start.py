# import qchem as qc
from qchem.Molecule import Molecule
from qchem.Parser import OrcaOutput
# from qchem.Calculation.BaseOrcaCalculation import BaseOrcaCalculation
from qchem.Calculation.OrcaCalculation import runOrcaCalculation
from qchem.Data.Enums import OrcaInputTemplate
from qchem.Calculation.OrcaInputFile import OrcaInputFile


import numpy as np
import pandas as pd
from pathlib import Path


#Solutes = ["Methanol", "Ethanol", "Benzene","Toluene"]
Solutes = ["Methanol"]
#Methods = ["ALPB", "CPCM", "COSMORS"]
Solvents = ["Water"]
Methods = ["ALPB"]


class SolvencyCalculation:
  """Class that performs solvency calculations. Calculates solvation energies and partition coefficients, storing them in a Dict. 
   e.g
   Calculation = SolvencyCalculation(Methods=["ALPB, CPCM], Solutes=["Methanol", "Benzene", "Toluene"], Solvents=["Water", "Woctanol"])
   NOTE: You must have xyz files of the solutes in the current dir. 
   """

  CPCM_SOLVENTS = ['chlorobenzene', 'chloroform', 'chcl3', 'a-chlorotoluene', 'o-chlorotoluene', 'conductor', 'm-cresol', 'mcresol', 'o-cresol', 
                    'cyclohexane', 'cyclohexanone', 'cyclopentane', 'cyclopentanol', 'cyclopentanone', 'decalin', 'cis-decalin', 'n-decane', 'decane', 
                    'dibromomethane', 'dibutylether', 'o-dichlorobenzene', 'odichlorobenzene', 'e-1,2-dichloroethene', 'z-1,2-dichloroethene', 'dichloromethane',
                      'ch2cl2', 'dcm', 'diethylether', 'diethylether', 'diethylsulfide', 'diethylamine', 'diiodomethane', 'diisopropylether', 'diisopropylether',
                        'cis-1,2-dimethylcyclohexane', 'dimethyldisulfide', 'n,n-dimethylacetamide', 'dimethylacetamide', 'n,n-dimethylformamide', 'dimethylformamide', 
                        'dmf', 'dimethylsulfoxide', 'dmso', 'diphenylether', 'dipropylamine', 'n-dodecane', 'dodecane', 'ethanethiol', 'ethanol', 'ethylacetate', 'ethylacetate', 
                        'ethanoate', 'ethylmethanoate', 'ethylphenylether', 'ethoxybenzene', 'ethylbenzene', 'fluorobenzene', 'formamide', 'formicacid', 'n-heptane', 
                        'heptane', 'n-hexadecane', 'hexadecane', 'n-hexane', 'hexane', 'hexanoicacid', 'iodobenzene', 'iodoethane', 'iodomethane', 'isopropylbenzene',
                          'p-isopropyltoluene', 'isopropyltoluene', 'mesitylene', 'methanol', 'methylbenzoate', 'methylbutanoate', 'methylethanoate', 'methylmethanoate', 
                          'methylpropanoate', 'n-methylaniline', 'methylcyclohexane', 'n-methylformamide', 'methylformamide', 'nitrobenzene', 'phno2', 'nitroethane', 
                          'nitromethane', 'meno2', 'o-nitrotoluene', 'onitrotoluene', 'n-nonane', 'nonane', 'n-octane', 'octane', 'n-pentadecane', 'pentadecane', 'pentanal',
                            'n-pentane','pentane', 'pentanoicacid', 'pentylethanoate', 'pentylamine', 'perfluorobenzene', 'hexafluorobenzene', 'phenol', 'propanal', 'propanoicacid', 
                      'propanonitrile', 'propylethanoate', 'propylamine', 'pyridine', 'tetrachloroethene', 'c2cl4', 'tetrahydrofuran', 'thf', 'tetrahydrothiophene-s,s-dioxide', '', '', 
                      'tetrahydrothiophenedioxide', 'sulfolane',  'tetralin', 'thiophene', 'thiophenol', 'toluene', 'trans-decalin', 'tributylphosphate', 'trichloroethene', 'triethylamine',
                        'n-undecane', 'undecane', 'water', 'h2o', 'xylene', 'm-xylene', 'o-xylene', 'p-xylene']
   
  SMD_SOLVENTS = ['chlorobenzene', 'chloroform', 'chcl3', 'a-chlorotoluene', 'o-chlorotoluene', 'm-cresol', 'mcresol', 'o-cresol', 
                    'cyclohexane', 'cyclohexanone', 'cyclopentane', 'cyclopentanol', 'cyclopentanone', 'decalin', 'cis-decalin', 'n-decane', 'decane', 
                    'dibromomethane', 'dibutylether', 'o-dichlorobenzene', 'odichlorobenzene', 'e-1,2-dichloroethene', 'z-1,2-dichloroethene', 'dichloromethane',
                      'ch2cl2', 'dcm', 'diethylether', 'diethylether', 'diethylsulfide', 'diethylamine', 'diiodomethane', 'diisopropylether', 'diisopropylether',
                        'cis-1,2-dimethylcyclohexane', 'dimethyldisulfide', 'n,n-dimethylacetamide', 'dimethylacetamide', 'n,n-dimethylformamide', 'dimethylformamide', 
                        'dmf', 'dimethylsulfoxide', 'dmso', 'diphenylether', 'dipropylamine', 'n-dodecane', 'dodecane', 'ethanethiol', 'ethanol', 'ethylacetate', 'ethylacetate', 
                        'ethanoate', 'ethylmethanoate', 'ethylphenylether', 'ethoxybenzene', 'ethylbenzene', 'fluorobenzene', 'formamide', 'formicacid', 'n-heptane', 
                        'heptane', 'n-hexadecane', 'hexadecane', 'n-hexane', 'hexane', 'hexanoicacid', 'iodobenzene', 'iodoethane', 'iodomethane', 'isopropylbenzene',
                          'p-isopropyltoluene', 'isopropyltoluene', 'mesitylene', 'methanol', 'methylbenzoate', 'methylbutanoate', 'methylethanoate', 'methylmethanoate', 
                          'methylpropanoate', 'n-methylaniline', 'methylcyclohexane', 'n-methylformamide', 'methylformamide', 'nitrobenzene', 'phno2', 'nitroethane', 
                          'nitromethane', 'meno2', 'o-nitrotoluene', 'onitrotoluene', 'n-nonane', 'nonane', 'n-octane', 'octane', 'n-pentadecane', 'pentadecane', 'pentanal',
                            'n-pentane','pentane', 'pentanoicacid', 'pentylethanoate', 'pentylamine', 'perfluorobenzene', 'hexafluorobenzene', 'propanal', 'propanoicacid', 
                      'propanonitrile', 'propylethanoate', 'propylamine', 'pyridine', 'tetrachloroethene', 'c2cl4', 'tetrahydrofuran', 'thf', 'tetrahydrothiophene-s,s-dioxide', '', '', 
                      'tetrahydrothiophenedioxide', 'sulfolane',  'tetralin', 'thiophene', 'thiophenol', 'toluene', 'trans-decalin', 'tributylphosphate', 'trichloroethene', 'triethylamine',
                        'n-undecane', 'undecane', 'water', 'h2o', 'xylene', 'm-xylene', 'o-xylene', 'p-xylene']

  COSMORS_SOLVENTS = ['chlorobenzene', 'chloroform', 'chcl3', 'm-cresol', 'mcresol', 'cyclohexane', 'cyclohexanone', 'decalin', 'n-decane', 'decane', 
                    'dibutylether', 'o-dichlorobenzene', 'odichlorobenzene', 'dichloromethane', 'ch2cl2', 'dcm', 'diethylether', 'diethylether', 'diisopropylether', 'diisopropylether',
                        'n,n-dimethylacetamide', 'dimethylacetamide', 'n,n-dimethylformamide', 'dimethylformamide', 
                        'dmf', 'dimethylsulfoxide', 'dmso', 'diphenylether', 'n-dodecane', 'dodecane', 'ethanol', 'ethylacetate', 'ethylacetate', 
                        'ethanoate', 'ethylphenylether', 'ethoxybenzene', 'ethylbenzene', 'fluorobenzene', 'n-heptane', 
                        'heptane', 'n-hexadecane', 'hexadecane', 'n-hexane', 'hexane', 'iodobenzene', 'isopropylbenzene','isopropyltoluene', 'mesitylene', 'methanol',
                          'n-methylformamide', 'methylformamide', 'nitrobenzene', 'phno2', 'nitroethane', 
                          'nitromethane', 'meno2', 'n-nonane', 'nonane', 'n-octane', 'octane', 'n-pentadecane', 'pentadecane', 'n-pentane','pentane', 'perfluorobenzene', 'hexafluorobenzene', 'phenol',
                     'pyridine', 'tetrachloroethene', 'c2cl4', 'tetrahydrofuran', 'thf', 
                      'tetralin', 'toluene', 'tributylphosphate', 'triethylamine',
                        'n-undecane', 'undecane', 'water', 'h2o']

  ALPB_SOLVENTS = ['chloroform', 'chcl3', 'dichloromethane', 'ch2cl2', 'dcm', ' diethyl ether', 'diethylether', 'n-dimethylformamide', 'dimethylformamide', 
        'dmf', 'dimethylsulfoxide', 'dmso', 'ethanol', ' ethylacetate', 'ethylacetate', 'ethanoate', 'furan', 'furane', 'n-hexadecane', 'hexadecane', 
         'n-hexane', 'hexane', 'methanol', 'nitromethane', 'meno2', 'octanol(wet)', 'wetoctanol', 'woctanol', 'phenol', 'tetrahydrofuran', 'thf', 'toluene',
          'water', 'h2o' ]
  
   
  """
   ADD ERROR CHECKS
   """
   
  def __init__(self, methods: list[str], solutes: list[str], solvents: list[str], path: str, cores='16'):
      self.Methods = methods
      self.Solutes = solutes
      self.Solvents = solvents
      self.path = path
      self.Cores = cores
      self.solvencyCalc = self.calculation()
      self.SolvationEnergies = self.getSolvationEnergies()
      # self.PartitionCoeffs = self.PartitionCoeffs

  def convert_energy(self, energy: float, initial_unit: str, final_unit: str) -> float:
    """Converts energy of unit1 to unit2. Can handle different units such as Eh, eV, kJ/mol, kcal/mol, and J/mol.
    
    Params:
        energy        (float): Initial Energy
        initial_unit  (str): Units of initial energy
        final_unit    (str): Desired units after conversion
        
    Returns:
        converted_energy (float): Initial energy expressed in units of final_unit
    """
    units = {"Eh": 1, "eV": 27.2107, "kJ/mol": 2625.5, "kcal/mol": 627.503, "J/mol": 2625500}
    factor = units[final_unit] / units[initial_unit]
    converted_energy = energy*factor
    return converted_energy


  def calculation(self):

    """
    Runs Geo Opt for all the xyz files provided
    """
    #Perform geo opt on solutes
    for solute in Solutes:
        mol  = Molecule(solute, f"{solute}.xyz")
        inp = OrcaInputFile(OrcaInputTemplate.BASICXYZPARALLEL
          , calculation='OPT'
          , basis = ''
          , functional = ''
          , xyz = mol.XYZBody()
          , cores = '16'
        )
        
        "Run Geo Opt Calculation"
        runOrcaCalculation(solute, inp, isLocal=True)
        # testcalc = OrcaCalculation(solute, inp, isLocal=True)
        # testcalc.RunCalculation()
        print(f"Completed {solute} geo opt")
  

    # print(f"Completed {solute} geo opt")

  #  def GenInps_RunCalculation(self):
    """
    Generates input files and runs Orca Calculation
    """

    for solvent in self.Solvents:
      for method in self.Methods:
          for solute in self.Solutes:
              mol  = Molecule(solute, f"{solute}.xyz") 
              if method == "ALPB":
                inp = OrcaInputFile(OrcaInputTemplate.BASICXYZPARALLEL,
                                      calculation=f"XTB ALPB({solvent})",
                                      basis='',
                                      functional='',
                                      cores=self.Cores,
                                      xyz=mol.XYZBody())
              if method == "CPCM":
                inp = OrcaInputFile(OrcaInputTemplate.BASICXYZPARALLEL,
                                      calculation=f"CPCM({solvent})",
                                      basis='',
                                      functional='',
                                      cores=self.Cores,
                                      xyz=mol.XYZBody())
              if method == "COSMORS":
                inp = OrcaInputFile(OrcaInputTemplate.BASICXYZPARALLEL,
                                      calculation=f"COSMORS({solvent})",
                                      basis='',
                                      functional='',
                                      cores=self.Cores,
                                      xyz=mol.XYZBody())
              if method == "SMD":
                  inp = OrcaInputFile(OrcaInputTemplate.BASICXYZPARALLEL,
                                      calculation=f"SMD({solvent})",
                                      basis='',
                                      functional='',
                                      cores=self.Cores,
                                      xyz=mol.XYZBody())
                
              #Run the Solvency Calculation
              runOrcaCalculation(f"{solute}_{method}", inp, isLocal=True)
              # calc = OrcaCalculation(f"{solute}_{method}", inp, isLocal=True)
              # calc.RunCalculation()
              print(f"Completed {solute} in {solvent} {method} calculation")

  def getSolvationEnergies(self) -> dict:
     
     """
     Extracts the solvation energies for all the calculations - storing them in a dictionnary
     """
     SolvationEnergies = {}

     for method in self.Methods:
       SolvationEnergies[method] = []    
     
     PATH = Path(self.path)
     for folder in PATH.iterdir():
      if folder.is_dir():
        for solutefolder in folder.iterdir():
          for file in solutefolder.iterdir():
            if file.name.endswith(".out"):
              output = OrcaOutput(file.name)

              if file.name.endswith("ALPB.out"):
                SolvationEnergies["ALPB"].append((f"{file.name[:-4]}", output.solvationEnergy))
              if file.name.endswith("CPCM.out"):
               """CPCM calculations return single point energies rather than Gsolv. 
               To calculate Gsolv, we must simply find the difference in energy when in presence of solvent vs without solvent
               i.e Gsolv = Ewith - Ewithout"""

                #extracting final energy in vacuum
               vacuum_out = OrcaOutput(solutefolder.name) 
                

      

        

     return SolvationEnergies


        
    


# #Perform geo opt on solutes
# for solute in Solutes:
#     mol  = Molecule(solute, f"{solute}.xyz")
#     inp = OrcaInputFile(OrcaInputTemplate.BASICXYZPARALLEL
#       , calculation='OPT'
#       , basis = ''
#       , functional = ''
#       , xyz = mol.XYZBody()
#       , cores = '16'
#     )
#     "Run Geo Opt Calculation"
#     runOrcaCalculation(solute, inp, isLocal=True)
#     # testcalc = OrcaCalculation(solute, inp, isLocal=True)
#     # testcalc.RunCalculation()
#     print(f"Completed {solute} geo opt")


# #Generate Input Files and Run Calculation
# for solvent in Solvents:
#   for method in Methods:
#       for solute in Solutes:
#           mol  = Molecule(solute, f"{solute}.xyz") 
#           if method == "ALPB":
#             inp = OrcaInputFile(OrcaInputTemplate.BASICXYZPARALLEL,
#                                   calculation=f"XTB ALPB({solvent})",
#                                   basis='',
#                                   functional='',
#                                   cores='16',
#                                   xyz=mol.XYZBody())
#           if method == "CPCM":
#             inp = OrcaInputFile(OrcaInputTemplate.BASICXYZPARALLEL,
#                                   calculation=f"CPCM({solvent})",
#                                   basis='',
#                                   functional='',
#                                   cores='16',
#                                   xyz=mol.XYZBody())
#           if method == "COSMORS":
#             inp = OrcaInputFile(OrcaInputTemplate.BASICXYZPARALLEL,
#                                   calculation=f"COSMORS({solvent})",
#                                   basis='',
#                                   functional='',
#                                   cores='16',
#                                   xyz=mol.XYZBody())
            
#           #Run the Calculation
#           runOrcaCalculation(f"{solute}_{method}", inp, isLocal=True)
#           # calc = OrcaCalculation(f"{solute}_{method}", inp, isLocal=True)
#           # calc.RunCalculation()
#           print(f"Completed {solute} in {solvent} {method} calculation")


#Extracting/Determining Solvation Energies







  
  # if file.name.endswith("CPCM.out"):
  

  # if file.name.endswith("COSMORS.out")
  


test = SolvencyCalculation(methods=Methods, solutes=Solutes, solvents=Solvents, path="Python_calc")

print("All Calculations Complete!")

# test = SolvencyCalculation(methods=Methods, solutes=Solutes, solvents=Solvents, path= "C:/Users/camil/uWaterloo/formula_nano/Code/Python_calc")
# print(test.SolvationEnergies)