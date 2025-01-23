import os
from qchem.Calculation.OrcaCalculation import OrcaCalcResult
from qchem.Calculation.OrcaInputFile import OrcaInputFile
from qchem.Calculation.OrcaCalculation import runOrcaCalculation
import multiprocessing
import time

class ClusterCalculation:
    """A Class that allows you to queue multiple Orca Calculations at the same time and maximize computer usage to finish multiple calculations in Parallel"""

    name: str
    """Name of the Cluster Calculation"""

    isLocal: bool
    """Boolean Flag to Indicate if the Calculation is Running Locally or in Docker (True = Local, False = Docker)"""

    STDOut: bool
    """Boolean Flag to Indicate if Standard Output logs should be printed"""

    maxCores: int
    """Maximum Number of Cores the Cluster can use at the same Time"""

    usedCores: int
    """The Current number of Cores being used"""

    calculations: list[OrcaInputFile]
    """List of Calculations that are to be Started"""

    completedCalculations : list[OrcaCalcResult]
    """List of Calculations that have been Completed"""

    index: int
    """Counter Identifying a Docker Container from another. Associated with the Calculations in Queue"""

    orcaCachePath: str
    """Path to the folder that stores temporary and resulting Calculation Files"""

    def __init__(
        self,
        calculations: list[OrcaInputFile],
        maxCores: int = 1,
        name: str = "ClusterCalculation",
        isLocal: bool = False,
        STDOut: bool = True,
    ):
        # Set the Variables
        self.name = name
        self.isLocal = isLocal
        self.STDOut = STDOut
        self.maxCores = maxCores
        self.usedCores = 0
        self.index = 0
        self.calculations = calculations
        self.completedCalculations = []
        self.orcaCachePath = os.path.join(os.getcwd(), "OrcaCache", name)
        

    def runCalculations(self):
        """Runs all Calculations that have been assigned to the Cluster"""

        processes:list[multiprocessing.Process] = []
        message_queue = multiprocessing.Queue()  # Create a message queue

        while self.calculations or any(p.is_alive() for p in processes):
            # Clean up finished processes
            for p in processes[:]:
                if not p.is_alive():
                    p.join()
                    self.completedCalculations.append(p.calculation)
                    processes.remove(p)
                    self.usedCores -= p.calculation.variables["cores"]

            # Start new calculations if there are available cores
            while self.calculations and self.usedCores < self.maxCores:
                calculation = self.calculations[0]
                # Check if we have Enough Cores to Spare for the Next Calculation
                if self.usedCores + calculation.variables["cores"] <= self.maxCores:
                    # Prepare and Start the Calculation
                    self.calculations.pop(0)
                    calculation.index = self.index
                    self.index += 1
                    p = multiprocessing.Process(target=self.runIndividualCalculation, args=(calculation,message_queue))
                    p.calculation = calculation  # Store calculation in the process
                    p.start()
                    p.is_alive()
                    processes.append(p)
                    self.usedCores += calculation.variables["cores"]
            
            # Check for messages from the processes
            while not message_queue.empty():
                message = message_queue.get()
                print(message)

            # Wait Half a Second before Checking Again
            time.sleep(0.5)
                
    def runIndividualCalculation(self, calculation: OrcaInputFile, message_queue: multiprocessing.Queue):
        """Runs an Individual Calculation and Provides Messages Saying it's Started and Finished"""
        message_queue.put(f"Starting Calculation #{self.index}")
        runOrcaCalculation(self.name + f"_{self.index}", calculation, self.index, self.isLocal, self.STDOut, self.orcaCachePath)
        self.completedCalculations.append(calculation)
        message_queue.put(f"Completed Calculation {self.index}")