from qchem.Calculation.OrcaCalculation import OrcaCalcResult
import multiprocessing
import time

class ClusterCalculation:
    """A Class that allows you to queue multiple Orca Calculations at the same time and maximize computer usage to finish multiple calculations in Parallel"""

    maxCores: int
    """Maximum Number of Cores the Cluster can use at the same Time"""

    usedCores: int
    """The Current number of Cores being used"""

    calculations: list[OrcaCalcResult]
    """List of Calculations that are to be Started"""

    completedCalculations : list[OrcaCalcResult]
    """List of Calculations that have been Completed"""

    index: int
    """Counter Identifying a Docker Container from another. Associated with the Calculations in Queue"""

    def __init__(
        self,
        calculations: list[OrcaCalcResult],
        maxCores: int = 1,
    ):
        self.calculations = calculations
        self.maxCores = maxCores
        self.usedCores = 0
        self.index = 1
        self.completedCalculations = []

    def runCalculations(self):
        """Runs all Calculations that have been assigned to the Cluster"""

        processes:list[multiprocessing.Process] = []
        message_queue = multiprocessing.Queue()  # Create a message queue

        while self.Calculations or any(p.is_alive() for p in processes):
            # Clean up finished processes
            for p in processes[:]:
                if not p.is_alive():
                    p.join()
                    self.completedCalculations.append(p.calculation)
                    processes.remove(p)
                    self.usedCores -= p.calculation.cores

            # Start new calculations if there are available cores
            while self.calculations and self.usedCores < self.maxCores:
                calculation = self.calculations[0]
                # Check if we have Enough Cores to Spare for the Next Calculation
                if self.usedCores + calculation.cores <= self.maxCores:
                    # Prepare and Start the Calculation
                    self.calculations.pop(0)
                    calculation.index = self.index
                    self.index += 1
                    p = multiprocessing.Process(target=self.runIndividualCalculation, args=(calculation,message_queue))
                    p.calculation = calculation  # Store calculation in the process
                    p.start()
                    p.is_alive()
                    processes.append(p)
                    self.usedCores += calculation.cores
            
            # Check for messages from the processes
            while not message_queue.empty():
                message = message_queue.get()
                print(message)

            # Wait Half a Second before Checking Again
            time.sleep(0.5)
                
    def runIndividualCalculation(self, calculation: OrcaCalcResult, message_queue: multiprocessing.Queue):
        """Runs an Individual Calculation and Provides Messages Saying it's Started and Finished"""
        message_queue.put(f"Starting Calculation {calculation.getInputFileName()}")
        calculation.runCalculation()
        self.completedCalculations.append(calculation)
        message_queue.put(f"Completed Calculation {calculation.getInputFileName()}")