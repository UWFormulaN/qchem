from qchem.OrcaCalculation import OrcaCalculation
import multiprocessing
import time

class ClusterCalculation:
    """A Class that allows you to queue multiple Orca Calculations at the same time and maximize computer usage to finish multiple calculations in Parallel"""

    MaxCores: int
    """Maximum Number of Cores the Cluster can use at the same Time"""

    UsedCores: int
    """The Current number of Cores being used"""

    Calculations: list[OrcaCalculation]
    """List of Calculations that are to be Started"""

    CompletedCalculations : list[OrcaCalculation]
    """List of Calculations that have been Completed"""

    Index: int
    """Counter Identifying a Docker Container from another. Associated with the Calculations in Queue"""

    def __init__(
        self,
        calculations: list[OrcaCalculation],
        maxCores: int = 1,
    ):
        self.Calculations = calculations
        self.MaxCores = maxCores
        self.UsedCores = 0
        self.Index = 1
        self.CompletedCalculations = []

    def RunCalculations(self):
        """Runs all Calculations that have been assigned to the Cluster"""

        processes:list[multiprocessing.Process] = []
        message_queue = multiprocessing.Queue()  # Create a message queue

        while self.Calculations or any(p.is_alive() for p in processes):
            # Clean up finished processes
            for p in processes[:]:
                if not p.is_alive():
                    p.join()
                    self.CompletedCalculations.append(p.calculation)
                    processes.remove(p)
                    self.UsedCores -= p.calculation.Cores

            # Start new calculations if there are available cores
            while self.Calculations and self.UsedCores < self.MaxCores:
                calculation = self.Calculations[0]
                # Check if we have Enough Cores to Spare for the Next Calculation
                if self.UsedCores + calculation.Cores <= self.MaxCores:
                    # Prepare and Start the Calculation
                    self.Calculations.pop(0)
                    calculation.Index = self.Index
                    self.Index += 1
                    p = multiprocessing.Process(target=self.RunIndividualCalculation, args=(calculation,message_queue))
                    p.calculation = calculation  # Store calculation in the process
                    p.start()
                    p.is_alive()
                    processes.append(p)
                    self.UsedCores += calculation.Cores
            
            # Check for messages from the processes
            while not message_queue.empty():
                message = message_queue.get()
                print(message)

            # Wait Half a Second before Checking Again
            time.sleep(0.5)
                
    def RunIndividualCalculation(self, calculation: OrcaCalculation, message_queue: multiprocessing.Queue):
        """Runs an Individual Calculation and Provides Messages Saying it's Started and Finished"""
        message_queue.put(f"Starting Calculation {calculation.GetInputFileName()}")
        calculation.RunCalculation()
        self.CompletedCalculations.append(calculation)
        message_queue.put(f"Completed Calculation {calculation.GetInputFileName()}")