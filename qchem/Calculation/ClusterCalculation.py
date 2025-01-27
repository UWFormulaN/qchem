import os
from qchem.Calculation.OrcaCalculation import OrcaCalcResult
from qchem.Calculation.OrcaInputFile import OrcaInputFile
from qchem.Calculation.OrcaCalculation import runOrcaCalculation
import multiprocessing
import time

class ClusterCalculation:
    """Class that manages multiple Orca Calculations in parallel. Will spawn / fill clusters to maximize core usage. Manages each calculations lifetimes and will fill in available calculation slots."""

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
        """Starts, Runs and Manages all Calculations assigned to the Cluster
        
        ## Parameters: \n
            self - ClusterCalculation: Default Parameter for the Class Instance
            
        ## Returns: \n
            None - No Return Value
        """
        processes:list[multiprocessing.Process] = []
        message_queue = multiprocessing.Queue()  # Create a message queue

        while self.calculations or any(p.is_alive() for p in processes):
            # Clean up finished processes
            for p in processes[:]:
                if not p.is_alive():
                    p.join()
                    processes.remove(p)
                    self.usedCores -= p.calculation.variables["cores"]

            # Start new calculations if there are available cores
            while self.calculations and self.usedCores < self.maxCores:
                calculation = self.calculations[0]
                # Check if we have Enough Cores to Spare for the Next Calculation
                if self.usedCores + calculation.variables["cores"] <= self.maxCores:
                    # Prepare and Start the Calculation
                    calculation.index = self.index
                    self.index += 1
                    p = multiprocessing.Process(target=self.runIndividualCalculation, args=(calculation,message_queue))
                    p.calculation = calculation  # Store calculation in the process
                    p.start()
                    p.is_alive()
                    processes.append(p)
                    self.usedCores += calculation.variables["cores"]
                    self.calculations.pop(0) # Maybe move this back to the top
            
            # Check for messages from the processes
            self.postMessages(message_queue)

            # Wait Half a Second before Checking Again
            time.sleep(0.5)
        
        self.postMessages(message_queue)
                
    def runIndividualCalculation(self, calculation: OrcaInputFile, messageQueue: multiprocessing.Queue):
        """Runs an Individual Calculation assigned to the Cluster. Spawns the Orca instance and waits until completion. Adds the results to the Message Queue to be released.
        
        ## Parameters: \n
            self - ClusterCalculation: Default Parameter for the Class Instance
            calculation - OrcaInputFile: The Calculations Input file
            messageQueue - multiprocessing.Queue: The Message Queue where results and the completion message will be added
            
        ## Returns: \n
            None - No Return Value
        """
        messageQueue.put(f"Starting Calculation #{self.index}")
        calcResults = runOrcaCalculation(self.name + f"_{self.index}", calculation, self.index, self.isLocal, self.STDOut, self.orcaCachePath)
        messageQueue.put(calcResults) # Store the Results in the Message Queue
        messageQueue.put(f"Completed Calculation {self.index}")
        
    def postMessages (self, messageQueue: multiprocessing.Queue):
        """Releases the Content from the Message Queue, adds Completed calculations to the appropriate property and prints the completion messages to the Terminal
        
         ## Parameters: \n
            self - ClusterCalculation: Default Parameter for the Class Instance
            messageQueue - multiprocessing.Queue: The Message Queue where results and the completion message will be added
            
        ## Returns: \n
            None - No Return Value
        """
        # Check for messages from the processes
        while not messageQueue.empty():
            message = messageQueue.get()
            if isinstance(message, OrcaCalcResult):
                self.completedCalculations.append(message)
            else:
                print(message)