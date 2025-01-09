import os
import subprocess
import time
from .OrcaInputFile import OrcaInputFile
from types import SimpleNamespace

class OrcaCalcResult:
    
    name: str
    
    orcaCachePath: str
    
    def __init__ (self, name, cachePath):
        self.name = name
        self.orcaCachePath = cachePath
        self.outputFilePath = os.path.join(self.orcaCachePath, GetOutputFileName(name))

def RunOrcaCalculation(name, inputFile: OrcaInputFile, index: int = 1, isLocal: bool = False, STDOut: bool = True):
    """Runs a Orca Calculation in a Docker Container"""

    orcaCache = "OrcaCache"
    orcaCachePath =  os.path.join(os.getcwd(), orcaCache, name)

    # Make Cache Folder if it doesn't Exist
    if not os.path.exists(orcaCachePath):
        os.makedirs(orcaCachePath)

    # Save the Input File to the folder
    inputFile.SaveInputFile(os.path.join(orcaCachePath, GetInputFileName(name)))

    # Get the Start Time of the Calculation
    startTimer = time.time()
    
    if STDOut:
        print(f"Running Calulation : {GetInputFileName(name)}")

    # Run the Calculation Locally or through a Docker Container
    if (isLocal):
        result = RunLocally(orcaCachePath, name)
    else:
        result = RunDockerContainer(orcaCachePath, name, index)
    
    # Get the Total Calculation time
    calculationTime = time.time() - startTimer
    
    # Post a message that an Error may have Occured
    if (result.stderr.__len__() > 0):
        print(f"WARNING Errors Maybe Occured : \n\n{result.stderr}")
    
    # If Standard Output Allowed post the Completion Message
    if STDOut:
        print(f"Calculation Complete ({ClockTime(calculationTime)}) : {GetInputFileName(name)}") 
    
    return OrcaCalcResult(name, orcaCachePath)

def RunLocally (cachePath: str, name: str):
    
    # Create the Command String
    command = ""
    
    # Windows OS
    if os.name == 'nt':
        orcaPath = subprocess.run("where orca", shell=True, text=True, capture_output=True).stdout.strip(" \n\"").removesuffix(".exe")
        command = f"cd /d {cachePath} && \"{orcaPath}\" \"{GetInputFileName(name)}\" > \"{GetOutputFileName(name)}\""
    else:
        # Unix based OS (Linux, Mac)
        command = f"cd \"{cachePath}\" && /Orca/orca {self.GetInputFileName()} > {self.GetOutputFileName()}"
    
    # Run the Orca Calculation locally
    return subprocess.run(command, shell=True, text=True, capture_output=True)
        
def RunDockerContainer ( cachePath, name, index):
    # Create the Command String
    command = f'docker run --name qchemorca{index} -v "{cachePath}":/home/orca mrdnalex/orca sh -c "cd /home/orca && /Orca/orca {GetInputFileName(name)} > {GetOutputFileName(name)}"'

    # Kill and Remove qchemorca container if it doesn't exist yet
    subprocess.run(f"docker kill qchemorca{index}" , shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
    subprocess.run(f"docker rm qchemorca{index}" , shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
    
    # Run the Calculation in a Container and wait
    result = subprocess.run(command, shell=True, text=True, capture_output=True)
    
    # Kill and Remove the Container
    subprocess.run(f"docker kill qchemorca{index}" , shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
    subprocess.run(f"docker rm qchemorca{index}" , shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
    
    return result

def ClockTime(seconds):
        """Converts Seconds to a Human Readable Time String"""
        # Convert Seconds to Hours, Minutes, and Seconds
        days = seconds // 86400
        hours = (seconds % 86400) // 3600
        minutes = (seconds % 3600) // 60
        remainingSeconds = seconds % 60

        # Generate the Time String
        parts = []
        if days > 0:
            parts.append(f"{int(hours)} day{'s' if days > 1 else ''}")
        if hours > 0:
            parts.append(f"{int(hours)} hour{'s' if hours > 1 else ''}")
        if minutes > 0:
            parts.append(f"{int(minutes)} minute{'s' if minutes > 1 else ''}")
        if remainingSeconds > 0:
            parts.append(
                f"{int(remainingSeconds)} second{'s' if remainingSeconds > 1 else ''}"
            )

        # Return the Time String
        return ", ".join(parts) if parts else "0 seconds"
    
def GetInputFileName (name):
    """Returns the Input File Name with it's extension"""
    return f"{name}.inp"

def GetOutputFileName (name):
    """Returns the Output File Name with it's extension"""
    return f"{name}.out"