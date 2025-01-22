import os
import subprocess
import time
from .OrcaInputFile import OrcaInputFile

class OrcaCalcResult:
    
    name: str
    """Name of the Calculation"""
    
    orcaCachePath: str
    """Path to the Calculations Directory"""
    
    def __init__ (self, name, cachePath):
        self.name = name
        self.orcaCachePath = cachePath
        self.outputFilePath = os.path.join(self.orcaCachePath, getOutputFileName(name))

def runOrcaCalculation(name, inputFile: OrcaInputFile, index: int = 1, isLocal: bool = False, STDOut: bool = True, cachePath: str = os.path.join(os.getcwd(), "OrcaCache")):
    """Runs a Orca Calculation in a Docker Container or Locally"""

    # The Cache Path for Storage
    orcaCachePath =  os.path.join(cachePath, name)

    # Make Cache Folder if it doesn't Exist
    if not os.path.exists(orcaCachePath):
        os.makedirs(orcaCachePath)

    # Save the Input File to the folder
    inputFile.saveInputFile(os.path.join(orcaCachePath, getInputFileName(name)))

    # Get the Start Time of the Calculation
    startTimer = time.time()
    
    if STDOut:
        print(f"Running Calulation : {getInputFileName(name)}")

    # Run the Calculation Locally or through a Docker Container
    if (isLocal):
        result = runLocally(orcaCachePath, name)
    else:
        result = runDockerContainer(orcaCachePath, name, index)
    
    # Get the Total Calculation time
    calculationTime = time.time() - startTimer
    
    # Post a message that an Error may have Occured
    if (result.stderr.__len__() > 0):
        print(f"WARNING Errors Maybe Occured : \n\n{result.stderr}")
    
    # If Standard Output Allowed post the Completion Message
    if STDOut:
        print(f"Calculation Complete ({clockTime(calculationTime)}) : {getInputFileName(name)}") 
    
    return OrcaCalcResult(name, orcaCachePath)

def runLocally (cachePath: str, name: str):
    """Runs the Orca Calculation Locally on your Personal Device (Requires Orca to be Installed (With all Dependencies and Extras))"""
    
    # Create the Command String
    command = ""
    
    # Windows OS
    if os.name == 'nt':
        orcaPath = subprocess.run("where orca", shell=True, text=True, capture_output=True).stdout.strip(" \n\"").removesuffix(".exe")
        command = f"cd /d {cachePath} && \"{orcaPath}\" \"{getInputFileName(name)}\" > \"{getOutputFileName(name)}\""
    else:
        # Unix based OS (Linux, Mac)
        command = f"cd \"{cachePath}\" && /Orca/orca {self.getInputFileName()} > {self.getOutputFileName()}"
    
    # Run the Orca Calculation locally
    return subprocess.run(command, shell=True, text=True, capture_output=True)
        
def runDockerContainer ( cachePath, name, index):
    """Runs the Orca Calculation inside a Docker Container (Requires Docker to be Installed (May need to update every once in a while))"""
    # Create the Command String
    command = f'docker run --name qchemorca{index} -v "{cachePath}":/home/orca mrdnalex/orca sh -c "cd /home/orca && /Orca/orca {getInputFileName(name)} > {getOutputFileName(name)}"'

    # Kill and Remove qchemorca container if it doesn't exist yet
    subprocess.run(f"docker kill qchemorca{index}" , shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
    subprocess.run(f"docker rm qchemorca{index}" , shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
    
    # Run the Calculation in a Container and wait
    result = subprocess.run(command, shell=True, text=True, capture_output=True)
    
    # Kill and Remove the Container
    subprocess.run(f"docker kill qchemorca{index}" , shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
    subprocess.run(f"docker rm qchemorca{index}" , shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
    
    return result

def clockTime(seconds):
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
    
def getInputFileName (name):
    """Returns the Input File Name with it's extension"""
    return f"{name}.inp"

def getOutputFileName (name):
    """Returns the Output File Name with it's extension"""
    return f"{name}.out"