import os
import subprocess
import time
from .OrcaInputFile import OrcaInputFile

class OrcaCalcResult:

    name: str
    """Name of the Calculation"""

    orcaCachePath: str
    """Path to the Calculations Directory"""

    def __init__(self, name, cachePath):
        self.name = name
        self.orcaCachePath = cachePath
        self.outputFilePath = os.path.join(self.orcaCachePath, getOutputFileName(name))


def runOrcaCalculation(
    name: str,
    inputFile: OrcaInputFile,
    index: int = 1,
    isLocal: bool = False,
    STDOut: bool = True,
    cachePath: str = os.path.join(os.getcwd(), "OrcaCache"),
):
    """Default Function that is exposed and Used to Run a Calculation using Orca. Will Dispatch the Calculation Locally or through Docker based off the provided parameters

    ## Parameters : \n
        name : str - Name of the Calculation, used for the Name of the Directory and the Input and Output File \n
        inputFile : OrcaInputFile - Input file describing the information for the Orca Calculation to run \n
        index : int - Number to identify individual Docker Orca Calculations running in parallel \n
        isLocal : bool - Boolean flag to indicate if the calculation runs locally or in Docker (True = Local, False = Docker) \n
        STDOut : bool - Boolean flag to indicate if Standard Output logs should be printed \n
        cachePath : str - Path to the folder that stores temporary and resulting Calculation Files

    ## Returns : \n
        None - No Return Value
    """

    # The Cache Path for Storage
    orcaCachePath = os.path.join(cachePath, name)

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
    if isLocal:
        result = runLocally(name, orcaCachePath)
    else:
        result = runDockerContainer(name, index, orcaCachePath)

    # Get the Total Calculation time
    calculationTime = time.time() - startTimer

    # Post a message that an Error may have Occured
    if result.stderr.__len__() > 0:
        print(f"WARNING Errors Maybe Occured : \n\n{result.stderr}")

    # If Standard Output Allowed post the Completion Message
    if STDOut:
        print(
            f"Calculation Complete ({clockTime(calculationTime)}) : {getInputFileName(name)}"
        )

    return OrcaCalcResult(name, orcaCachePath)


def runLocally(name: str, cachePath: str):
    """Runs the Orca Calculation using a Local Installation of Orca Quantum Computing Software, Requires Orca to be Installed with all Additional Dependencies
    
    ## Parameters : \n
        name : str - Name of the Calculation, used for the Name of the Directory and the Input and Output File \n
        cachePath : str - Path to the folder that stores temporary and resulting Calculation Files
        
    ## Returns : \n
        subprocess.CompletedProcess - Resulting Completed Subprocess Object of the Calculation Execution
    """

    # Create the Command String
    command = ""

    # Windows OS
    if os.name == "nt":
        orcaPath = (
            subprocess.run("where orca", shell=True, text=True, capture_output=True)
            .stdout.strip(' \n"')
            .removesuffix(".exe")
        )
        command = f'cd /d {cachePath} && "{orcaPath}" "{getInputFileName(name)}" > "{getOutputFileName(name)}"'
    else:
        # Unix based OS (Linux, Mac)
        command = f'cd "{cachePath}" && /Orca/orca {getInputFileName()} > {getOutputFileName()}'

    # Run the Orca Calculation locally
    return subprocess.run(command, shell=True, text=True, capture_output=True)


def runDockerContainer(name: str, index: int, cachePath: str):
    """Runs the Orca Calculation using a Docker Container, Requires Docker to be Installed. If the Image is not Downloaded, it will automatically be Downloaded before the Calculation Starts
    
    ## Parameters : \n
        name : str - Name of the Calculation, used for the Name of the Directory and the Input and Output File \n
        index : int - Number to identify individual Docker Orca Calculations running in parallel \n
        cachePath : str - Path to the folder that stores temporary and resulting Calculation Files
    
    ## Returns : \n
        subprocess.CompletedProcess - Resulting Completed Subprocess Object of the Calculation Execution
    """
    # Create the Command String
    command = f'docker run --name qchemorca{index} -v "{cachePath}":/home/orca mrdnalex/orca sh -c "cd /home/orca && /Orca/orca {getInputFileName(name)} > {getOutputFileName(name)}"'

    # Kill and Remove qchemorca container if it doesn't exist yet
    subprocess.run(
        f"docker kill qchemorca{index}",
        shell=True,
        stderr=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
    )
    subprocess.run(
        f"docker rm qchemorca{index}",
        shell=True,
        stderr=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
    )

    # Run the Calculation in a Container and wait
    result = subprocess.run(command, shell=True, text=True, capture_output=True)

    # Kill and Remove the Container
    subprocess.run(
        f"docker kill qchemorca{index}",
        shell=True,
        stderr=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
    )
    subprocess.run(
        f"docker rm qchemorca{index}",
        shell=True,
        stderr=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
    )

    return result


def clockTime(seconds):
    """Converts Calculation Time Seconds to Human Readable Clock Format

    ## Parameters : \n
        self - Default Parameter for the Class Instance \n
        seconds : int - Number of Seconds to Convert to Clock Format

    ## Returns : \n
        str - The time in a clock format (x days : y hours : z mins : a sec)
    """
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


def getInputFileName(name):
    """Gives the Name of the Input File 
    
    ## Parameters : \n
        name : str - Name of the Input file (Without File Extension)
    
    ## Returns : \n
        str - Name of the Input File with the Extension
    """
    return f"{name}.inp"


def getOutputFileName(name):
    """Gives the Name of the Output File 
    
    ## Parameters : \n
        name : str - Name of the Output file (Without File Extension)
    
    ## Returns : \n
        str - Name of the Output File with the Extension"""
    return f"{name}.out"
