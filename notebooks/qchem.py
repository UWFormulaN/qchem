import os
import sys
from dotenv import load_dotenv


def setup_project_environment():
    # Load environment variables from .env file
    load_dotenv()

    # Retrieve the value of the PROJECT_ROOT environment variable
    project_root = os.getenv("PROJECT_ROOT")

    # Check if PROJECT_ROOT is None
    if project_root is None:
        print("Error: PROJECT_ROOT environment variable is not set.")
        sys.exit(1)

    # Add the data_tools directory to Python path
    sys.path.insert(0, project_root)


# Call the setup_project_environment function when qchem.py is imported
setup_project_environment()
