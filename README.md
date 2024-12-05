QChem is a project used for Quantum Chemistry Simulations and Data Anaylsis/Visualization, developed by the Formula Nano student design team at the University of Waterloo.

Usage:

To use the package at the moment, you need to install Poetry and then install the dependencies contained in the pushed lock file. This ensures consistency among users. Ideally work should be done in Jupyter Notebooks or Python Scripts within the notebooks folder. This will allow you to import the package, allowing access to any submodules (e.g., "from qchem import utils").

Development:

Work on seperate branches for changes and submit PR's to main; squash commits when appropriate to do so. When rebasing, squashing, or changing history in any way use --force-with-lease instead of --force as an extra caution and a stale checker.

## How to use poetry:
To start, you will need to have poetry installed. DO NOT do this through pip, it needs to be installed into it's own virtual environment. 
On Windows (Powershell), the recommended command is `(Invoke-WebRequest -Uri https://install.python-poetry.org -UseBasicParsing).Content | py -`. On Linux, macOS and Windows (WSL), the recommended command is `curl -sSL https://install.python-poetry.org | python3 -`. See more at https://python-poetry.org/docs/#installing-with-the-official-installer.

When entering the project, run `poetry shell` to enter a virtual environment. This environment is built with dependencies from the lock file, which ensures consistency across users. 
