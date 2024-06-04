QChem is a project used for Quantum Chemistry Simulations and Data Anaylsis/Visualization, developed by the Formula Nano student design team at the University of Waterloo.

Usage:

To use the package at the moment, you need to install Poetry and then install the dependencies contained in the pushed lock file. This ensures consistency among users. Ideally work should be done in Jupyter Notebooks or Python Scripts within the notebooks folder. This will allow you to import qchem, allowing access to any submodules (e.g., "from qchem import utils").

Development:

Work on seperate branches for changes and submit PR's to main; squash commits when appropriate to do so. When rebasing, squashing, or changing history in any way use --force-with-lease instead of --force as an extra caution and a stale checker.
