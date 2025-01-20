from typing import Any
from ..Data.Enums import OrcaInputTemplate


class OrcaInputFile:
    """
    Represents and Creates the Input Files for Orca Calculations. Uses a Template File describing the Files form and fills out the Information
    
    Template is an input file with missing variables in the form of &{variable_name}
    kwargs is a dictionary with variable names as keys and their values as values
    ## Examples:

    !SP &{basis} PBE
    *xyzfile 0 1 aspirin.xyz

    with variables={'basis': 'def2-SVP'}
    will be converted to

    !SP def2-SVP PBE
    *xyzfile 0 1 aspirin.xyz

    ## Example with template file:

    tester = OrcaInputFile(OrcaInputTemplate.BASIC
      , calculation='OPT'
      , basis='def2-SVP'
      , functional='PBE'
      , xyzfile='aspirin.xyz'
    )
    """

    inputFileContents: str
    """Input File that’s been generated as a Single String"""
    
    template: str | OrcaInputTemplate
    """Template that descibes the Calculation to be run and the Variables to be filled in"""
    
    variables: dict[str, Any]
    """Dictionary of Variables to be filled in the Template"""

    def __init__(self, template: str | OrcaInputTemplate, **variables):
        self.template = template
        self.variables = variables
        self.inputFileContents = self.generateInputFile()

    def generateInputFile(self) -> str:
        """Generates the Input File by replacing the placeholders with the defines Values
        
        ## Parameters : \n
            self - Default Parameter for the Class Instance
            
        ## Returns : \n
            str - The Input File as a Single String
        """
        if isinstance(self.template, OrcaInputTemplate):
            inputContent = self.template.value
        elif self.template[-4:] == ".inp":
            with open(self.template, "r") as file:
                inputContent = file.read()
        else:
            inputContent = self.template

        for key, value in self.variables.items():
            placeholder = f"&{{{key}}}"
            inputContent = inputContent.replace(placeholder, str(value))
        return inputContent

    def saveInputFile(self, filePath: str):
        """Saves the Generated Input files content to the specified path
        
        ## Parameters : \n
            self - Default Parameter for the Class Instance \n
            filePath : str - Path the Input File is saved to
        
        ## Returns : \n
            None - No Return Value
        """
        with open(filePath, "w") as file:
            file.write(self.inputFileContents)
