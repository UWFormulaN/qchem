from ..Data.Enums import OrcaInputTemplate


class OrcaInputFile:
    """
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

    InputFileContents: str
    """The Generated Input File as a String"""

    def __init__(self, template: str | OrcaInputTemplate, **variables):
        self.template = template
        self.variables = variables
        self.InputFileContents = self.generateInputFile()

    def generateInputFile(self) -> str:
        """Generates the input file content by replacing placeholders with actual values."""
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
        """Saves the generated input file content to a specified path."""
        with open(filePath, "w") as file:
            file.write(self.InputFileContents)
