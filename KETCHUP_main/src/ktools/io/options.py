"""
options function

collects basic functions that work with reading options

"""

def read_options_file(fn: str) -> dict:
    """
    Reads option which are in stored in YAML file. Only performs basic checking, and each module needs to
    validate the items therein.

    Parameters
    ----------
    fn : str
        The filename of the YAML file to be read and processed.

    Returns
    -------
    dict
        A dictionary containing the options read from the file.
        Returns an empty dictionary if the file is empty or an error occurs.
    """

    import yaml, os

    with open(fn, 'r') as file:
        try:
            data = yaml.safe_load(file)
            # print(f"{data = }")
        except yaml.YAMLError as error_str:
            print(f"Syntax error in YAML file: {error_str}")
            data = {}

   # note: each specific module needs to validate the items in the options dictionary.

    return data
