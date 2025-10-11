"""
Flat file import.

Importation of flat data file for model and mechanism data.
Currently only supports experimental data import.

"""

import pandas as pd
import os
from os.path import join


def read_flat_data(filename_data: str, data_type: str = 'dynamic', debug: bool = False) -> pd.DataFrame:
    """
    Reads a flat data file, which uses tsv delimiters and specialized headers, into a pandas DataFrame.
    Currently only dynamic data are supported.

    Parameters
    ----------
    filename_data : str
        The filename (possibly including path info) of the flat file.
    data_type : str, optional
        The type of data in the file ('static' or 'dynamic').
        Defaults to 'dynamic'.
    debug : bool, optional
        If True, prints debugging information.

    Returns
    -------
    pd.DataFrame
        Pandas DataFrame containing the processed data from the flat file.
    None
        Returns None if there's an error.
    """

    if data_type == 'static':
        # for right now, flat files only support dynamic data
        import warnings
        warnings.warn("Warning: flat files only currently support dynamic data.")
    elif data_type == 'dynamic':
        df_dat = _tsv_to_df(filename_data, debug=debug)

        # first three lines are header and meta data, with row 1 set to header
        # Row 1: Metadata, Time, and Column identifiers, using ids from model enclosed in [] with
        #        initial conditions identified by _0 appended outside the brackets (e.g., [atp]_0)
        # Row 2: data types: t (time), c (compounds), e (enzyme), id (experiment id)
        # Row 3: problem usage info: meta, time or t, independent or i, dependent or d, ignore or g.
        #        independent, dependent, ignore are used with c (compounds) and e (enzymes).
        #        time is used with t. Only should occur once
        #        meta is for all metadata such as id and any user defined items, which processing skips over
        # Row 4+: data

        #print(df_dat)

        # find data columns that contain numeric values (even if ignored)
        col_values = [col for col in df_dat.columns if df_dat.loc[1, col] in ["time", "t", "independent", "i",
                                                                              "dependent", "d", "ignore", "g"]]
        #print (col_values)

        # make sure all items that are values are floats
        for col in col_values:
            df_dat.loc[df_dat.index > 1, col] = (
                df_dat.loc[df_dat.index > 1, col].astype("float64"))
        #print(df_dat)

    else:
        import warnings
        warnings.warn("Warning: Unsupported data type.")

    return df_dat


def _tsv_to_df(filename_data: str, debug: bool = False) -> pd.DataFrame or None:
    """
    Internal utility function that reads a tsv file into a Pandas DataFrame.
    Performs basic processing to remove comments and blank lines.

    Parameters
    ----------
    filename_data : str
        The filename (possibly including path info) of the flat file.
    debug : bool, optional
        If True, prints debugging information.

    Returns
    -------
    pd.DataFrame
        Pandas DataFrame containing the minimally processed data from the flat file.
    None
        Returns None if there's an error.
    """
    try:
        if debug:
            print(f"Beginning to read file {filename_data}")
        with open(filename_data, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"Error: File not found {filename_data}")
        return None
    except Exception as e:
        print(f"Error {e}")
        return None

    # ignore comment lines that start with # or % and lines that only contain whitespace
    comment_signifier = ('#', '%')
    data_lines = [line for line in lines if not any(line.startswith(char) for char in comment_signifier) and
                                              line.strip()]

    if not data_lines:
        print("File {filename_data} is empty or contains only comments and/or blank lines.")
        return None

    # process the data lines into a dataframe
    try:
        import io
        df = pd.read_csv(io.StringIO("".join(data_lines)), sep='\t', header=0, engine='python')
    except Exception as e:
        print(f"Error {e} when processing {filename_data}")
        return None
    return df
