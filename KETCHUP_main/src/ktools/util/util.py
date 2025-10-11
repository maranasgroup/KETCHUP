""" 
util functions

collects basic functions that provide general utility 

""" 

def extract_from_brackets(text: str) -> str or None:
    """
    Extracts the substring between the first "[" and the last "]" in a 
    string. Used for example to get id strings.

    Parameters
    ----------
    text : str
        Input string.

    Returns
    -------
    str
        The substring between the first "[" and the last "]" (if found).
    None
        Otherwise returns None.
    """
    start_bracket = text.find("[")
    end_bracket = text.rfind("]")
    
    if (start_bracket == -1) or (end_bracket == -1) or \
       (start_bracket >= end_bracket):
        return None
    else:
        return text[start_bracket + 1:end_bracket]


