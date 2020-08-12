def int_or_float(val):
    """Used for filtering sklearn arguments that should be treated as float if between 0 and 1 and as int otherwise
    """
    if val < 1:
        return float(val)
    else:
        return int(val)