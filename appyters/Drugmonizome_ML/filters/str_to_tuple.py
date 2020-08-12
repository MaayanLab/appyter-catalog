def str_to_tuple(s):
    return tuple(int(val) for val in s.strip('()').split(',') if len(val.strip()) > 0)