import gzip


def openfile(file_name, mode='rt', **kwargs):
    """open gzip or plain file"""
    if file_name.endswith('.gz'):
        file_obj = gzip.open(file_name, mode=mode, **kwargs)
    else:
        file_obj = open(file_name, mode=mode, **kwargs)
    return file_obj


def read_one_col(fn):
    """read one column file into list"""
    with openfile(fn) as f:
        return [x.strip() for x in f]