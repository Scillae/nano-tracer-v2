import pickle
import os.path


def assignment_parser(s):
    assert s is not None
    s = [i.strip() for i in s.split() if i.strip()]
    eq_sign = s.index('=')
    left = s[:eq_sign]
    right = s[eq_sign + 1:]
    return left[0] if left else None, right


def parse_single_nucleotide(s):
    return [float(i.strip()) for i in s.split() if i.strip()]


def nextline(cursor):
    line = cursor.readline()
    while line == '\n':
        line = cursor.readline()
    if line == '':
        return None
    return line


def formatter(type_list, s):
    """ elements in s must be splited by space

    :param type_list: list of type convert function/ctor
    :param s: source string
    :return: list of converted elements
    """
    if not (type(s) == list or type(s) == tuple):
        s = [i.strip() for i in s.split() if i.strip()]
    total = len(s)
    idx = 0

    def recurrent(tl):
        nonlocal idx
        temp = []
        for element in tl:
            if type(element) == tuple:
                temp.append(recurrent(element))
            else:
                if idx >= total:
                    raise Exception('unmatch parameter count')
                temp.append(element(s[idx]))
                idx += 1
        return tuple(temp)

    res = recurrent(type_list)
    assert idx == len(s)
    return res


def save_load(p, obj):
    """
    Save an object if not saved, or load an object if not loaded.
    If both saved and loaded, update the saved one
    :param p: path of the saved
    :param obj: object, the loaded one
    """
    print(f'TM path: {p}')
    # Noella_v2: may use `is not` instead?
    if p != None:
        chkdir(os.path.dirname(p))  # os.path.split(p)[0]
    if p is None:
        print('TM skipping!')
        return obj
    if (obj is not None) and (not os.path.isfile(p)):
        print('TM saving!')
        pickle.dump(obj, open(p, "wb"))
        r_obj = obj
    elif (obj is None) and os.path.isfile(p):
        print('TM loading!')
        r_obj = pickle.load(open(p, "rb"))
    elif (obj is not None) and os.path.isfile(p):
        print('TM updating savepoint!')
        pickle.dump(obj, open(p, "wb"))
        r_obj = obj
    else:
        print('TM save_load both empty')
        r_obj = False
    return r_obj


def chkdir(d):
    """
    Create a directory if not existing
    :param d: path
    """
    if not os.path.isdir(d):
        os.makedirs(d)
        print(f'Created Directory: {d}')
    return


def dims_adjust(dims_ls, conf_suffix, single=True, sp_suffix=''):
    """
    Adjust the dims_ls based on conf_suffix
    :param dims_ls: the dimensions of the nanostar
    :param conf_suffix: the suffix describing nanostar's configuration
    :param single: if it's for a single nanostar
    :param sp_suffix: special suffix
    """
    dims_ls = [20, 2, 7]
    if conf_suffix.split('_')[0] == '':
        return
    if conf_suffix.split('_')[0] == '-jun':
        dims_ls[1] = int(conf_suffix.split('_')[1])
    else:
        assert 0 == 1
    return
