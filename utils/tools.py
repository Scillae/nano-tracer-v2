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