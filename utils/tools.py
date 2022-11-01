import numpy as np
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
    """ 
    elements in s must be splited by space
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


def dims_adjust(ns_dims, conf_suffix, single=True, sp_suffix=''):
    """
    Adjust the ns_dims based on conf_suffix
    :param ns_dims: the dimensions of the nanostar
    :param conf_suffix: the suffix describing nanostar's configuration
    :param single: if it's for a single nanostar
    :param sp_suffix: special suffix
    """
    ns_dims = [20, 2, 7]
    if conf_suffix.split('_')[0] == '':
        return
    if conf_suffix.split('_')[0] == '-jun':
        ns_dims[1] = int(conf_suffix.split('_')[1])
    else:
        raise Exception("Unknown Conformation Suffix!")

def data_query(data, query_vals, vals_types=None):
    default_vals = {'int': 0, 'str': '', 'dict': None, 'list': None}
    return_vals = []
    for i, val_name in enumerate(query_vals):
        if val_name in data:
            val = data[val_name] 
        else:
            val = default_vals[vals_types[i]] if vals_types is not None else None
        return_vals.append(val)
    return return_vals

def dist(t1, t2):
    return np.sqrt(np.sum(np.square(np.array(t1) - np.array(t2))))

def get_ns_params(arm_num: int):
    '''
    Parameters for generating time-averaged patch angle vs time plots.
    All are empirical values currently.
    '''
    if arm_num == 3:
        time_window_width = 11 # 15
        stacking_min_length = 10 # 
        stacking_crit_ang = 120 # 
        stacking_crit_rmsd = 12 # 
        nonstacking_min_length = 5 # 
        nonstacking_crit_ang = 105 # 
        nonstacking_crit_rmsd = 15 # 
        ns_struc = {'arm_number':3, 'linked_PA': [(0,1),(0,2),(1,2)], 'PA' : [(0,1),(0,2),(1,2)], 'pairing_linked':[((0,1),(2,3)),((0,3),(1,2))], 'pairing_unlinked':[((0,2),(1,3))]}
    elif arm_num == 4:
        time_window_width = 11 # 15
        stacking_min_length = 10 # 
        stacking_crit_ang = 120 # 
        stacking_crit_rmsd = 12 # 
        nonstacking_min_length = 5 # 
        nonstacking_crit_ang = 65 # 
        nonstacking_crit_rmsd = 15 # 
        ns_struc = {'arm_number':4, 'linked_PA': [(0,1),(0,3),(1,2),(2,3)], 'PA': [(0,1),(0,3),(1,2),(2,3),(0,2),(1,3)], 'pairing_linked':[((0,1),(2,3)),((0,3),(1,2))], 'pairing_unlinked':[((0,2),(1,3))]}
    elif arm_num == 5:
        time_window_width = 11 # 15
        stacking_min_length = 10 # useful when CoM is at the center of PJ.
        stacking_crit_ang = 120 #  # should be 130?
        stacking_crit_rmsd = 12 # 
        nonstacking_min_length = 5 # 
        nonstacking_crit_ang = 65 # 
        nonstacking_crit_rmsd = 15 # 
        ns_struc = {'arm_number':5, 'linked_PA': [(0,1),(0,4),(1,2),(2,3),(3,4)], 'PA' : [(0,1),(0,4),(1,2),(2,3),(3,4),(0,2),(1,3),(2,4),(0,3),(1,4)], 'pairing_linked':[((0,1),(2,3)),((0,3),(1,2))], 'pairing_unlinked':[((0,2),(1,3))]}
    elif arm_num == 6:
        time_window_width = 11 # 15
        stacking_min_length = 10 # useful when CoM is at the center of PJ.
        stacking_crit_ang = 120 # 
        stacking_crit_rmsd = 12 # 
        nonstacking_min_length = 5 # 
        nonstacking_crit_ang = 50 # 
        nonstacking_crit_rmsd = 15 # 
        ns_struc = {'arm_number':6, 'linked_PA': [(0,1),(0,5),(1,2),(2,3),(3,4),(4,5)], 'PA' : [(0,1),(0,5),(1,2),(2,3),(3,4),(4,5),(0,2),(1,3),(2,4),(3,5),(0,3),(1,4),(2,5),(0,4),(1,5)], 'pairing_linked':[((0,1),(2,3)),((0,3),(1,2))], 'pairing_unlinked':[((0,2),(1,3))]}
    elif arm_num == 2:
        time_window_width = 11 # 15
        stacking_min_length = 10 # 
        stacking_crit_ang = 120 # 
        stacking_crit_rmsd = 12 # 
        nonstacking_min_length = 5 # 
        nonstacking_crit_ang = 105 # 
        nonstacking_crit_rmsd = 15 # 
        ns_struc = {'arm_number':2, 'linked_PA': [(0,1)], 'PA' : [(0,1)], 'pairing_linked':[], 'pairing_unlinked':[]}

    else:
        assert 0==1    
    assert time_window_width % 2 == 1 # must be odd
    return time_window_width, stacking_min_length, stacking_crit_ang, stacking_crit_rmsd, nonstacking_min_length, nonstacking_crit_ang, nonstacking_crit_rmsd, ns_struc
