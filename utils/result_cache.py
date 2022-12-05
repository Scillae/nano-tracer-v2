import pickle
import os


def SL_result_cache(data: dict, varname: str, SL_mode: str):
    """
    result cache interface
    :workspace data: expected to contain enough information for specifying path (and to contain the content to be saved if SL_mode is save)
    :param varname: the variable (abbr) name. For example:
            'ns' --> nanostars';
            'pa' --> 'patch angle'.
    :param SL_mode: indicating how SL works.
            'load': load cached data into data['SL_content']. False if not cached.
            'save': save data['SL_content'] into path generated by generate_path.
    """

    data['varname'] = varname

    var_path = generate_path(data, 'var_path')
    if SL_mode == 'load':
        data['SL_content'] = save_load(var_path, None)
    elif SL_mode == 'save':
        _ = save_load(var_path, data['SL_content'])
        data['SL_content'] = None
    else:
        raise Exception("Unknown SL_mode!")


def generate_path(data: dict, path_name: str):
    """
    generate path according data
    :param data: expected to contain enough information for specifying path
    :param path_name: name of the path being acquired
    :return: path specified by `path_name`
    """
    arms = data['arm_number'] if 'arm_number' in data else 0
    temp = data['temp'] if 'temp' in data else 0
    conc = data['conc'] if 'conc' in data else 0
    varname = data['varname'] if 'varname' in data else 0
    sp_suffix = data['sp_suffix'] if 'sp_suffix' in data else ''
    conf_suffix = data['conf_suffix'] if 'conf_suffix' in data else ''
    flag_suffix = data['flag_suffix'] if 'flag_suffix' in data else ''

    label = f'{arms}arms@({temp}C,{conc}M){conf_suffix}{flag_suffix}{sp_suffix}'
    loose_lbl = f'{temp}C-{conc}M-GPU{sp_suffix}'
    top_path = f'D:/Streaming-tmp/ox-sync/simul-inputs-{arms}arms{conf_suffix}{flag_suffix}/{arms}arm-rods-clustered{conf_suffix}.top'
    traj_path = f'D:/Streaming-tmp/ox-sync/simul-inputs-{arms}arms{conf_suffix}{flag_suffix}/{loose_lbl}/trajectory.dat'
    if flag_suffix in ['-cenT', '-cenToxDNA2']:
        top_path = f'D:/Streaming-tmp/ox-sync/simul-inputs-{arms}arms{conf_suffix}{flag_suffix}/{arms}arm-rods{conf_suffix}-cenT.top'
    # savepath = f'data/result_cache_files/{flag_suffix}/{arms}arms{conf_suffix}/{loose_lbl}/{label}'
    savepath = f'data/result_cache_files{flag_suffix if flag_suffix else "/"}{arms}arms{conf_suffix}/{loose_lbl}/{label}'
    if varname != 0:
        plotpath = f'results/{arms}arms{conf_suffix}{flag_suffix}/{varname}/{varname}_hist-{label}.png'
        var_path = f'{savepath}.{varname}'

    if path_name == 'var_path' and varname != 0:
        return var_path
    elif path_name == 'plot_path' and varname != 0:
        return plotpath
    elif path_name == 'top_path':
        return top_path
    elif path_name == 'traj_path':
        return traj_path
    elif path_name == 'label':
        return label
    elif path_name == 'save_path':
        return savepath
    else:
        assert 0 == 1


def save_load(p: str, obj):
    """
    Save an object if not saved, or load an object if not loaded.
    If both saved and loaded, update the saved one
    :param p: path of the saved
    :param obj: object, the loaded one
    :return: object, either passed in or loaded
    """
    print(f'save_load: path is {p}')
    if p is not None:
        chkdir(os.path.dirname(p))
    if p is None:
        raise Exception("save_load: no path indicated!")
    if (obj is not None) and (not os.path.isfile(p)):
        print('save_load: saving!')
        pickle.dump(obj, open(p, "wb"))
        r_obj = obj
    elif (obj is None) and os.path.isfile(p):
        print('save_load: loading!')
        r_obj = pickle.load(open(p, "rb"))
    elif (obj is not None) and os.path.isfile(p):
        print('save_load: updating savepoint!')
        pickle.dump(obj, open(p, "wb"))
        r_obj = obj
    else:
        print('save_load: save / load both empty.')
        r_obj = False
    return r_obj


def chkdir(d: str):
    """
    Create a directory if not existing.
    """
    if not os.path.isdir(d):
        os.makedirs(d)
        print(f'Created Directory: {d}')
    return