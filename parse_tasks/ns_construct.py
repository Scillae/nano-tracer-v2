from utils.result_catch import generate_path
from infra.NanoConstructor import NanoConstructor
from parse_tasks.strands_construct import strands_construct

def ns_construct(data):
    varname = 'ns'
    # load data from result_catch
    SL_result_catch(data, varname, 'load')
    if data['SL_content']:
        catched_result = data['SL_content']
        data['SL_content'] = None
        return catched_result
    

    # not catched
    # load strands
    strands_series = strands_construct(data, generate_path(data, 'top_path'), generate_path(data, 'traj_path'))


    nc = NanoConstructor(strands_series, data['ns_dims'], data['arm_number'])
    # box_dim hacking
    import re
    with open(path_traj,'r') as f:
        f.readline()
        ret=re.match('^b = ([0-9]+) ([0-9]+) ([0-9]+)\n',f.readline())
    box_dim = np.array((int(ret.group(1)),int(ret.group(2)),int(ret.group(3))))
    # hacking ends

    ns_series = nc.construct(box_dim=box_dim)
    # save in result_catch
    data['SL_content'] = ns_series
    SL_result_catch(data, varname, 'save')
    return ns_series