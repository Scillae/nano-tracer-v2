from utils import generate_path, SL_result_catch
from infra import NanoConstructor
from parse_tasks.strands_construct import strands_construct
import numpy as np


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
    
    ns_series = nc.construct(box_dim=box_dim)
    # save in result_catch
    data['SL_content'] = ns_series
    SL_result_catch(data, varname, 'save')
    return ns_series
