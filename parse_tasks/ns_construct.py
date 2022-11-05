from utils import generate_path, SL_result_cache
from infra import NanoConstructor
from parse_tasks.strands_construct import strands_construct
import numpy as np


def ns_construct(data):
    """
    Construct nanostar series.
    data: Must indicate data['arm_number'], data['ns_dims'];
        (For result cache to work properly) Should indicate data['temp'] for temperature, data['conc'] for salt concentration.
        May include data['flag_suffix'] for simulation configuration flags, data['conf_suffix'] for nanostar topology suffix, data['sp_suffix'] for special suffix.
    """
    varname = 'ns'
    # load data from result_cache
    SL_result_cache(data, varname, 'load')
    if data['SL_content']:
        cached_result = data['SL_content']
        data['SL_content'] = None
        return cached_result

    # not cached
    # load strands
    strands_series = strands_construct(data)

    nc = NanoConstructor(strands_series, data['ns_dims'], data['arm_number'])

    ns_series = nc.construct(box_dim=strands_series.params['box_dim'])
    # save in result_cache
    data['SL_content'] = ns_series
    SL_result_cache(data, varname, 'save')
    return ns_series
