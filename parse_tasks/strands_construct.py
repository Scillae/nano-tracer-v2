from infra import StrandConstructor
from utils import SL_result_cache, generate_path


def strands_construct(data):
    """
    Construct nanostar series.
    :workspace data: (For result cache to work properly) Should indicate data['arm_number'] for arm number, data['temp'] for temperature, data['conc'] for salt concentration.
        May include data['flag_suffix'] for simulation configuration flags, data['conf_suffix'] for nanostar topology suffix, data['sp_suffix'] for special suffix.
    """
    varname = 'strands'
    # load data from result_cache
    SL_result_cache(data, varname, 'load')
    if data['SL_content']:
        cached_result = data['SL_content']
        data['SL_content'] = None
        return cached_result

    # not cached
    sc = StrandConstructor(generate_path(data, 'top_path'), generate_path(data, 'traj_path'))
    strands_series = sc.read_data()

    # save in result_cache
    data['SL_content'] = strands_series
    SL_result_cache(data, varname, 'save')
    return strands_series
