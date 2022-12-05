from parse_tasks import ns_construct
import numpy as np
from infra import TimeSeries
from utils.result_cache import SL_result_cache

def k2_calc(data:dict):
    '''

    '''
    varname = 'k2'
    # load data from result_cache
    SL_result_cache(data, varname, 'load')
    if data['SL_content']:
        cached_result = data['SL_content']
        data['SL_content'] = None
        return cached_result

    # not cached
    # load ns series
    ns_series = ns_construct(data)

    # construct result series of k2
    k2_results = TimeSeries()  # {t_stamp: (k2,(l1, l2, l3),(v1, v2, v3))}

    for t_stamp, ns in ns_series.items():
        i_arr33 = np.zeros((3,3)) # the density matrix
        CoM_pos = np.zeros(3) # Center of Mass of the *whole nanostar*
        base_cnt = 0

        # obtain CoM
        for strand in ns.strands.values():
            for base in strand.base_sequence.values():
                CoM_pos = np.add(CoM_pos, np.array(base.position))
                base_cnt += 1
        CoM_pos = np.divide(CoM_pos, base_cnt)
        
        # calculate I, the density matrix
        for strand in ns.strands.values():
            for base in strand.base_sequence.values():
                m = base.mass
                x, y, z = base.position - CoM_pos
                i_arr33[0][0] += (y**2 + z**2)*m
                i_arr33[1][1] += (x**2 + z**2)*m
                i_arr33[2][2] += (y**2 + x**2)*m
                i_arr33[0][1] += -x*y*m
                i_arr33[1][0] += -x*y*m
                i_arr33[0][2] += -x*z*m
                i_arr33[2][0] += -x*z*m
                i_arr33[1][2] += -z*y*m
                i_arr33[2][1] += -z*y*m
        e_vals, e_vecs = np.linalg.eig(i_arr33)
        l1, l2, l3 = e_vals
        v1, v2, v3 = e_vecs
        k2 = 1 - (27*l1*l2*l3)/((l1+l2+l3)**3)
        k2_results[t_stamp] = (k2,(l1, l2, l3),(v1, v2, v3))
    
    # save in result_cache
    data['SL_content'] = k2_results
    SL_result_cache(data, varname, 'save')
    return k2_results