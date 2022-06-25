from parse_tasks import ns_construct
import numpy as np
from collections import OrderedDict

def patch_angle_calc(data, arm_number, ns_dims):
    varname = 'pa'
    # load data from result_catch
    SL_result_catch(data, varname, 'load')
    if data['SL_content']:
        catched_result = data['SL_content']
        data['SL_content'] = None
        return catched_result


    # not catched
    # load ns series
    ns_series = ns_construct(data, data['arm_number'], ns_dims)


    p_angs_vtime_dic = OrderedDict() #{t_stamp: angle_results}
    for t_stamp, ns in ns_series.time_capsule.items():
        # Change of definition: origin~arm_end --> arm_start~arm_end
        # junction_pos = np.average(np.array([base.position for base in ns.center.values()]), 0) # assuming mass of nucs are the same.
        arms = ns.arms
        arms_idx = list(arms.keys())
        is_sharing_strand = True
        angle_results = {}   # {(ia1, ia2):(ang, (vec_1, vec_2), is_sharing)}
        for idx_1, ia1 in enumerate(arms_idx):
            arm1 = arms[ia1]
            for idx_2 in range(idx_1+1, len(arms_idx)):
                ia2 = arms_idx[idx_2]
                arm2 = arms[ia2]
                if dims_ls[0] > 12:
                    first_pair_a1 = arm1.base_pairs[dims_ls[0]-2-10]
                    first_pair_a2 = arm2.base_pairs[dims_ls[0]-2-10]
                    last_pair_a1 = arm1.base_pairs[dims_ls[0]-2]
                    last_pair_a2 = arm2.base_pairs[dims_ls[0]-2]
                else: # arm too short!!
                    first_pair_a1 = arm1.base_pairs[dims_ls[0]-10] # b-form DNA, loop length == 10
                    first_pair_a2 = arm2.base_pairs[dims_ls[0]-10]
                    last_pair_a1 = arm1.base_pairs[dims_ls[0]]
                    last_pair_a2 = arm2.base_pairs[dims_ls[0]]                    
                base_ls = list(last_pair_a1)
                base_ls.extend(last_pair_a2)
                if len(set([base.strand_id for base in base_ls])) < 4:
                    is_sharing_strand = True
                else:
                    is_sharing_strand = False
                fp_a1_pos = np.average(np.array([base.position for base in first_pair_a1]),0)
                fp_a2_pos = np.average(np.array([base.position for base in first_pair_a2]),0)
                lp_a1_pos = np.average(np.array([base.position for base in last_pair_a1]),0) # assuming mass of nucs are the same.
                lp_a2_pos = np.average(np.array([base.position for base in last_pair_a2]),0)
                # vec_1 = lp_a1_pos - junction_pos
                # vec_2 = lp_a2_pos - junction_pos
                vec_1 = lp_a1_pos - fp_a1_pos
                vec_2 = lp_a2_pos - fp_a2_pos
                ang_cos = obtain_cos(vec_1, vec_2) # 0~180
                # ang_cross = obtain_cross(vec_1, vec_2) # -90~90
                ang = ang_cos # if ang_cross >= 0 else (360 - ang_cos) # expanding angle range from 0~180 to 0~360
                angle_results[(ia1, ia2)] = (ang, (vec_1, vec_2), is_sharing_strand)
        p_angs_vtime_dic[t_stamp] = angle_results
    
    # save in result_catch
    data['SL_content'] = p_angs_vtime_dic
    SL_result_catch(data, varname, 'save')
    return p_angs_vtime_dic


def obtain_cos(v1,v2):
    return np.degrees(np.arccos(np.sum(v1*v2)/(np.sqrt(np.sum(np.square(v1))) * np.sqrt(np.sum(np.square(v2)))))) # cos(x) = n1 * n2 / (|n1|*|n2|), angle only 0~pi, negative -> abs

def obtain_cross(v1,v2):
    return np.cross(v1,v2)/(np.sqrt(np.sum(np.square(v1))) * np.sqrt(np.sum(np.square(v2)))) # sin(x) = (n1 x n2) / (|n1|*|n2|), angle only -90~90

