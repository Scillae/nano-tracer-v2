from parse_tasks import ns_construct
import numpy as np
from infra import TimeSeries
from utils import SL_result_catch


def patch_angle_calc(data: dict):
    """
    Calculate Inter-Arm Angles (Patch Angles) of a nanostar series.
    :workspace data: Must indicate data['arm_number'], data['ns_dims'];
        Should indicate data['temp'] for temperature, data['conc'] for salt concentration.
        May include data['flag_suffix'] for simulation configuration flags, data['conf_suffix'] for nanostar topology suffix, data['sp_suffix'] for special suffix.
    """
    varname = 'pa'
    # load data from result_catch
    SL_result_catch(data, varname, 'load')
    if data['SL_content']:
        catched_result = data['SL_content']
        data['SL_content'] = None
        return catched_result

    # not catched
    # load ns series
    ns_series = ns_construct(data)

    p_angs_vtime_dic = TimeSeries()  # {t_stamp: angle_results}

    ns_dims = data['ns_dims']
    for t_stamp, ns in ns_series.items():
        # Change of definition: origin~arm_end --> arm_start~arm_end
        # junction_pos = np.average(np.array([base.position for base in ns.center.values()]), 0)
        # assuming mass of nucs are the same.
        arms = ns.arms
        arms_idx = list(arms.keys())
        angle_results = {}  # {(ia1, ia2):(ang, (vec_1, vec_2), is_sharing)}

        # select two arms to form an arm pair
        for idx_1, ia1 in enumerate(arms_idx):
            arm1 = arms[ia1]
            for idx_2 in range(idx_1 + 1, len(arms_idx)):
                ia2 = arms_idx[idx_2]
                arm2 = arms[ia2]

                # select bps forming the vectors representing the direction of arms
                if ns_dims[0] > 12:
                    first_pair_a1 = arm1.base_pairs[ns_dims[0] - 2 - 10]
                    first_pair_a2 = arm2.base_pairs[ns_dims[0] - 2 - 10]
                    last_pair_a1 = arm1.base_pairs[ns_dims[0] - 2]
                    last_pair_a2 = arm2.base_pairs[ns_dims[0] - 2]
                else:  # arm too short!!
                    first_pair_a1 = arm1.base_pairs[ns_dims[0] - 10]  # b-form DNA, loop length == 10
                    first_pair_a2 = arm2.base_pairs[ns_dims[0] - 10]
                    last_pair_a1 = arm1.base_pairs[ns_dims[0]]
                    last_pair_a2 = arm2.base_pairs[ns_dims[0]]

                    # check if the two arms share a strand
                is_sharing_strand = True if len(
                    set([arm1.strand_id_0, arm1.strand_id_1, arm2.strand_id_0, arm2.strand_id_1])) == 3 else False

                # calculate patch angle
                vec_1 = CoM_bp(last_pair_a1) - CoM_bp(first_pair_a1)
                vec_2 = CoM_bp(last_pair_a2) - CoM_bp(first_pair_a2)
                ang_cos = obtain_cos(vec_1, vec_2)  # 0~180
                # ang_cross = obtain_cross(vec_1, vec_2) # -90~90
                ang = ang_cos  # if ang_cross >= 0 else (360 - ang_cos) # expanding angle range from 0~180 to 0~360
                angle_results[(ia1, ia2)] = (ang, (vec_1, vec_2), is_sharing_strand)
        p_angs_vtime_dic[t_stamp] = angle_results
    p_angs_vtime_dic.params['Arm_Pairs'] = [(ia1, arms_idx[idx_2]) for idx_1, ia1 in enumerate(arms_idx) for idx_2 in
                                            range(idx_1 + 1, len(arms_idx))]

    # save in result_catch
    data['SL_content'] = p_angs_vtime_dic
    SL_result_catch(data, varname, 'save')
    return p_angs_vtime_dic


def obtain_cos(v1, v2):
    return np.degrees(np.arccos(np.sum(v1 * v2) / (np.sqrt(np.sum(np.square(v1))) * np.sqrt(
        np.sum(np.square(v2))))))  # cos(x) = n1 * n2 / (|n1|*|n2|), angle only 0~pi, negative -> abs


def obtain_cross(v1, v2):
    return np.cross(v1, v2) / (np.sqrt(np.sum(np.square(v1))) * np.sqrt(
        np.sum(np.square(v2))))  # sin(x) = (n1 x n2) / (|n1|*|n2|), angle only -90~90


def CoM_bp(bp):
    # assuming mass of nucs are the same.
    return (np.array(bp[0].position) + np.array(bp[1].position)) / 2
