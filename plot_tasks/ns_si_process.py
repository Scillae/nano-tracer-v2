import numpy as np

from plot_tasks.ns_pa_plot_vstime import create_ia_idx_ls


def data_process_func(stacking_res: dict, data: dict):
    '''
    Obtain information regarding stacked arm pairs vs time. (Generate stack_info)
    (If stacked info is to be plotted, use this as the data_process_func for that script.)
    :stacking_res: result of stack information
    :data: no dependence yet. remember to register if added later.
    '''
    # stacking_res: OrderedDict{t_stamp: is_stacking, stack_ls}
    arms, temp, conc, sp_suffix, conf_suffix, flag_suffix, ns_dims = data
    # init
    stack_info = {}
    for ia1,ia2 in create_ia_idx_ls(arms):
        stack_info[(ia1,ia2)] = {}
        stack_info[(ia1,ia2)]['bool'] = np.zeros(len(stacking_res.keys()),dtype=bool)
        stack_info[(ia1,ia2)]['t'] = np.zeros(len(stacking_res.keys()),dtype=int)
        stack_info[(ia1,ia2)]['adj_bps'] = np.zeros(len(stacking_res.keys()),dtype=bool)
    # fill data in
    for i, (t, (is_stacking, stack_ls)) in enumerate(stacking_res.items()):
        
        if is_stacking == False:
            continue
        for adj_bp, adj_bp2, arm_id, arm2_id in stack_ls:
            stack_info[(arm_id,arm2_id)]['bool'][i] = True
            stack_info[(arm_id,arm2_id)]['adj_bps'][i] = adj_bp, adj_bp2 # deprecated
    # smoothing
    window_hw = 5 # window_hw *2 +1 == window_width
    for ia1,ia2 in create_ia_idx_ls(arms):
        bool_ls = stack_info[(ia1,ia2)]['bool']
        smooth_ls = [True if sum(bool_ls[i:i+window_hw*2 +1]) > window_hw+1 else False for i in range(len(bool_ls)-window_hw*2)] # threshold: 50%
        r_ls = [True if sum(bool_ls[0:i+window_hw +1]) > window_hw or bool_ls[i] else False for i in range(window_hw)] # still 50% of full window width?
        r_ls.extend(smooth_ls)
        r_ls.extend([True if sum(bool_ls[i-window_hw:]) > window_hw or bool_ls[i] else False for i in range(-window_hw,0)])
        assert len(r_ls) == len(bool_ls)
        stack_info[(ia1,ia2)]['bool'] = r_ls
        stack_info[(ia1,ia2)]['t'] = list(stacking_res) # filling in time indices as well
    return stack_info
