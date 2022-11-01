from utils.tools import get_ns_params

def ns_pa_plot_summary(data):
    '''
    Plot the nanostars' means vs stds under different conditions
    :data: Expected to contain:
        data['temp_list']; data['conc_list']; data['arm_number']; data['flag_suffix']; 
        data['conf_suffix']; data['sp_suffix']; data['ns_dims']; data['varname']; 
    '''
    # for each nanostar to be summarized, produce a separate data dictionary
    data_package_list = auto_generate_summary_range(data)
    # iterate through all the individual nanostar
    result_data_list = []
    for child_data in data_package_list: # can be parallel if needed, just make ns_pa_plot return data
        result_data_list.append(ns_pa_plot(child_data)) # add summary into child_data
    
    # plot the result

    
    return

def ns_pa_plot(data):
    '''
    Calculate the mean and the std (& skewness) of the patch angle.
    If wishing to plot, consider defining a ns_plot function in utils.tools module
    '''
    # retrieve result from calc_func.stacking_local_identify_calc
    varname = 'sijpa'
    # load data from result_catch
    SL_result_catch(data, varname, 'load')
    if not data['SL_content']:
        print(f'The results trying to be plotted is not catched. Run plot_tasks.ns_pa_plot_vstime first.')
        return False
    stack_info = data_process_func(data['SL_content'], data)
    data['SL_content'] = None
    
    
    stacking_options_ls = [('stacking', 'nonstacking', 'unlinked'),('no-stacking,unlinked','no-stacking,linked')] # 'stacking', 'nonstacking', 'unlinked', 'no-stacking,unlinked','no-stacking,linked', transiting?
    _,_,_,_,_,_,_,ns_struc = get_ns_params(data['arm_number'])

    # determine propensity to stack
    stk_cnt = 0
    for idx in ns_struc['PA']:
        stk_cnt += np.sum(np.array(stack_info[idx]['bool'],dtype=int))
    prop_stacking = stk_cnt/((ns_struc['arm_number']//2) * len(stack_info[(0,1)]['t']))
    summ_dic = {'prop_stacking':prop_stacking} 
    for stacking_option in stacking_options_ls[0]:
        var_ls = create_var_ls(stacking_option,ns_struc,pa_vtime_dic,stack_info)
        summ_dic[stacking_option] = hist_summ(var_ls, 36) if var_ls is not None else None # n, m1, std, m3_s
    for stacking_option in stacking_options_ls[1]:
        var_ls = create_var_ls(stacking_option,ns_struc,pa_vtime_dic,stack_info)
        summ_dic[stacking_option] = hist_summ(var_ls, 36) if var_ls is not None else None # n, m1, std, m3_s
    return summ_dic


def create_var_ls(stacking_option, ns_struc, stack_info):
    var_ls = []
    mask = np.zeros(len(stack_info[(0,1)]['t'])).astype(bool)
    if stacking_option in ['unlinked', 'stacking', 'nonstacking']: # w/ stk
        for idx in ns_struc['PA']:
            mask += np.array(stack_info[idx]['bool'],dtype=bool)
        if np.sum(mask.astype(int)) < 100: # empirical!
            return None
        if stacking_option == 'unlinked':
            for idx in ns_struc['unlinked_PA']:
                var_ls.extend(stack_info[idx]['ang'] * mask.astype(int)) # assuming patch angle cannot be 0
        elif stacking_option == 'stacking':
            for idx in ns_struc['linked_PA']:
                var_ls.extend(stack_info[idx]['ang'] * mask.astype(int) * np.array(stack_info[idx]['bool']).astype(int))
        elif stacking_option == 'nonstacking':
            for idx in ns_struc['linked_PA']:
                var_ls.extend(stack_info[idx]['ang'] * mask.astype(int) * (~np.array(stack_info[idx]['bool'])).astype(int))
    elif stacking_option in ['no-stacking,linked','no-stacking,unlinked']: # w/o stk
        for idx in ns_struc['PA']:
            mask += np.array(stack_info[idx]['bool'],dtype=bool)
        mask = ~mask
        if np.sum(mask.astype(int)) < 100: # empirical!
            return None
        if stacking_option == 'no-stacking,linked':
            for idx in ns_struc['linked_PA']:
                var_ls.extend(stack_info[idx]['ang'] * mask.astype(int))
        elif stacking_option == 'no-stacking,unlinked':
            for idx in ns_struc['unlinked_PA']:
                var_ls.extend(stack_info[idx]['ang'] * mask.astype(int))
    var_ls = list(filter(None, var_ls)) # truncate 0s
    if len(var_ls) == 0:
        return None
    return var_ls


def data_process_func(p_ang_res, data):
    '''
    Convert stored results into what suitable for plotting.
    :p_ang_res: patch angle result
    :data: expected to have 'Arm_Pairs', 'Arm_IDs'
    '''
    import plot_tasks.ns_pa_plot_vstime.data_process_func
    angle_dic = plot_tasks.ns_pa_plot_vstime.data_process_func(p_ang_res, data)
    # pool down the data into list
    t_ls = [t for t in angle_dic[(0,1)] if type(t) == int]
    ang_ls = [t_dic[t][0] for t_dic in angle_dic.values() if t_dic['is_sharing']==True for t in t_ls] # pooled for is_sharing
    return ang_ls

def auto_generate_summary_range(data):
    data_package_list = []
    for arm_num in data['arm_number']:
        for conc in data['conc_list']:
            for temp in data['temp_list']:
                child_data = {}
                child_data['arm_number'] = arm_num
                child_data['temp'] = temp
                child_data['conc'] = conc
                child_data['flag_suffix'] = data['flag_suffix']
                child_data['conf_suffix'] = data['conf_suffix']
                child_data['sp_suffix'] = data['sp_suffix']
                child_data['ns_dims'] = data['ns_dims']
                child_data['varname'] = data['varname']
                data_package_list.append(child_data)
    return data_package_list

def hist_summ(var_ls, bin_num):
    x_lim = (0,180)
    y_lim = (0,0.2)
    n,bin_edges = np.histogram(var_ls,bins = bin_num, range = x_lim)
    n, m1, std, m3_s = moments_calc(n, var_ls)
    return n, m1, std, m3_s

def moments_calc(n, var_ls):
    '''
    Calculate the 0th raw, 1st raw, 2nd central, and 3rd standardized moment of a given distribution.
    '''
    n = np.array(n)
    m0 = np.sum(n) # 0th unitless raw moment: integration
    m1 = np.sum(var_ls)/m0
    m2_c = stats.moment(var_ls, moment=2) # 2nd central moment: variance
    std = m2_c ** 0.5 # standard deviation
    m3_c = stats.moment(var_ls, moment=3) # 3rd central moment
    m3_s = m3_c / (std ** 3) # 3rd standardized moment: skewness
    return n, m1, std, m3_s

def ns_pa_hist_plot(data):
    '''
    Calculate the mean and the std (& skewness) of the patch angle.
    If wishing to plot the histogram, consider defining a ns_plot function in utils.tools module
    '''
    # retrieve result from calc_func.patch_angle_calc
    varname = 'pa'
    # load data from result_catch
    SL_result_catch(data, varname, 'load')
    if not data['SL_content']:
        print(f'The results trying to be plotted is not catched. Run calc_func.patch_angle_calc first.')
        return False
    var_vals = data_process_func(data['SL_content'], data)
    data['SL_content'] = None

    # get ready for plotting histogram
    n,bin_edges = np.histogram(var_ls,bins = bin_num, range = x_lim)
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

    # moments
    n, m1, std, m3_s = moments_calc(n, var_ls)
    data['Summary'] = [m1, std]
    
    return
