from utils.tools import get_ns_params, auto_generate_summary_range

stacking_options_ls = [('stacking', 'nonstacking', 'unlinked'),('no-stacking,unlinked','no-stacking,linked')] # 'stacking', 'nonstacking', 'unlinked', 'no-stacking,unlinked','no-stacking,linked', transiting?

def ns_pa_summary_plot(data): # entry point for generating the overall patch angle summary plot of selected nanostars
    '''
    Plot the nanostars' means vs stds under different conditions
    :data: Expected to contain:
        data['temp_list']; data['conc_list']; data['arm_number']; data['flag_suffix']; 
        data['conf_suffix']; data['sp_suffix']; data['ns_dims']; data['figure_savepath']
    '''
    # for each nanostar to be summarized, produce a separate data dictionary
    data_package_list = auto_generate_summary_range(data)
    # iterate through all the individual nanostar
    from joblib import Parallel, delayed # in Parallel
    result_data_list = Parallel(n_jobs=12)(delayed(ns_pa_plot)(child_data) for child_data in data_package_list)

    # plot the result
    plt = stacking_scatter_plot(result_data_list)
    
    # save the plot
    plt.savefig(data['figure_savepath'], dpi = 800)
    plt.close()
    return

def ns_pa_hist_plot(data): # entry point for generating the patch angle histogram plot of single nanostar
    '''
    Calculate the mean and the std (& skewness) of the patch angle.
    If wishing to plot the histogram, consider defining a ns_plot function in utils.tools module, and plot using 'n' in # moments
    '''
    # retrieve result from calc_func.patch_angle_calc
    varname = 'pa'
    # load data from result_cache
    SL_result_cache(data, varname, 'load')
    if not data['SL_content']:
        print(f'The results trying to be plotted is not cached. Run calc_func.patch_angle_calc first.')
        return False
    var_vals = data_process_func(data['SL_content'], data)
    data['SL_content'] = None

    # Plotting Parameter
    x_lim = (0,180)
    y_lim = (0,0.2)
    bin_num = 36

    # get ready for plotting histogram
    n,bin_edges = np.histogram(var_ls,bins = bin_num, range = x_lim)
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

    # moments
    n, m1, std, m3_s = moments_calc(n, var_ls)
    data['Summary'] = [m1, std]
    
    return

def ns_pa_plot(data):
    '''
    Calculate the mean and the std (& skewness) of the patch angle.
    If wishing to plot, consider defining a ns_plot function in utils.tools module
    '''
    # retrieve result from calc_func.stacking_local_identify_calc
    varname = 'sijpa'
    # load data from result_cache
    SL_result_cache(data, varname, 'load')
    if not data['SL_content']:
        print(f'The results trying to be plotted is not cached. Run plot_tasks.ns_pa_plot_vstime first.')
        return False
    stack_info = data_process_func(data['SL_content'], data)
    data['SL_content'] = None
    
    _,_,_,_,_,_,_,ns_struc = get_ns_params(data['arm_number'])

    # determine propensity to stack
    stk_cnt = 0
    for idx in ns_struc['PA']:
        stk_cnt += np.sum(np.array(stack_info[idx]['bool'],dtype=int)) # get the number of stack
    prop_stacking = stk_cnt/((ns_struc['arm_number']//2) * len(stack_info[(0,1)]['t'])) # divided by the maximum possible number of stack
    summ_dic = {'prop_stacking':prop_stacking} 
    for stacking_option in stacking_options_ls[0]: # angles in nanostars with stacked arm pairs
        var_ls = create_var_ls(stacking_option,ns_struc,pa_vtime_dic,stack_info)
        summ_dic[stacking_option] = hist_summ(var_ls, 36) if var_ls is not None else None # hist_summ returns n, m1, std, m3_s of the provided distribution
    for stacking_option in stacking_options_ls[1]: # angles in nanostars without stacked arm pairs
        var_ls = create_var_ls(stacking_option,ns_struc,pa_vtime_dic,stack_info)
        summ_dic[stacking_option] = hist_summ(var_ls, 36) if var_ls is not None else None # hist_summ returns n, m1, std, m3_s of the provided distribution
    return summ_dic


def create_var_ls(stacking_option, ns_struc, stack_info):
    '''
    Mask off the undesired angles according to the chosen stacking option.
    '''
    var_ls = []
    mask = np.zeros(len(stack_info[(0,1)]['t'])).astype(bool) # slots of undesired angles corresponds to 0
    if stacking_option in stacking_options_ls[0]: # angles in nanostars with stacked arm pairs
        for idx in ns_struc['PA']:
            mask += np.array(stack_info[idx]['bool'],dtype=bool)
        if np.sum(mask.astype(int)) < 100: # empirical!
            return None
        if stacking_option == 'unlinked':
            for idx in ns_struc['unlinked_PA']:
                var_ls.extend(stack_info[idx]['ang'] * mask.astype(int)) # assuming patch angle cannot be 0
        elif stacking_option == 'stacking':
            for idx in ns_struc['linked_PA']:
                var_ls.extend(stack_info[idx]['ang'] * mask.astype(int) * np.array(stack_info[idx]['bool']).astype(int)) # pool down the multiple lists
        elif stacking_option == 'nonstacking':
            for idx in ns_struc['linked_PA']:
                var_ls.extend(stack_info[idx]['ang'] * mask.astype(int) * (~np.array(stack_info[idx]['bool'])).astype(int))
        else:
            raise Exception('Corresponding angle not coded.')
    elif stacking_option in stacking_options_ls[1]: # angles in nanostars without stacked arm pairs
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
        else:
            raise Exception('Corresponding angle not coded.')
    else:
        raise Exception('Corresponding angle not coded.')
    var_ls = list(filter(None, var_ls)) # truncate 0s
    if len(var_ls) == 0:
        return None
    return var_ls # the pooled-down list of all desired angles in a nanostar's all time steps


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

def hist_summ(var_ls, bin_num):
    '''
    For patch angles.
    Calculate the 0th raw, 1st raw, 2nd central, and 3rd standardized moment of a given distribution.
    n can be used for plotting histograms if desired.
    '''
    x_lim = (0,180)
    y_lim = (0,0.2)
    n,bin_edges = np.histogram(var_ls,bins = bin_num, range = x_lim) # unused. for compatibility only
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

def stacking_scatter_plot(plot_summ_dic): # legacy code
    # f=open('tmp\\tmp.txt','a')
    row,col = (2,1)
    fig = plt.figure(figsize=(3*col+1,3.5*row))
    gs = fig.add_gridspec(row, col, hspace=0.3, wspace=0.1) # hspace=0.4, wspace=0.1
    axs = gs.subplots(sharex='col')
    # fig, axs = plt.subplots()
    jun_list = [0, 1, 2, 5, 10]
    conc_list = [0.05, 0.1, 0.3, 0.5] # 0.05, 0.1, 0.3, 0.5
    temp_list = [20, 23, 27, 30, 40, 50] # 20, 23, 27, 30, 40, 50
    markers = ['v', 's',  '*', 'o', 'p', '8', 'h', 'H', 'D', 'd', 'P', 'X']
    plotting_stacking_type = {'with_stacking':('stacking', 'nonstacking', 'unlinked'), 'linked': ('stacking', 'nonstacking', 'no-stacking,linked')}
    # different stacking types?
    cm_sub = np.linspace(0.2, 0.7, 6)
    cmap = plt.get_cmap('hot')
    colors = [cmap(x) for x in cm_sub]
    for cond, vals in plot_summ_dic.items():
        arm,jun,temp,conc = cond
        prop_stacking = vals['prop_stacking']
        del vals['prop_stacking']
        for stacking_type, summ in vals.items():
            if summ == None:
                continue
            freq, m1, std, m3_s = summ
            n = sum(freq)
            if jun == 2:
                conf_suffix = ''
            else:
                conf_suffix = f'-jun_{jun}'
            sp_suffix = ''
            # label = f'{arm}arms@({temp}C,{conc}M){conf_suffix}{sp_suffix}'
            # f.write(f'{label} , Type:{stacking_type} ~ Mean:{m1:.3f} ; STD:{std:.3f} ; #total:{n} ; is_stacking: {prop_stacking}\n')
            if stacking_type in plotting_stacking_type['with_stacking']:
                lw = 0.5
                # title = f'{arm}Arms, Mean and STD of Patch Angles when NS is stacking'
                if stacking_type in plotting_stacking_type['linked']:
                    axs[0].scatter(m1, std, color=colors[temp_list.index(temp)], marker=markers[conc_list.index(conc)], s=5*(jun+1), linewidth=lw, edgecolor='#000000')
                else:
                    axs[0].scatter(m1, std, color=colors[temp_list.index(temp)], marker=markers[conc_list.index(conc)], s=5*(jun+1), linewidth=lw, edgecolor=colors[temp_list.index(temp)], facecolors='none')
            else:
                lw = 0.5 # 0.2
                # title = f'{arm}Arms, Mean and STD of Patch Angles when NS is NOT stacking'
                if stacking_type in plotting_stacking_type['linked']:
                    axs[1].scatter(m1, std, color=colors[temp_list.index(temp)], marker=markers[conc_list.index(conc)], s=5*(jun+1), linewidth=lw, edgecolor='#000000')
                else:
                    axs[1].scatter(m1, std, color=colors[temp_list.index(temp)], marker=markers[conc_list.index(conc)], s=5*(jun+1), linewidth=lw, edgecolor=colors[temp_list.index(temp)], facecolors='none')
            # if stacking_type in plotting_stacking_type['linked']:
            #     ax.scatter(m1, std, color=colors[temp_list.index(temp)], marker=markers[conc_list.index(conc)], s=5*(jun+1), linewidth=lw, edgecolor='#000000')
            # else:
            #     ax.scatter(m1, std, color=colors[temp_list.index(temp)], marker=markers[conc_list.index(conc)], s=5*(jun+1), linewidth=lw, edgecolor=colors[temp_list.index(temp)], facecolors='none')
    axs[0].set_xlabel('Mean Patch Angle', fontsize=14) # r'$\mu(^\circ)$'
    axs[0].set_ylabel('Patch Angle Standard Deviation', fontsize=8) # r'$\sigma(^\circ)$'
    axs[0].xaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter(u"{x:.0f}°"))
    axs[0].yaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter(u"{x:.0f}°"))
    axs[0].tick_params(bottom=True,top=True,left=True,right=True,direction='in')
    title = f'{arm}Arms-oxrna2, Mean and STD of Patch Angles when NS is stacking'
    axs[0].set_title(title, fontsize=8)
    # axs[0].set_xlim((70,160))
    # axs[0].set_ylim((15,50))
    axs[0].set_xlim((80,150))
    axs[0].set_ylim((15,45))

    axs[1].set_xlabel('Mean Patch Angle', fontsize=14) # r'$\mu(^\circ)$'
    axs[1].set_ylabel('Patch Angle Standard Deviation', fontsize=8) # r'$\sigma(^\circ)$'
    axs[1].xaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter(u"{x:.0f}°"))
    axs[1].yaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter(u"{x:.0f}°"))
    axs[1].tick_params(bottom=True,top=True,left=True,right=True,direction='in')
    title = f'{arm}Arms-oxrna2, Mean and STD of Patch Angles when NS is NOT stacking'
    axs[1].set_title(title, fontsize=8)
    # axs[1].set_xlim((70,160))
    # axs[1].set_ylim((15,50))
    axs[1].set_xlim((80,150))
    axs[1].set_ylim((15,45))
    
    return plt