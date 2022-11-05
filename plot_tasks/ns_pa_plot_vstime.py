import numpy as np
from utils.result_cache import generate_path
from utils.tools import chkdir, data_query, dims_adjust, get_ns_params
import plot_tasks.ns_si_process.data_process_func

def ns_pa_vstime_plot(data: dict):
    '''
    Plot nanostar's patch angles vs time per arm pair. Raw values, Running avgs, and running rmsds are plotted.
    Single trajectory!
    WIP: what's expected in data
    '''
    varname = 'pa'
    # remember to truncate unused unpacked values later.
    arm_number, temp, conc, sp_suffix, conf_suffix, flag_suffix, dims_ls = data_query(data, ['arm_number', 'temp', 'conc', 'sp_suffix', 'conf_suffix', 'flag_suffix', 'dims_ls'], ['int','double','double','str','str','str','ls'])
    dims_adjust(dims_ls, conf_suffix, single, sp_suffix) # ns_dims: [arm, center, single_end] #lengths

    # retrieve result from calc_func.stacking_local_identify_calc
    varname = 'si'
    # load data from result_cache
    SL_result_cache(data, varname, 'load')
    if not data['SL_content']:
        raise Exception('The results trying to be plotted is not cached. Run calc_tasks.stacking_local_identify_calc first.')
    stack_info = plot_tasks.ns_si_process.data_process_func(data['SL_content'], data)
    data['SL_content'] = None

    # retrieve result from calc_func.patch_angle_calc
    varname = 'pa'
    # load data from result_cache
    SL_result_cache(data, varname, 'load')
    if not data['SL_content']:
        raise Exception('The results trying to be plotted is not cached. Run calc_func.patch_angle_calc first.')
    var_vals = data_process_func(data['SL_content'], data)
    data['SL_content'] = None

    #### plot confs ####
    x_var = rf'Patch Angles ($^\circ$)'
    plot_path = generate_path(data, 'plot_path')
    label = generate_path(data, 'label')
    #### conf ends ####

    plot_confs = {'x_var':x_var, 'plot_path':plot_path, 'label':label}
    ns_pa_plot(var_vals, stack_info, plot_confs, data) # move var_vals, stack_info into data if they might be accessed anywhere except this function
    # can register 'stack_info' into data if needed. ns_time_pa_plot has modified it.
    # plot saved inside the function
    
    # save in result_cache
    varname = 'sijpa' # stack_info joining patch angle
    data['SL_content'] = stack_info
    SL_result_cache(data, varname, 'save')
    return

def ns_pa_plot(var_vals: dict, stack_info: dict, plot_confs: dict, data: dict):
    '''
    Plot the value (patch angle vtime) vs. time plot of a single trajectory.
    :var_vals: value vs time.
    :stack_info: the dictionary that logged stack info
    :plot_confs: params of the plot.
    :data: in which the descriptions of nanostars (trajectory) are stored. NOT the data to be plotted.
    '''
    # unpack configurations & obtain parameters
    x_var, plot_path, label = data_query(data, ['x_var','plot_path','label'],['str','str','str'])
    time_window_width, stacking_min_length, stacking_crit_ang, stacking_crit_rmsd, _, _, _, ns_struc = get_ns_params(int(label[0]))
    
    # subplots layout
    row,col = (len((var_vals.items())),3)
    fig = plt.figure(figsize=(3*3*col+1,3.5*row))
    gs = fig.add_gridspec(row, col, hspace=0.3, wspace=0.1) # hspace=0.4, wspace=0.1
    axs = gs.subplots(sharex='col')
    
    # plot running avg & rmsd
    for i, (ia1,ia2), ang_vs_time in enumerate(var_vals.items()): # ang_vs_time: {time:ang}
        # i: loop of arm pairs.
        # colorcode the linked & unlinked
        if (ia1,ia2) in ns_struc['linked_PA']:
            legend_color='#F7AED0'# #E55050
        else:
            legend_color='#4994FF'
        
        ang_ls = [ang_vs_time[t][0] for t in ang_vs_time if type(t) == int] # there might be misc keys mixed in t. Need to fix.
        time_ls = [t for t in ang_vs_time if type(t) == int]
        assert len(stack_info[(ia1,ia2)]['t']) == len(time_ls) # frames may be dropped in calc_func.patch_angle_calc. If so, write lines that drop time frames in stack_info which is absent in patch angle result
        time_idx = time_ls[time_window_width//2:-time_window_width//2+1] # truncate because runningavg & runningrmsd have less datapoints
        
        # data preparation
        stack_info[(ia1,ia2)]['avg'] = [np.sum(ang_ls[j:j+time_window_width])/time_window_width for j in range(len(time_idx))]
        stack_info[(ia1,ia2)]['rmsd'] = [np.std(ang_ls[j:j+time_window_width]) for j in range(len(time_idx))]
        stack_info[(ia1,ia2)]['raw'] = ang_ls[time_window_width//2:-time_window_width//2+1] # short length, use with running avg & rmsd
        stack_info[(ia1,ia2)]['ang'] = ang_ls # full length

        avg_running_is_stacking = stack_info[(ia1,ia2)]['bool'][5:-5] # truncated because running avg & rmsd have less datapoints
        ang_running_avg = stack_info[(ia1,ia2)]['avg'] 
        ang_running_rmsd = stack_info[(ia1,ia2)]['rmsd']
        
        # plotting
        c = ['#4994FF' if is_stacking == True else '#FF55FF' for is_stacking in avg_running_is_stacking]
        axs[i,0].scatter(time_idx,ang_running_avg, label = f'Arm{ia1}~Arm{ia2}', s=4, c=c)
        axs[i,1].scatter(time_idx,ang_running_rmsd, label = f'Arm{ia1}~Arm{ia2}', s=4, c=c)
        # the 'fancy' shaded running avg w/ rmsd overlay on the raw values
        axs[i,2].fill_between(time_idx,np.array(ang_running_avg)-np.array(ang_running_rmsd), np.array(ang_running_avg)+np.array(ang_running_rmsd), where=np.array(c)=='#4994FF', color='#4994FF', alpha=0.4, linewidth=0)
        axs[i,2].fill_between(time_idx,np.array(ang_running_avg)-np.array(ang_running_rmsd), np.array(ang_running_avg)+np.array(ang_running_rmsd), where=np.array(c)=='#FF55FF', color='#00FF00', alpha=0.2, linewidth=0)
        axs[i,2].plot(time_idx,ang_running_avg, label = f'Arm{ia1}~Arm{ia2}', c='#4994FF', linewidth=0.4, zorder=-100)

        # setting parameters
        axs[i,0].set_ylim((0,180))
        axs[i,0].yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(base=60))
        axs[i,0].yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(base=30))
        axs[i,0].set_ylabel(f'{x_var} (Avg, wid:{time_window_width})') #, fontsize=18
        
        l = axs[i,0].legend(loc = 'lower right')
        for text in l.get_texts():
            text.set_color(legend_color)
        
        axs[i,1].set_xlabel('Time') #, fontsize=18
        axs[i,1].set_ylabel(f'{x_var} (STD, wid:{time_window_width})') #, fontsize=18
        
        l = axs[i,1].legend(loc = 'lower right')
        for text in l.get_texts():
            text.set_color(legend_color)
    
    # plot raw patch angle values
    for i, ((ia1,ia2), ang_vs_time) in enumerate(var_vals.items()):
        ang_ls = [ang_vs_time[t][0] for t in ang_vs_time if type(t) == int]
        time_ls = [t for t in ang_vs_time if type(t) == int]
        axs[i,2].scatter(time_ls,ang_ls, label = f'Arm{ia1}~Arm{ia2}', s=4, c= ['#4994FF' if is_stacking == True else '#00FF00' for is_stacking in stack_info[(ia1,ia2)]['bool']]) # ,c=patch_color_list[i]
        
        axs[i,2].set_ylabel(x_var) #, fontsize=18
        axs[i,2].yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(base=60))
        axs[i,2].yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(base=30))
        axs[i,2].set_ylim((0,180))

    # briefing the ratios of interest in title
    total_count = counting_stacking(stack_info, 100, True, False, ns_struc)
    stacking_0_count = counting_stacking(stack_info, 0, False, False, ns_struc)
    stacking_not0_count = counting_stacking(stack_info, 0, False, True, ns_struc)
    stacking_1_count = counting_stacking(stack_info, 1, False, False, ns_struc)
    stacking_2_count = counting_stacking(stack_info, 2, False, False, ns_struc)
    stacking_3_count = counting_stacking(stack_info, 3, False, False, ns_struc)
    plt.suptitle(f'Total:{total_count}, 0_stk:{stacking_0_count}, Not0_stk:{stacking_not0_count}, 1_stk: {stacking_1_count}, 2_stk: {stacking_2_count}, 3_stk: {stacking_3_count}')
    
    chkdir(os.path.dirname(plot_path))
    plt.savefig(plot_path,dpi=400)
    plt.close()
    return

def counting_stacking(stack_info: dict, number_stack: int, is_count_all: bool, is_reverse_select_number_stack: bool, ns_struc: dict):
    '''
    Count in how many frames the ns contains the designated number of stack.
    :stack_info: the dictionary that logged stack info
    :number_stack: the designated number of stack
    :is_count_all: override #stack. count all frames
    :is_reverse_select_number_stack: count frames that have #stack not being the designated value
    :ns_struc: the dictionary that describe the topology layout of the nanostar's arms.
    '''
    if is_count_all:
        # is_containing_stack = np.array(stack_info[(0,1)]['bool']) + np.array(stack_info[(1,2)]['bool']) + np.array(stack_info[(2,3)]['bool']) + np.array(stack_info[(0,3)]['bool']) # selecting w/ stacking
        is_containing_stack_arr = np.ones(len(stack_info[(0,1)]['bool']),dtype=bool) # draw everything
    else:
        is_containing_stack_arr = np.zeros(len(stack_info[(0,1)]['bool']),dtype=int)
        for iaidx in ns_struc['PA']: # not 'linked_PA': no longer assume linked
            is_containing_stack_arr += np.array(stack_info[iaidx]['bool'],dtype=int)
        is_containing_stack_arr = is_containing_stack_arr == number_stack
    if not is_reverse_select_number_stack:
        return np.sum(is_containing_stack_arr)
    else:
        return np.sum(~is_containing_stack_arr)     

    stacking_result = {}
    # stacking in one series
    for iaidx in ns_struc['linked_PA']:
        s = stack_info[iaidx]['bool']
        stacking_result[iaidx] = [sum(s),len(s),sum(s)/len(s)]
    # simutaneous stacking
    if ns_struc['arm_number'] == 4:
        for iaidx1, iaidx2 in ns_struc['pairing_linked']:
            s1 = stack_info[iaidx1]['bool']
            s2 = stack_info[iaidx2]['bool']
            compare_result = [True if s1[i]==True and s2[i]==True else False for i in range(len(s1))]
            stacking_result[(iaidx1,iaidx2)] = [sum(compare_result),len(compare_result),sum(compare_result)/len(compare_result)]
    return stacking_result



def data_process_func(p_ang_res: dict, data: dict): # should be trimmed.
    '''
    Convert stored results into what suitable for plotting.
    :p_ang_res: patch angle result
    '''
    arm_pairs = p_ang_res.params['Arm_Pairs']
    arm_IDs = p_ang_res.params['Arm_IDs']

    # tracing of specific p-angles setup
    angle_dic = OrderedDict() # {(ia1, ia2):{t_stamp:angle}} 
    for ia1, ia2 in arm_pairs:
        angle_dic[(ia1, ia2)] = OrderedDict()
    
    drop_t_ls = [] # frames discarded
    for t_stamp, ang_res_per_armpair in p_ang_res.items():
        linked_cnt = 0
        for arm_pair, (ang, vec_tp, is_linked) in ang_res_per_armpair.items():
            angle_dic[arm_pair][t_stamp] = [ang, vec_tp] # collecting patch angles of a specific arm pair
            if 'is_linked' not in angle_dic[ia_tp].keys():
                angle_dic[ia_tp]['is_linked'] = is_linked # log if the arm pair is linked
            if is_linked: # sanity check: count the total linked arm pairs in each frame
                linked_cnt += 1
        # the number of linked arm pairs should equals to the number of arms
        if linked_cnt != len(arm_IDs): 
            drop_t_ls.append(t_stamp)
            angle_dic[arm_pair][t_stamp] = None
            print(f'Drop frame due to wrong number of linked patch angle: {t_stamp}')
    print(f'Total time steps dropped: {len(drop_t_ls)}')
    return angle_dic

def create_ia_idx_ls(arm_num: int):
    return [(ia1,ia2) for ia1 in range(arm_num) for ia2 in range(ia1+1,arm_num)]
