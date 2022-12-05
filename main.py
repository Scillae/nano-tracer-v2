

def generate_data_for_ns_pa_summary_plot():
    data = {
        'temp_list': [20, 23, 27, 30, 40, 50], 
        'conc_list': [0.05, 0.1,0.3,0.5], 
        '#jxn_bases_list': [0, 1, 2, 5, 10],
        'arm_number': [3, 4, 5, 6], 
        'flag_suffix': '', 
        'conf_suffix': '', 
        'sp_suffix': '', 
        'ns_dims': [20,2,7], 
        'figure_savepath': './',
        'varname': 'pa',
    }
    return data

def generate_data_for_batch_issue_k2_heatmap_plot():
    data = {
        'temp_list': [20, 23, 27, 30, 40, 50], 
        'conc_list': [0.05, 0.1,0.3,0.5], 
        '#jxn_bases_list': [0, 1, 2, 5, 10],
        'arm_number': [3, 4, 5, 6], 
        'flag_suffix': '', 
        'conf_suffix': '', 
        'sp_suffix': '', 
        'ns_dims': [20,2,7],
        'varname': 'k2' 
    }
    return data

def generate_data_for_batch_issue_pa_vstime_plot():
    data = {
        'temp_list': [20, 23, 27, 30, 40, 50], 
        'conc_list': [0.05, 0.1,0.3,0.5], 
        '#jxn_bases_list': [0, 1, 2, 5, 10],
        'arm_number': [3, 4, 5, 6], 
        'flag_suffix': '', 
        'conf_suffix': '', 
        'sp_suffix': '', 
        'ns_dims': [20,2,7],
        'varname': 'pa' 
    }
    return data

def generate_data_for_ns_pa_vstime_plot():
    data = {
        'arm_number' : 4,
        'temp' : 30,
        'conc' : 0.5,
        'flag_suffix' : '',
        'conf_suffix' : '',
        'sp_suffix' : '',
        'ns_dims' : [20,2,7],
    }
    return data

def run_ns_pa_plot_summary():
    data = generate_data_for_ns_pa_summary_plot()
    from plot_tasks.ns_pa_plot_summary import ns_pa_summary_plot
    ns_pa_summary_plot(data)

def run_ns_k2_heatmap_plot_batch_issue():
    data = generate_data_for_batch_issue_k2_heatmap_plot()
    from plot_tasks.plot_tasks_batch_issue import batch_issue_k2_heatmap
    batch_issue_k2_heatmap(data)

def run_ns_pa_vstime_plot_batch_issue():
    data = generate_data_for_batch_issue_pa_vstime_plot()
    from plot_tasks.plot_tasks_batch_issue import batch_issue_pa_vstime
    batch_issue_pa_vstime(data)

def run_ns_pa_vstime_plot_single_shot():
    data = generate_data_for_ns_pa_vstime_plot()
    from calc_tasks.stacking_local_identify_calc import stacking_local_identify_calc
    stacking_local_identify_calc(data)
    # from plot_tasks.ns_pa_plot_vstime import ns_pa_vstime_plot
    # ns_pa_vstime_plot(data)

if __name__ == '__main__':
    run_ns_pa_plot_summary()
    run_ns_k2_heatmap_plot_batch_issue()
    run_ns_pa_vstime_plot_batch_issue()
    run_ns_pa_vstime_plot_single_shot()