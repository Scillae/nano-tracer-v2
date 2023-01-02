from utils.tools import auto_generate_summary_range

def batch_issue_k2_heatmap(data:dict): # entry point for generating the overall patch angle summary plot of selected nanostars
    '''
    Issue a batch of k2 heatmap plot tasks
    :data: Expected to contain:
        data['temp_list']; data['conc_list']; data['arm_number']; data['flag_suffix']; 
        data['conf_suffix']; data['sp_suffix']; data['ns_dims']
    '''
    # for each nanostar to be summarized, produce a separate data dictionary
    data_package_list = auto_generate_summary_range(data)
    # iterate through all the individual nanostar
    # from joblib import Parallel, delayed # in Parallel
    # from plot_tasks.ns_k2_plot_heatmap import ns_k2_heatmap_plot
    # _ = Parallel(n_jobs=12)(delayed(ns_k2_heatmap_plot)(child_data) for child_data in data_package_list)
    from plot_tasks.ns_k2_plot_heatmap import ns_k2_heatmap_plot
    for child_data in data_package_list:
        ns_k2_heatmap_plot(child_data)
    return


def batch_issue_pa_vstime(data): # entry point for generating the overall patch angle summary plot of selected nanostars
    '''
    Issue a batch of patch angle vs time plot tasks
    :data: Expected to contain:
        data['temp_list']; data['conc_list']; data['arm_number']; data['flag_suffix']; 
        data['conf_suffix']; data['sp_suffix']; data['ns_dims']
    '''
    # for each nanostar to be summarized, produce a separate data dictionary
    data_package_list = auto_generate_summary_range(data)
    # iterate through all the individual nanostar
    # from joblib import Parallel, delayed # in Parallel
    # from plot_tasks.ns_pa_plot_vstime import ns_pa_vstime_plot
    # _ = Parallel(n_jobs=12)(delayed(ns_pa_vstime_plot)(child_data) for child_data in data_package_list)
    from plot_tasks.ns_pa_plot_vstime import ns_pa_vstime_plot
    for child_data in data_package_list:
        ns_pa_vstime_plot(child_data)
    return
