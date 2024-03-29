from utils import get_ns_params, SL_result_cache, generate_path, chkdir, data_query
# import plot_tasks.ns_si_process.data_process_func
import plot_tasks.ns_si_process
import scipy.stats as stat
import matplotlib.pyplot as plt

import numpy as np
import os

def ns_k2_heatmap_plot(data:dict):
    '''
    Plot single nanostar's oblacity vs prolacity heatmap.
    :data: Expected to contain: (WIP)
        data['temp']; data['conc']; data['arm_number']; data['flag_suffix']; 
        data['conf_suffix']; data['sp_suffix']; data['ns_dims']; data['figure_savepath']
    '''
    # retrieve result from calc_func.stacking_local_identify_calc
    varname = 'si'
    # load data from result_cache
    SL_result_cache(data, varname, 'load')
    if not data['SL_content']:
        from calc_tasks.stacking_local_identify_calc import stacking_local_identify_calc as slic
        slic(data)
        SL_result_cache(data, varname, 'load')
        # raise Exception('The results trying to be plotted is not cached. Run calc_tasks.stacking_local_identify_calc first.')
    stack_info = plot_tasks.ns_si_process.data_process_func(data['SL_content'], data)
    data['SL_content'] = None
    
    # retrieve result from calc_func.k2_calc
    varname = 'k2'
    # load data from result_cache
    SL_result_cache(data, varname, 'load')
    if not data['SL_content']:
        from calc_tasks.k2_calc import k2_calc as k2c
        k2c(data)
        SL_result_cache(data, varname, 'load')
        # raise Exception('The results trying to be plotted is not cached. Run calc_tasks.k2_calc first.')
    var_vals = data_process_func(data['SL_content'], data)
    data['SL_content'] = None
    
    #### plot confs ####
    x_lim = (1,3)
    y_lim = (1,2)
    plot_path = generate_path(data, 'plot_path')
    label = generate_path(data, 'label')
    #### conf ends ####

    plot_confs = {'x_lim':x_lim, 'y_lim':y_lim, 'plot_path':plot_path, 'label':label}
    heatmap_single_ns_all_steps(var_vals, stack_info, plot_confs, data) # plot saved inside the function
    
    return


def heatmap_single_ns_all_steps(var_vals:dict, stack_info:dict, plot_confs:dict, data:dict):
    '''
    To be documented
    '''
    x_lim, y_lim, plot_path, label = data_query(plot_confs, ['x_lim','y_lim','plot_path','label'],['tuple','tuple','str','str'])

    _,_,_,_,_,_,_,ns_struc = get_ns_params(data['arm_number'])
    var_dic = data_process_func(var_vals, data)

    # load data
    lbd_axs_ls = var_dic['axes']
    l2_arr = np.array([conf_lbd[1][0]/conf_lbd[0][0] for conf_lbd in lbd_axs_ls])
    l3_arr = np.array([conf_lbd[2][0]/conf_lbd[1][0] for conf_lbd in lbd_axs_ls])
    
    # draw heatmap
    draw_k2_heatmap(l2_arr, l3_arr, x_lim, y_lim, stack_info, plot_path)
    return

def draw_k2_heatmap(l2_arr, l3_arr, x_lim, y_lim, stack_info, plot_path):
    '''
    (Truncated Version) Draw heatmap.
    '''
    title_prompt = 'FWHM'
    title_not_reverse = '' # legacy dumb var
    xmin, xmax = x_lim
    ymin, ymax = y_lim

    # initialize meshgrid and values
    X, Y = np.mgrid[xmin:xmax:3000j, ymin:ymax:3000j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([l2_arr, l3_arr])
    
    # try to generate gaussian kernel
    try:
        kernel = stat.gaussian_kde(values,bw_method=0.2)
    except:
        raise Exception('Singular matrix! Maybe too few points.')
    
    # calculate Z
    Z = np.reshape(kernel(positions).T, X.shape) # local density of kernel, not normalized, total should be the number of points
    Z = np.rot90(Z) # maybe I wrote the shape wrong somehow? It works now.
    Z_FWHM = generate_Z_layer(Z, 0, is_FWHM=True).astype(bool) # a ring with specific criteria. Number indicating percentile. FWHM=True: Full Width at Half Maximum
    Z_rings = Z_FWHM # add multiple selected rings together
    Z_draw = Z*((~Z_rings).astype(int)) # erase Zs the rings

    # plot
    fig, ax = plt.subplots()
    ax.imshow(Z_draw, cmap=plt.cm.gist_earth_r, extent=[xmin, xmax, ymin, ymax])
    # ax.scatter(l2_arr, l3_arr, s=0.4,c='#FFFFFF')
    
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])
    ax.xaxis.set_major_locator(plt.MultipleLocator(0.5))
    ax.yaxis.set_major_locator(plt.MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.1))
    ax.tick_params(bottom=True,top=True,left=True,right=True,direction='in')
    ax.tick_params(bottom=True,top=True,left=True,right=True,direction='in',which='minor')
    
    chkdir(os.path.dirname(plot_path))
    plt.savefig(os.path.splitext(plot_path)[0]+'-Heatmap'+'.png',dpi=400)
    plt.close()
    return

def generate_Z_layer(Z,prob,delta=0.05,is_FWHM=False):
    '''
    Generate a layer with specified criteria. (Not generating a uniform ring, which is hard.)
    :prob: percentile/100
    :delta: (+-) width of the circle in prob. Actual width = delta*2
    :is_FWHM: if overwrite the condition by FWHM
    '''
    if is_FWHM:
        Z_threshold = (np.max(Z))/2 # delta at center now
    else:
        prob = 1-prob
        Z_sum = np.sum(Z)
        Z_thresholds_arr = np.linspace(np.min(Z),np.max(Z),1000)
        Z_threshold = 0
        for Z_th in Z_thresholds_arr:
            if np.sum(Z*(Z>Z_th)) / Z_sum <= prob: # integral of values above Z_th drops to 1-prob
                Z_threshold = Z_th
                break
    # Produce Z
    Z_sub = Z*(Z>Z_threshold-delta/2)*(Z<Z_threshold+delta/2)
    return Z_sub

def data_process_func(k2_ls_res, data):
    t_ls = [i for i in k2_ls_res]
    k2_ls = [i[0] for i in k2_ls_res.values()]
    lbd_ls = [sorted(i[1]) for i in k2_ls_res.values()]
    del_ls = [(i[1]-i[0],i[2]-i[1]) for i in lbd_ls]
    axs_ls = [sorted(zip(i[1],i[2])) for i in k2_ls_res.values()] # [((lbds),(axes))] instead of [(axes)]
    return {'t':t_ls, 'k2':k2_ls, 'lambdas':lbd_ls, 'deltas':del_ls, 'axes':axs_ls}