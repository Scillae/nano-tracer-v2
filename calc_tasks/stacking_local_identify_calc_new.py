from parse_tasks import ns_construct
import numpy as np
from infra import TimeSeries
from utils.result_cache import SL_result_cache

from collections import OrderedDict

def stacking_local_identify_calc(data:dict):
    '''
    
    '''
    varname = 'si'
    # load data from result_cache
    SL_result_cache(data, varname, 'load')
    if data['SL_content']:
        cached_result = data['SL_content']
        data['SL_content'] = None
        return cached_result

    # not cached
    # load ns series
    ns_series = ns_construct(data)
    # determine if the nanostar is RNA-based
    is_RNA = bool('RNA' in data)
    
    # construct result series of stack event identifying
    stacking_results = TimeSeries() # OrderedDict{t_stamp: is_stacking, stack_ls}
    for t_stamp, ns in ns_series.items():
        
        is_stacking = False # default to be False
        stack_ls = [] # [(adj_bp, adj_bp2, arm_id, arm2_id)] stores the 'unpaired' central bases adjacent to the 1st bps.
        shift_dic = {} # can be int or tuple. int --> which bp in the arm is the 1st; tuple: use the provided
        mismatch_b_ls = [] # pool of all mis-matched bases (unpaired bases that belong to an arm and are expected to be paired)
        arms_IDs = list(ns.arms.keys())
        





        # fix the nanostar's broken pairs and mis-matched bps. (Could be a part of the models: re-matched_ns_construct.py?)
        # detect the 1st paired bp in arm & create the mis-matched bases pool
        for ia0 in arms_IDs:
            arm0 = ns.arms[ia0]
            shift_dic[ia0] = 0
            for i in range(data['ns_dims'][0]):
                if check_if_bp(arm0.base_pairs[i+1][0],arm0.base_pairs[i+1][1]): # idx of base_pairs starts from 1
                    mismatch_b_ls.extend(arm0.base_pairs[i+1])
                    shift_dic[ia0] = i + 1 # stores where the arm0 'truly' start ---- excluding the mis-matched base pairs.
        







        # pair up the mis-matched bps if possible; it NO LONGER follow the topology of the nanostar! Only depending on the bases' standalone characteristics
        if mismatch_b_ls: # empty is False; enter here if not empty
            mismatch_bp_ls = [] # where the 'paired' mis-matched bps are stored
            used_p_ls = [] # where the examined mis-matched bases are stored, preventing redundant calculation
            for b in mismatch_b_ls:
                if b in used_p_ls: # skip if already examined
                    continue
                dist_ls = [np.linalg.norm(np.array(b.position)-np.array(p.position)) for p in mismatch_b_ls] # b is the base working on, p is the base we iterate
                p_ls = sorted(zip(dist_ls,mismatch_b_ls)) # sort to iterate from the closest base
                p_ls.pop(0) # discard the first one: dist = 0, itself
                for d, p in p_ls: # distance d is named for possible future use
                    if check_if_bp(b, p): # store 'paired' mis-matched bases in pair; register in used base pool
                        mismatch_bp_ls.append((b,p))
                        used_p_ls.append(p)

            
            
            
            
            # find which arm the mismatched should append to
            fbp_dic = {} # {(first base pair) : index of arm}
            for ia in arms_IDs:
                fbp_dic[ns.arms[ia].base_pairs[shift_dic[ia]]] = ia # all are the 'true' start (broken bps excluded)
                shift_dic[ia] = {} #  shift_dic[ia] = {dist: bp}
            for bp in mismatch_bp_ls:
                # pool all 1st bps, find the closest, check if normal aligned w/ disp
                fbp_dist_ls = [CoM_bp(fbp)-CoM_bp(bp) for fbp in fbp_dic.keys()] # find the distance between the 'true' start and the re-matched bp
                match_cand_ls = sorted(zip(fbp_dist_ls, fbp_dic.keys()),key=lambda x: np.linalg.norm(x[0])) # sort the candidate re-matched bps according to distance
                bp_normal = (np.array(bp[0].normal)-np.array(bp[1].normal))/2 # find the basepair's normal vector. Notice that the direction may be flipped!
                for dist, fbp in match_cand_ls:
                    # i = -1 # this i stores where the re-matched bp locates in the arm: -1 is right before the 1st bp; -2 is one even before -1, etc.
                    fbp_normal = (np.array(fbp[0].normal)-np.array(fbp[1].normal))/2 # Notice that the direction may be flipped! use absolute value of dot to find the angle
                    if (np.linalg.norm(dist) < 0.5 * 1.5) and np.abs(np.dot(bp_normal,fbp_normal)) > 0.7: # typical_dist_stk_bps(0.5) * 150% threshold
                        shift_dic[fbp_dic[fbp]][np.linalg.norm(dist)] = bp
        
        
        
        
        
        
        
        
        # modify the nanostar object according to the 'true' start
        for ia, fbp_idx in shift_dic.items():
            if type(fbp_idx) is int:
                if fbp_idx == 1: # intact arm; no modification required
                    continue
                arm = ns.arms[ia]
                arm.base_pairs = generate_new_arm_base_pairs(arm,fbp_idx, tmp_time=t_stamp) # generate a new arm that starts from the 'true' start
                ns.arms[ia] = arm
            else: # fbp_idx == {dist: bp}
                arm = ns.arms[ia]
                arm.base_pairs = generate_new_arm_base_pairs(arm,-1*len(fbp_idx), bp_append = [d_bp[1] for d_bp in sorted(list(fbp_idx.items()), key=lambda x:x[0], reverse=True)], tmp_time=t_stamp) # generate a new arm that starts from the re-matched bps; the farthest re-matched bp should be the 1st bp of the arm
                ns.arms[ia] = arm
       
        
        
        
        
        
        
        # identify the stacked arm pairs
        # generate potential stk_arm_pairs
        ptt_aps = [] # [(arm1,arm2)]
        arms_IDs = list(ns.arms.keys())
        for idx_0, ia0 in enumerate(arms_IDs):
            arm0_dir = get_arm_dir(ns.arms[ia0], is_RNA) # get the arm's direction
            for idx_1 in range(idx_0+1, len(arms_IDs)):
                arm1_dir = get_arm_dir(ns.arms[arms_IDs[idx_1]],is_RNA)
                # find if angle~ {arm_dir0:arm_dir1} < threshold
                if obtain_cos(arm0_dir, arm1_dir*(-1)) < 30: # one arm reversed. 30deg --> ang==150deg. </EMPIRICAL THRESHOLD/>
                    ptt_aps.append((arm0, ns.arms[arms_IDs[idx_1]]))
        
        
        
        
        
        
        
        # check if no potential arm pairs; skip if so
        if not ptt_aps:
            stacking_results[t_stamp]=(is_stacking, stack_ls)
            continue
        
        # identify & create junction base pairs (jxn_bps)
        # why only junction base pairs matter: if bps are broken, they have no effect on calculation; if bps are paired, either it's good arm pairs or re-matched arm pairs. Both are taken care in the previous nanostar-fixing stage
        jxn_bps = [(b0,b1) for b0 in ns.center.values() for b1 in ns.center.values() if check_if_bp(b0,b1) and b0.base_id < b1.base_id]
        
        
        
        
        
        
        
        # check if potential stackings are truly stackings
        for arm0, arm1 in ptt_aps:
            armbp0 = arm0.base_pairs[1]
            armbp1 = arm1.base_pairs[1]
            incyl_cnt = 0 # how many bps are in between the 1st bps of the selected two arms (armbps)

            
            
            
            
            for jxn_bp in jxn_bps:
                # check if the jxn_bp does not lie between the 1st bps
                ang0 = obtain_cos(v_normalize(CoM_bp(jxn_bp)-CoM_bp(armbp0)),v_normalize(CoM_bp(armbp1)-CoM_bp(armbp0))) # angle between {jxn_bp:armbp0, armbp1:armbp0}
                ang1 = obtain_cos(v_normalize(CoM_bp(jxn_bp)-CoM_bp(armbp1)),v_normalize(CoM_bp(armbp0)-CoM_bp(armbp1))) # angle between {jxn_bp:armbp1, armbp0:armbp1}
                if ang0 > 90 or ang1 > 90:
                    continue

                
                
                
                
                # check if dist_ax < threshold. dist_ax: dist of CoM_bp from line connecting the two armbps. (perpendicular to the line connecting the armbps)
                dist_ax = np.sin(ang0) * np.linalg.norm(CoM_bp(jxn_bp)-CoM_bp(armbp0))
                # </EMPIRICAL THRESHOLD/> 
                # 1.5: 150% scaling of the entire threshold
                # 1.3: typical value of the distance between the stacked junction base pairs and the central line (totally empirical)
                # optional multipliers: 
                # np.sin(15*np.pi/180): the (dumb) correction responding to the fact that the stacked arms form an ~150 deg angle (instead of 180 deg). A jxn_bp perfectly aligned with one armbp would be 15 deg off from the central line.
                # incyl_cnt: though never actually used, this term reminds me that the interactions between junction base pairs themselves should be considered as well.
                if dist_ax < 1.5 * 1.3: 
                    incyl_cnt += 1
            
            
            
            
            
            
            
            # check if distance between the armbps < threshold
            if np.linalg.norm(CoM_bp(armbp0) - CoM_bp(armbp1)) < 0.4 * (incyl_cnt + 1) * 1.5 * 2: # do we want this arc_approx?
                # </EMPIRICAL THRESHOLD/> 
                # typical distance between stacked base pairs: 0.5 (0.4 for RNA)
                # (incyl_cnt + 1): not empirical, the total number of segments
                # 1.5: 150% tolerance of the threshold settings
                # 2: 200% tolerance of the (thermal) vibration of the junction base pairs: they are not stacked as being in a true dsDNA; the force binding them stacked is much weaker.
                stack_ls.append((None, None, arm0.arm_id,arm1.arm_id)) # it was stack_ls.append((adj_bp, adj_bp2, arm.arm_id,arm2.arm_id)), now set to None for compatibility & future use
        
        # log the result
        is_stacking = True if stack_ls else False
        stacking_results[t_stamp]=(is_stacking, stack_ls)
    
    # save in result_cache
    data['SL_content'] = stacking_results
    SL_result_cache(data, varname, 'save')
    return stacking_results

def generate_new_arm_base_pairs(arm,starting_loc, tmp_time=0, bp_append = None): 
    '''
    Generate a new arm with its length adjusted (truncated because of broken bps, or amended because of re-matched bps).
    :arm: the original arm; 
    :starting_loc: the starting location of the arm (positive results in a shorter arm, negative make the arm elongated); 
    :bp_append: bps to be appended (append in order: 1st in bp_append becomes the 1st in the new_arm)
    '''
    new_arm = OrderedDict()
    if bp_append is not None:
        assert starting_loc < 0
        assert starting_loc*(-1) == len(bp_append)
        for i in range(starting_loc*(-1)):
            new_arm[i+1] = bp_append[i] # arm's base pair id starts from 1
        for ibp, bp in arm.base_pairs.items():
            new_arm[ibp-starting_loc] = bp
    else:
        assert starting_loc > 0
        for ibp, bp in arm.base_pairs.items():
            if ibp < starting_loc:
                continue
            new_arm[ibp-starting_loc+1] = bp
    return new_arm

def check_if_bp(b1,b2):
    if b1.base_id == b2.base_id:
        return False
    is_oppos_dir = np.dot(np.array(b1.backbone),np.array(b2.backbone)) < 0
    is_normal_align = np.dot(np.array(b1.normal),np.array(b2.normal)) < 0 # anti-aligned
    disp = (np.array(b1.position)-np.array(b2.position))
    is_same_plane_b1 = np.abs(np.dot(v_normalize(disp),np.array(b1.normal))) < 0.3
    is_same_plane_b2 = np.abs(np.dot(v_normalize(disp),np.array(b2.normal))) < 0.3
    if np.linalg.norm(disp) == 0:
        return False
    is_disp_small = np.linalg.norm(disp) < 1.3 * 1.5 # normal value: ~1.3 SU. 2 is ~150%*1.3
    disp /= np.linalg.norm(disp)
    is_align_disp_b1 = np.abs(np.dot(disp, np.array(b1.backbone))) > 0.5
    is_align_disp_b2 = np.abs(np.dot(disp, np.array(b2.backbone))) > 0.5
    return is_oppos_dir and is_normal_align and is_disp_small and (is_align_disp_b1 or is_align_disp_b2) and (is_same_plane_b1 or is_same_plane_b2)

def arc_approx_corr(armbp0,armbp1,norm): # this function will approximate the distance between two arms' 1st bps as an arc instead of a straight line, according to the angle formed by the two arms
    ang = obtain_cos(np.array(armbp0[0].normal),np.array(armbp1[0].normal))
    central_ang = (180-ang)*2
    r_ang = (ang/2-(180-ang)/2)
    r = norm/np.sin(central_ang*(np.pi/180)) * np.sin(r_ang*(np.pi/180)) # law of sines
    return r

def get_arm_dir(arm,is_RNA):
    dir_bp1_bp2 = v_normalize(CoM_bp(arm.base_pairs[2])-CoM_bp(arm.base_pairs[1]))
    is_b1_rev = np.dot(np.array(arm.base_pairs[1][1].normal),dir_bp1_bp2) > 0 # check if b1 of bp1 points in the direction of arm
    # define arm_dir as the avg of normal of bp1
    arm_dir = v_normalize((np.array(arm.base_pairs[1][0].normal) * (is_b1_rev*2-1) * (-1) + np.array(arm.base_pairs[1][1].normal) * (is_b1_rev*2-1))/2) # if b1 not rev, b0 is rev.
    if is_RNA: # rotate to get true normal vector
        # get correct base CoM direction
        fbp_CoM_v = v_normalize(np.array(arm.base_pairs[1][0].position)-np.array(arm.base_pairs[1][1].position))
        fbp_CoM_v = fbp_CoM_v if np.dot(fbp_CoM_v, dir_bp1_bp2) > 0 else -fbp_CoM_v
        # ensuring orthogonality
        fbp_CoM_v = v_normalize(np.cross(np.cross(arm_dir,fbp_CoM_v),arm_dir))
        # rotate by 19deg
        arm_dir = np.cos(np.radians(19))*arm_dir + np.sin(np.radians(19))*fbp_CoM_v
    return arm_dir

def CoM_bp(bp):
    return (np.array(bp[0].position)+np.array(bp[1].position))/2

def v_normalize(v):
    return v/np.linalg.norm(v)

def obtain_cos(v1,v2):
    return np.degrees(np.arccos(np.sum(v1*v2)/(np.sqrt(np.sum(np.square(v1))) * np.sqrt(np.sum(np.square(v2)))))) # cos(x) = n1 * n2 / (|n1|*|n2|), angle only 0~pi, negative -> abs
