from models import Strand, Base, Arm
from collections import OrderedDict
import numpy as np
import copy


# Noella: Accepted_sci
# dist function should be defined before using
# Noella_v2:
# better moving to `utils`
def dist(t1, t2):
    return np.sqrt(np.sum(np.square(np.array(t1) - np.array(t2))))


class NanoStar:
    def __init__(self, strands_dic, dims_ls, arm_num, box_dim=None):  # box_dim hacking
        """
        init NanoStar
        :param strands_dic:
        :param dims_ls:
        :param arm_num:
        :param box_dim:
        """
        # revert all strands: done in TimeMachine.add_strands
        # strands_dic = OrderedDict()
        # sanity check of strand amount
        assert len(strands_dic) == arm_num
        strands = copy.deepcopy(strands_dic)  # avoid modification to the ref var.
        self.strands = strands_dic
        self.arms = self.binding(strands, dims_ls, box_dim)
        self.center = self.center_gen(self.strands,
                                      dims_ls)  # PBCC performed in self.binding, so self.center must goes afterwards. Also, 'strands' is modified by self.binding, so use self.strands here.
        self.box_dim = box_dim
        self.dims = dims_ls

    def center_gen(self, strands_dic, dims_ls) -> dict:
        """
        calculate center of strands
        :param strands_dic:
        :param dims_ls:
        :return:
        """
        len_arm, len_cen, len_end = dims_ls
        cen_dic = OrderedDict()
        st = 0
        # Noella: Accepted_sci
        # for strand in strands_dic:
        for _, strand in strands_dic.items():
            for i in range(len_cen):
                # Noella: Questioned_sci(base sequence is a dic, and thus should not be accessed by index.)
                # `strand.base_sequence` is a list defined in line 15, `base_seq_r_ls` (a list) is passed to contructor as required parameter `base_seq`. All `strand.base_sequence.values()` should be replaced
                # around lines: 31 44 48 68 71 (90) [may not be exactly because of new added comments lines]
                # cen_dic[st*len_cen+i] = strand.base_sequence.values()[len_arm+i]
                # dic or list: dict
                sel_b = list(strand.base_sequence.values())[len_arm + i]
                cen_dic[sel_b.base_id] = sel_b
            st += 1
        return cen_dic

    def binding(self, strands_dic, dims_ls, box_dim):
        """
        bind strands
        :param strands_dic:
        :param dims_ls:
        :param box_dim:
        :return:
        """
        strand_id_idx = list(strands_dic.keys())
        # find 3' head of an arbitrary strand
        s0 = strands_dic[min(strand_id_idx)]
        # pool all bases with their locations # now only binding points
        strand_id_idx, strands_dic, bind_base = self.pairing_single_step(strand_id_idx, strands_dic, s0, box_dim,
                                                                         dims_ls, True)
        # construct arm
        s1 = strands_dic[bind_base.strand_id]
        arm_idx = 0
        arm_dic = OrderedDict()
        arm_dic[arm_idx] = Arm(arm_idx, dims_ls, s0.strand_id, {s0.strand_id: s0, s1.strand_id: s1})
        # loop out all other strands
        while len(strands_dic) != 1:
            s0 = s1  # deepcopy?
            strand_id_idx.remove(s0.strand_id)
            strand_id_idx, strands_dic, bind_base = self.pairing_single_step(strand_id_idx, strands_dic, s0, box_dim,
                                                                             dims_ls, False)
            s1 = strands_dic[bind_base.strand_id]
            arm_idx += 1
            arm_dic[arm_idx] = Arm(arm_idx, dims_ls, s0.strand_id, {s0.strand_id: s0, s1.strand_id: s1})
            del strands_dic[s0.strand_id]
        return arm_dic

    def periodic_box_cond_correction(self, bind_base, strands_dic, box_dim):
        """
        pbc
        :param bind_base:
        :param strands_dic:
        :param box_dim:
        :return:
        """
        # Assumtions:
        # 'intended diffusion across box' only acting on continuous part, likely a whole strand.
        # 'intended diffusion across box' won't cast bases across box boundary.
        # safety check: is_diffusion. Disabled: may have 2 adjacent diffused strands.
        # dis_arr = np.array([dist(bind_base.position, base.position) if base is not bind_base else None for base in pool_base_ls])
        # dis_arr = dis_arr[dis_arr != np.array(None)]
        # assert all(dis_arr > 4) # confirming diffused strand: all greater than 5 (empirical value, binding:1~3)
        # PBCC
        # PBC-CoM Centering # no need for base_cnt since scaling sine and cosine simultaneously does not change the angle.
        phasors = np.zeros((2, 3))  # NOT normalized to 1. Instead, noted in sines/cosines
        for strand in strands_dic.values():
            for base in strand.base_sequence.values():
                angles = np.array(base.position) / (box_dim / (2 * np.pi))
                phasors += np.vstack((np.sin(angles), np.cos(angles)))
        CoM_pos = np.arctan2(-phasors[0], -phasors[1]) + np.pi  # Now in 0~2pi. '-' counterbalancing '+np.pi'.
        CoM_pos *= (box_dim / (2 * np.pi))
        # Update Position.
        for strand in strands_dic.values():
            for base in strand.base_sequence.values():
                old_pos = np.array(base.position)
                new_pos = old_pos - CoM_pos + box_dim / 2  # old_pos - CoM_pos = r_CoM. shift to mid-box afterward.
                _ = base.set_position(tuple(new_pos))
                _ = self.strands[base.strand_id].base_sequence[base.base_id].set_position(tuple(new_pos))
        # PBCC
        base = bind_base
        old_pos = np.array(base.position)
        # traverse backward, skipping bind_base
        while strands_dic[base.strand_id].base_sequence.get(base.prev_id) is not None and dist(old_pos, strands_dic[
            base.strand_id].base_sequence[base.prev_id].position) < 4:
            base = strands_dic[base.strand_id].base_sequence[base.prev_id]
            old_pos = np.array(base.position)
            new_pos = ((old_pos % box_dim) + box_dim) % box_dim  # NOT equavalent to (old_pos + box_dim) % box_dim
            _ = base.set_position(tuple(new_pos))
            # fix self.strands as well
            _ = self.strands[base.strand_id].base_sequence[base.base_id].set_position(tuple(new_pos))
        # modifying bind_base
        base = bind_base
        old_pos = np.array(base.position)
        new_pos = ((old_pos % box_dim) + box_dim) % box_dim
        _ = base.set_position(tuple(new_pos))
        _ = self.strands[base.strand_id].base_sequence[base.base_id].set_position(tuple(new_pos))
        # traverse forward, skipping bind_base
        while strands_dic[base.strand_id].base_sequence.get(base.next_id) is not None and dist(old_pos, strands_dic[
            base.strand_id].base_sequence[base.next_id].position) < 4:
            base = strands_dic[base.strand_id].base_sequence[base.next_id]
            old_pos = np.array(base.position)
            new_pos = ((old_pos % box_dim) + box_dim) % box_dim
            _ = base.set_position(tuple(new_pos))
            # fix self.strands as well
            _ = self.strands[base.strand_id].base_sequence[base.base_id].set_position(tuple(new_pos))
        return strands_dic  # still, bases not in the same strand may be close to each other!
        # return self.pbc_CoM_centering(strands_dic,box_dim) 

    def pairing_single_step(self, strand_id_idx, strands_dic, s0, box_dim, dims_ls, is_checking_pbcc):
        """

        :param strand_id_idx:
        :param strands_dic:
        :param s0:
        :param box_dim:
        :param dims_ls:
        :param is_checking_pbcc:
        :return:
        """
        binding_location = 4  # controlling how far the binding location is from the end of arm. (4:= quadrisection point. 2:= bisection point)
        if is_checking_pbcc:
            CoM_pos = self.get_CoM_pos()
            criteria = (dims_ls[0] + dims_ls[1]) / (
                    20 - binding_location) * 8.3 * 1.2  # 20 bp : ~8.3 SU. Setting tolerance as 20%
            bp_ls = [
                list(strands_dic[idx].base_sequence.values())[-1 - dims_ls[2] - dims_ls[0] // binding_location] if dist(
                    np.array(list(strands_dic[idx].base_sequence.values())[
                                 -1 - dims_ls[2] - dims_ls[0] // binding_location].position),
                    CoM_pos) > criteria else None for idx in strand_id_idx]  # distances from binding point to the CoM
            bp_ls.extend([list(strands_dic[idx].base_sequence.values())[
                              dims_ls[2] + dims_ls[0] // binding_location] if dist(np.array(
                list(strands_dic[idx].base_sequence.values())[dims_ls[2] + dims_ls[0] // binding_location].position),
                CoM_pos) > criteria else None for idx
                          in strand_id_idx])
            if any(bp is not None for bp in bp_ls):
                for base in bp_ls:
                    if base is not None:
                        # send strands with largest distance to pbcc
                        strands_dic = self.periodic_box_cond_correction(base, strands_dic, box_dim)
                # recursion
                strand_id_idx, strands_dic, bind_base = self.pairing_single_step(strand_id_idx, strands_dic, s0,
                                                                                 box_dim, dims_ls,
                                                                                 is_checking_pbcc)  # careful -- recursion
                return strand_id_idx, strands_dic, bind_base
        test_length = 3  # adjust if arm is not long enough. direction: center?dedicated to current NS designs.
        test_distance = 4  # paired bases cannot be more than 'test_distance' apart.
        s0_head = list(s0.base_sequence.values())[dims_ls[0] // binding_location]
        pool_base_ls = []
        for idx in strand_id_idx:
            if idx == s0.strand_id:
                continue
            bind_point_1 = list(strands_dic[idx].base_sequence.values())[
                -1 - dims_ls[2] - dims_ls[0] // binding_location]
            bind_point_2 = list(strands_dic[idx].base_sequence.values())[dims_ls[2] + dims_ls[0] // binding_location]
            # if is_checking_pbcc:
            #     # Centering bind_points w/o pbcc
            #     bp1_CoM_pos = np.array(bind_point_1.position) - CoM_pos
            #     bp2_CoM_pos = np.array(bind_point_2.position) - CoM_pos
            #     # Check if need PBCC. 20 bp : ~8.3 SU. Setting tolerance as 50%
            #     is_b1_pbcc = np.linalg.norm(bp1_CoM_pos) > (dims_ls[0]+dims_ls[1])/(20-binding_location)*8.3*1.5
            #     is_b2_pbcc = np.linalg.norm(bp2_CoM_pos) > (dims_ls[0]+dims_ls[1])/(20-binding_location)*8.3*1.5
            #     if is_b1_pbcc or is_b2_pbcc:
            #         # check if binding points 
            #         if is_b1_pbcc:
            #             strands_dic = self.periodic_box_cond_correction(bind_point_1, pool_base_ls, strands_dic, box_dim)
            #         else:
            #             strands_dic = self.periodic_box_cond_correction(bind_point_2, pool_base_ls, strands_dic, box_dim)
            #         strand_id_idx, strands_dic, bind_base = self.pairing_single_step(strand_id_idx, strands_dic, s0, box_dim, dims_ls, is_checking_pbcc) # careful -- recursion
            #         return strand_id_idx, strands_dic, bind_base
            pool_base_ls.extend([bind_point_1,
                                 bind_point_2])  # selecting the pairing node, no longer the whole strand. But Not Assuming Direction.
        # get the closest base
        min_dist = 1000000
        bind_base = None
        for candidate_base in pool_base_ls:
            dis = dist(candidate_base.position, s0_head.position)
            orient_dot = np.dot(np.array(s0_head.backbone), np.array(candidate_base.backbone))
            direction = np.array(candidate_base.position) - np.array(s0_head.position)
            direction_dot = np.dot(direction / np.linalg.norm(direction),
                                   s0_head.backbone)  # check if candidate locates near the direction that s0_head's backbone vector points at.
            if dis < min_dist and orient_dot < 0 and dis < 4 and direction_dot > 0:  # and orient_dot < -0.7 and dis < 4 and direction_dot > 0.7
                min_dist = dis
                bind_base = candidate_base
        assert bind_base is not None
        # if bind_base is None:
        #     # select the very far one w/ backbone more negative as bind_base.
        #     cand_base_ls = []
        #     for bind_base in pool_base_ls:
        #         dis_arr = np.array([dist(bind_base.position, base.position) if base is not bind_base else None for base in pool_base_ls])
        #         dis_arr = dis_arr[dis_arr != np.array(None)]
        #         if sum(dis_arr > 0.9*box_dim[0]) >= 6: # being far from at least 3 strands simultaneously
        #             cand_base_ls.append(bind_base)
        #     orient_dot_arr = np.array([np.dot(np.array(s0_head.backbone),np.array(cand_base.backbone)) for cand_base in cand_base_ls])
        #     bind_base = cand_base_ls[np.argmin(orient_dot_arr)] # the little angle between backbones is neglected. use np.abs to do the reflection.      
        # check if periodic box condition correction should be applied.
        test_dist_ls = [dist(s0.base_sequence[s0_head.base_id - i].position,
                             strands_dic[bind_base.strand_id].base_sequence[bind_base.base_id + i].position) for i in
                        range(test_length)]
        assert all(np.array(test_dist_ls) < test_distance)
        return strand_id_idx, strands_dic, bind_base

    def get_CoM_pos(self):
        """

        :return:
        """
        strands_dic = self.strands
        # CoM
        CoM_pos = np.zeros(3)
        base_cnt = 0
        for strand in strands_dic.values():
            for base in strand.base_sequence.values():
                CoM_pos = np.add(CoM_pos, np.array(base.position))
                base_cnt += 1
        CoM_pos = np.divide(CoM_pos, base_cnt)
        return CoM_pos

        # def pbc_CoM_centering(self,strands_dic,box_dim):
        #     # Lazy box-centering: only enabled when PBCCing. Consider making it global?
        #     # Assuming PBC-centering unrelated to PBCC. Perform after PBCC.
        #     # CoM
        #     CoM_pos = np.zeros(3)
        #     base_cnt = 0
        #     for strand in strands_dic.values():
        #         for base in strand.base_sequence.values():
        #             CoM_pos = np.add(CoM_pos, np.array(base.position))
        #             base_cnt += 1
        #     CoM_pos = np.divide(CoM_pos, base_cnt)
        #     # PBC & CoM centering
        #     for strand in strands_dic.values():
        #         for base in strand.base_sequence.values():
        #             old_pos = np.array(base.position)
        #             new_pos = old_pos - CoM_pos + box_dim/2 # old_pos - CoM_pos = r_CoM. shift to mid-box afterward.
        #             ret_pos = base.set_position(tuple(new_pos))
        #             assert all(ret_pos == new_pos)
        #             ret_pos = self.strands[base.strand_id].base_sequence[base.base_id].set_position(tuple(new_pos))
        #             assert all(ret_pos == new_pos)
        #     return strands_dic

        '''
        len_arm, len_cen, len_end = dims_ls
        len_tot = 2*len_arm + len_cen + len_end
        binding_pts_2dls = []
        for strand in strands_dic.values():
            str_ls = list(strand.base_sequence)
            binding_pts_2dls.append(list(str_ls[0],str_ls[len_end],str_ls[len_tot-1-len_end],str_ls[len_tot-1]))
        # now find the min-dist
        #TODO: THE LOOP HERE MUST BE REDESIGNED!!!
        bind_tp_ls = []
        for i in range(len(binding_pts_2dls)):
            idx_arr = np.arange(len(binding_pts_2dls))
            idx_arr = idx_arr[idx_arr != i] # double-counting if '!=' ? Must Double Count! No Need to double count!
            min_dist_dic = OrderedDict()
            for j in idx_arr: # still need to compare 'ALL' pairs of DNAs
                s0 = binding_pts_2dls[i]
                s1 = binding_pts_2dls[j]
                min_dist_pair_dic = OrderedDict()
                for k in range(len(s0)):
                    for l in range(k+1, len(s1)):
                        min_dist_pair_dic[dist(s0[k], s1[l])] = [s0[k], s1[l]]
                        # min_dist_pair_ls.append((dist(s0[k], s1[l]), s0[k], s1[l]))
                # min_dist_ls.append(min(min_dist_pair_ls))
                m_d = min(list(min_dist_pair_dic.keys()))
                m_v = min_dist_pair_dic[m_d]
                min_dist_dic[j] = m_v.insert(0, m_d) # Find min in 16 pts. This generate tot_i * tot_j min. In each i, 2 comparable min out of tot_j should be found.
                # min_dist_dic: 5 items, each j: [s0_b, s1_b, dist]
            min_d_0 = min([m_v_ls[2] for m_v_ls in list(min_dist_dic.values())])
            min_item_0 = [item for item in min_dist_dic.items() if item[1][2] == min_d_0][0] # select tuples (j,[s0,s1,d]) if d==d_0
            min_dist_dic.pop(min_item_0[0])
            min_d_1 = min([m_v_ls[2] for m_v_ls in list(min_dist_dic.values())])
            min_item_1 = [item for item in min_dist_dic.items() if item[1][2] == min_d_0][0] # 2 contact points for each strand with all other ones
        '''
        '''
            min_dist_cp = min_dist_ls[:] # avoid potential reference/shallow copy issue
            min_2dists_ls = [tp for tp in min_dist_cp if not any(tp[0] > t[0] for t in min_dist_ls)] # unsure about this line
            min_dist_cp = [tp for tp in min_dist_ls if tp not in min_2dists_ls] # drop lowest
            min_2dists_ls.extend([tp for tp in min_dist_cp if not any(tp[0] < t[0] for t in min_dist_cp)]) # extend only 1 element.
        '''
        '''
                s0 = binding_pts_2dls[i]
                s1 = binding_pts_2dls[j]
                dist = [dist(s0[0].position, s1[2].position),
                        dist(s0[0].position, s1[1].position),
                        dist(s0[1].position, s1[0].position),
                        dist(s0[1].position, s1[3].position),
                        dist(s0[2].position, s1[3].position),
                        dist(s0[2].position, s1[0].position),
                        dist(s0[3].position, s1[1].position),
                        dist(s0[3].position, s1[2].position),
                ]
                min_dist_val = min(dist)
                min_dist_idx = dist.index(min_dist_val)
                # check if actually binding #
                dist.sort()
                if dist[-2] > dist[-1]*2:
                    continue # break out to the next
                # end check #
                # 2 binding, possiblilities: 0~0, 0~41, 41~0, 41~41
                # NOTE if strands are reversed here, they would be reversed again since each strand counts twice? No! No double count! But I cannot reverse it here...
                if s0[min_dist_idx].base_id == s0[0].base_id: #s0 leading normal
                    if s1[min_dist_idx].base_id == s1[1].base_id: #s1 normal
                        bind_tp_ls.append(i, (s0[min_dist_idx], s1[min_dist_idx]), (strands_dic[s0[min_dist_idx].strand_id], strands_dic[s1[min_dist_idx].strand_id]), (False, False))
                    else: # s1[3], s1 reversed
                        bind_tp_ls.append(i, (s0[min_dist_idx], s1[min_dist_idx]), (strands_dic[s0[min_dist_idx].strand_id], OrderedDict(reversed(list(strands_dic[s1[min_dist_idx].strand_id].items())))), (False, True))
                    
                elif s0[min_dist_idx].base_id == s0[3].base_id: #s0 leading reversed
                    if s1[min_dist_idx].base_id == s1[1].base_id: #s1 normal
                        bind_tp_ls.append(i, (s0[min_dist_idx], s1[min_dist_idx]), (OrderedDict(reversed(list(strands_dic[s0[min_dist_idx].strand_id].items()))), strands_dic[s1[min_dist_idx].strand_id]), (True, False))
                    else: # s1[3]
                        bind_tp_ls.append(i, (s0[min_dist_idx], s1[min_dist_idx]), (OrderedDict(reversed(list(strands_dic[s0[min_dist_idx].strand_id].items()))), OrderedDict(reversed(list(strands_dic[s1[min_dist_idx].strand_id].items())))), (True, True))

                elif s1[min_dist_idx].base_id == s0[2].base_id: #s0 ending normal

                    











        tp_base_ls = []
        bt_base_ls = []
        for strand_tp in strands_dic.items():
            (strand_id, strand) = strand_tp
            bt_base_ls.append(strand.base_sequence[max(list(strand.base_sequence.keys()))])
            tp_base_ls.append(strand.base_sequence[min(list(strand.base_sequence.keys()))])
        bind_ls = []
        for tp_base in tp_base_ls:
            pair_tp = (0,0,0)
            min_dist_dbl = 100000000000000000.00 #inf
            for bt_base in bt_base_ls:
                dist_dbl = np.sqrt(np.sum(np.square(np.array(tp_base.position) - np.array(bt_base.position))))
                assert min_dist_dbl != dist_dbl
                if dist_dbl < min_dist_dbl:
                    min_dist_dbl = dist_dbl
                    pair_tp = (tp_base,bt_base,dist_dbl)
            bind_ls.append(pair_tp)
        print(bind_ls)
        '''


'''
import numpy as np
from models import Arm
class NanoStar:
    def __init__(self, strands, dims_ls, bind_ls):
        self.arms = self.construct_arms(strands, dims_ls, bind_ls)
        self.single_ends = self.get_single_ends()

    # construct arms, main work is to determine strand pair,
    # then construct by Arm's ctor
    # TODO how to present centers: move Arm.pair() here, and construct Center class.
    def construct_arms(self, strands, dims_ls, bind_ls):
        arms_ls = []
        for i in range(len(bind_ls)):
            tp_base, bt_base, _ = bind_ls[i]
            a = Arm(i, dims_ls,[strands[tp_base.strand_id], strands[bt_base.strand_id]])
            arms_ls.append(a)
        return arms_ls


    # get single ends from self.arms
    def get_single_ends(self):
        pass
'''
