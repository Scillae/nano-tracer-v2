from typing import Optional

from models import Strand, Base, Arm
from utils.tools import dist
from collections import OrderedDict
import numpy as np
import copy


class NanoStar:
    def __init__(self, strands: dict, ns_dims: list, arm_num: int, box_dim=None):  # box_dim hacking
        """
        init NanoStar
        :param strands: all strands used to construct the nanostar instance
        :param ns_dims: a list describing the dimensions of the nanostars, [len_arm, len_cen, len_end]
        :param arm_num: how many arms the nanostar has
        :param box_dim: size of the simulation box
        """
        assert len(strands) == arm_num
        strands_editing = copy.deepcopy(strands)  # avoid modification to the ref var.
        self.strands = strands
        self.arms = self.binding(strands_editing, ns_dims, box_dim)
        self.center = self.center_gen(self.strands,
                                      ns_dims)  # PBCC performed in self.binding, so self.center must goes afterwards. Also, 'strands' is modified by self.binding, so use self.strands here.
        self.box_dim = box_dim
        self.dims = ns_dims

    def center_gen(self, strands: dict, ns_dims: list) -> dict:
        """
        identify the center of strands
        :param strands: all strands
        :param ns_dims: a list describing the dimensions of the nanostars, [len_arm, len_cen, len_end]
        :return: central bases
        """
        len_arm, len_cen, len_end = ns_dims
        central_bases = OrderedDict()
        st = 0
        # Noella: Accepted_sci
        # for strand in strands:
        for _, strand in strands.items():
            for i in range(len_cen):
                sel_b = list(strand.base_sequence.values())[len_arm + i]
                central_bases[sel_b.base_id] = sel_b
            st += 1
        return central_bases

    def binding(self, strands: dict, ns_dims: list, box_dim) -> dict:
        """
        bind strands to form arms
        :param strands: all strands
        :param ns_dims: a list describing the dimensions of the nanostars, [len_arm, len_cen, len_end]
        :param box_dim: size of the simulation box
        :return: arms
        """
        strand_id_idx = list(strands.keys())
        # find 3' head of an arbitrary strand
        s0 = strands[min(strand_id_idx)]
        # pool all bases with their locations # now only binding points
        strand_id_idx, strands, bind_base = self.pairing_single_step(strand_id_idx, strands, s0, box_dim,
                                                                     ns_dims, True)
        # construct arm
        s1 = strands[bind_base.strand_id]
        arm_idx = 0
        arms = OrderedDict()
        arms[arm_idx] = Arm(arm_idx, ns_dims, s0.strand_id, {s0.strand_id: s0, s1.strand_id: s1})
        # loop out all other strands
        while len(strands) != 1:
            s0 = s1  # deepcopy?
            strand_id_idx.remove(s0.strand_id)
            strand_id_idx, strands, bind_base = self.pairing_single_step(strand_id_idx, strands, s0, box_dim,
                                                                         ns_dims, False)
            s1 = strands[bind_base.strand_id]
            arm_idx += 1
            arms[arm_idx] = Arm(arm_idx, ns_dims, s0.strand_id, {s0.strand_id: s0, s1.strand_id: s1})
            del strands[s0.strand_id]
        return arms

    def periodic_box_cond_correction(self, bind_base, strands: dict, box_dim):
        """
        Perform Periodic Box Condition Correction near a given base. The PBCC algorithm is similar to that implemented by oxView
        :param bind_base: the binding site of a strand for pairing two strands
        :param strands: all strands
        :param box_dim: size of the simulation box
        :return: modified strands
        """
        # Assumtions:
        # 'intended diffusion across box' only acting on continuous part, likely a whole strand.
        # 'intended diffusion across box' won't cast bases across box boundary.
        # safety check: is_diffusion. Disabled: may have 2 adjacent diffused strands.
        # dis_arr = np.array([dist(bind_base.position, base.position) if base is not bind_base else None for base in pool_base_ls])
        # dis_arr = dis_arr[dis_arr != np.array(None)]
        # assert all(dis_arr > 4) # confirming diffused strand: all greater than 5 (empirical value, binding:1~3)

        # PBC-CoM Centering # no need for base_cnt since scaling sine and cosine simultaneously does not change the angle.
        phasors = np.zeros((2, 3))  # NOT normalized to 1. Instead, noted in sines/cosines
        for strand in strands.values():
            for base in strand.base_sequence.values():
                angles = np.array(base.position) / (box_dim / (2 * np.pi))
                phasors += np.vstack((np.sin(angles), np.cos(angles)))
        CoM_pos = np.arctan2(-phasors[0], -phasors[1]) + np.pi  # Now in 0~2pi. '-' counterbalancing '+np.pi'.
        CoM_pos *= (box_dim / (2 * np.pi))
        # Update Position.
        for strand in strands.values():
            for base in strand.base_sequence.values():
                old_pos = np.array(base.position)
                new_pos = old_pos - CoM_pos + box_dim / 2  # old_pos - CoM_pos = r_CoM. shift to mid-box afterward.
                _ = base.set_position(tuple(new_pos))
                _ = self.strands[base.strand_id].base_sequence[base.base_id].set_position(tuple(new_pos))
        # PBCC
        base = bind_base
        old_pos = np.array(base.position)
        # traverse backward, skipping bind_base
        while strands[base.strand_id].base_sequence.get(base.prev_id) is not None and dist(old_pos, strands[
            base.strand_id].base_sequence[base.prev_id].position) < 4:
            base = strands[base.strand_id].base_sequence[base.prev_id]
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
        while strands[base.strand_id].base_sequence.get(base.next_id) is not None and dist(old_pos, strands[
            base.strand_id].base_sequence[base.next_id].position) < 4:
            base = strands[base.strand_id].base_sequence[base.next_id]
            old_pos = np.array(base.position)
            new_pos = ((old_pos % box_dim) + box_dim) % box_dim
            _ = base.set_position(tuple(new_pos))
            # fix self.strands as well
            _ = self.strands[base.strand_id].base_sequence[base.base_id].set_position(tuple(new_pos))
        return strands  # still, bases not in the same strand may be close to each other!
        # return self.pbc_CoM_centering(strands,box_dim) 

    def pairing_single_step(self, strand_id_idx, strands, s0, box_dim, ns_dims, is_checking_pbcc):
        """
        Single step of pairing strands for producing arms
        :param strand_id_idx: ids of all strands
        :param strands: all strands
        :param s0: the current strand
        :param box_dim: size of the simulation box
        :param ns_dims: a list describing the dimensions of the nanostars, [len_arm, len_cen, len_end]
        :param is_checking_pbcc: whether periodic boundary condition correction should be checked
        :return: strand_id_idx (updated), strands (updated), bind_base (used to identify the paired strand)
        """
        binding_location = 4  # controlling how far the binding location is from the end of arm. (4:= quadrisection point. 2:= bisection point)
        if is_checking_pbcc:
            CoM_pos = self.get_CoM_pos()
            criteria = (ns_dims[0] + ns_dims[1]) / (
                    20 - binding_location) * 8.3 * 1.3  # 20 bp : ~8.3 SU. Setting tolerance as 30% because of simul-inputs-6arms-jun_10-oxrna2/20C-0.1M-GPU, was 20%
            bp_ls = [
                list(strands[idx].base_sequence.values())[-1 - ns_dims[2] - ns_dims[0] // binding_location] if dist(
                    np.array(list(strands[idx].base_sequence.values())[
                                 -1 - ns_dims[2] - ns_dims[0] // binding_location].position),
                    CoM_pos) > criteria else None for idx in strand_id_idx]  # distances from binding point to the CoM
            bp_ls.extend([
                list(strands[idx].base_sequence.values())[ns_dims[2] + ns_dims[0] // binding_location] if dist(np.array(
                    list(strands[idx].base_sequence.values())[ns_dims[2] + ns_dims[0] // binding_location].position),
                    CoM_pos) > criteria else None for idx in strand_id_idx])
            if any(bp is not None for bp in bp_ls):
                for base in bp_ls:
                    if base is not None:
                        # send strands with largest distance to pbcc
                        strands = self.periodic_box_cond_correction(base, strands, box_dim)
                # recursion
                strand_id_idx, strands, bind_base = self.pairing_single_step(strand_id_idx, strands, s0,
                                                                             box_dim, ns_dims,
                                                                             is_checking_pbcc)  # careful -- recursion
                return strand_id_idx, strands, bind_base
        test_length = 3  # adjust if arm is not long enough. direction: center?dedicated to current NS designs.
        test_distance = 4  # paired bases cannot be more than 'test_distance' apart.
        s0_head = list(s0.base_sequence.values())[ns_dims[0] // binding_location]
        pool_base_ls = []
        for idx in strand_id_idx:
            if idx == s0.strand_id:
                continue
            bind_point_1 = list(strands[idx].base_sequence.values())[
                -1 - ns_dims[2] - ns_dims[0] // binding_location]
            bind_point_2 = list(strands[idx].base_sequence.values())[ns_dims[2] + ns_dims[0] // binding_location]
            # if is_checking_pbcc: # recursively
            #     # Centering bind_points w/o pbcc
            #     bp1_CoM_pos = np.array(bind_point_1.position) - CoM_pos
            #     bp2_CoM_pos = np.array(bind_point_2.position) - CoM_pos
            #     # Check if need PBCC. 20 bp : ~8.3 SU. Setting tolerance as 50%
            #     is_b1_pbcc = np.linalg.norm(bp1_CoM_pos) > (ns_dims[0]+ns_dims[1])/(20-binding_location)*8.3*1.5
            #     is_b2_pbcc = np.linalg.norm(bp2_CoM_pos) > (ns_dims[0]+ns_dims[1])/(20-binding_location)*8.3*1.5
            #     if is_b1_pbcc or is_b2_pbcc:
            #         # check if binding points 
            #         if is_b1_pbcc:
            #             strands = self.periodic_box_cond_correction(bind_point_1, pool_base_ls, strands, box_dim)
            #         else:
            #             strands = self.periodic_box_cond_correction(bind_point_2, pool_base_ls, strands, box_dim)
            #         strand_id_idx, strands, bind_base = self.pairing_single_step(strand_id_idx, strands, s0, box_dim, ns_dims, is_checking_pbcc) # careful -- recursion
            #         return strand_id_idx, strands, bind_base

            # selecting the pairing node, no longer the whole strand. But Not Assuming Direction.
            pool_base_ls.extend([bind_point_1,
                                 bind_point_2])
        # get the closest base
        min_dist = 1000000
        bind_base: Optional[Base] = None
        for candidate_base in pool_base_ls:
            dis = dist(candidate_base.position, s0_head.position)
            orient_dot = np.dot(np.array(s0_head.backbone), np.array(candidate_base.backbone))
            direction = np.array(candidate_base.position) - np.array(s0_head.position)
            # check if candidate locates near the direction that s0_head's backbone vector points at.
            direction_dot = np.dot(direction / np.linalg.norm(direction),
                                   s0_head.backbone)
            # and orient_dot < -0.7 and dis < 4 and direction_dot > 0.7
            if dis < min_dist and orient_dot < 0 and dis < 4 and direction_dot > 0:
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
                             strands[bind_base.strand_id].base_sequence[bind_base.base_id + i].position) for i in
                        range(test_length)]
        assert all(np.array(test_dist_ls) < test_distance)
        return strand_id_idx, strands, bind_base

    def get_CoM_pos(self):
        """
        Obtain the center of mass of the whole nanostar
        :return: CoM position in ``np.array`` .
        """
        strands = self.strands
        # CoM
        CoM_pos = np.zeros(3)
        base_cnt = 0
        for strand in strands.values():
            for base in strand.base_sequence.values():
                CoM_pos = np.add(CoM_pos, np.array(base.position))
                base_cnt += 1
        CoM_pos = np.divide(CoM_pos, base_cnt)
        return CoM_pos
