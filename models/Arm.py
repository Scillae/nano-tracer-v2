import numpy as np
from collections import OrderedDict


class Arm:
    def __init__(self, arm_id, dims_ls, leading_strid, strands):
        """
        init Arm
        :param arm_id: number id of the arm
        :param dims_ls: ?
        :param leading_strid:
        :param strands: list of strands
        """
        self.arm_id = arm_id
        assert len(strands) == 2
        # Noella: Accepted_sci
        # self.strand_id_0, self.strand_id_1 = [strand.strand_id for strand in strands]
        self.strand_id_0, self.strand_id_1 = [strand for strand in strands]
        # Noella: Questioned_sci(True. 2*arm.)
        # Maybe Scillae want to compare length of strand.base_sequence with arm+cen+end, but they are not equal in given example
        # dims_ls = 20 + 2 + 7 = 29, len(strands[self.strand_id_0].base_sequence) = 49
        # Or (arm + arm + center + end) ? 
        # assert np.sum(np.array(dims_ls)) == len(strands[self.strand_id_0]) and np.sum(np.array(dims_ls)) == len(strands[self.strand_id_1])
        self.base_pairs, self.single_end = self.pair(dims_ls, leading_strid, strands)
        # TODO how to present centers

    # pair bases
    def pair(self, dims_ls, leading_strid, strands):
        len_arm, len_cen, len_end = dims_ls
        s0 = list(strands[leading_strid].base_sequence.values())
        del strands[leading_strid]
        s1 = list(list(strands.values())[0].base_sequence.values())
        s1.reverse()  # reverse s1 so that the single_end is at the start.

        pair_tp_dic = OrderedDict()  # {pair_id:(pair_base0, pair_base1)}
        single_end_dic = OrderedDict()  # {single_end_id:single_end_base}
        for i in range(len_end):
            # index starts from center, beginning from 1.
            single_end_dic[len_end - i] = s1[i]
        del s1[0:len_end]
        for i in range(len_arm):
            # index starts from center, beginning from 1.
            pair_tp_dic[len_arm - i] = (s0[i], s1[i])
        return pair_tp_dic, single_end_dic

        '''
        #single-end recog. TODO:'=='
        len_arm, len_cen, len_end = dims_ls
        s0 = strands[self.strand_id_0].base_sequence.items() # may use values() if not desiring absolute locations
        s1 = strands[self.strand_id_1].base_sequence.items()
        s0_tp_ls = list(s0)[:len_end]
        s0_bt_ls = list(s0)[len(s0)-len_end:]
        s1_tp_ls = list(s1)[:len_end]
        s1_bt_ls = list(s1)[len(s1)-len_end:]
        single_end = []
        s0_flip_flag = False
        s1_flip_flag = False
        s0_end = False
        s0e = s1e = s0bbs = s1bbs = 0  #base_before_singleend, base_other_end
        if(s0_tp_ls == s1_tp_ls):
            single_end = s0_tp_ls
            _, s0bbs = s0[len_end]
            _, s1bbs = s1[len_end]            
            _, s0e = s0[len(s0)-1]
            _, s1e = s1[len(s1)-1]
            if dist(s0e.position, s1bbs.position) < dist(s1e.position, s0bbs.position):# s0[48] is head, s0 needs flip
                s0_flip_flag = True
            else:# s1[48] is head, s1 needs flip
                s1_flip_flag = True
                s0_end = True
        elif(s0_bt_ls == s1_bt_ls):
            single_end = s0_bt_ls
            _, s0bbs = s0[len(s0)-1-len_end]
            _, s1bbs = s1[len(s0)-1-len_end]            
            _, s0e = s0[0]
            _, s1e = s1[0]
            if dist(s0e.position, s1bbs.position) < dist(s1e.position, s0bbs.position):# s0[0] is head, s1 needs flip
                s1_flip_flag = True
            else:# s1[0] is head, s0 needs flip
                s0_flip_flag = True
                s0_end = True
        elif(s0_tp_ls == s1_bt_ls):
            single_end = s0_tp_ls
            _, s0bbs = s0[len_end]
            _, s1bbs = s1[len(s0)-1-len_end]             
            _, s0e = s0[len(s0)-1]
            _, s1e = s1[0]
            if dist(s0e.position, s1bbs.position) < dist(s1e.position, s0bbs.position):# s0[48] is head, both need flip
                s1_flip_flag = True
                s0_flip_flag = True
            else:# s1[0] is head, no flip
                s0_end = True
        elif(s0_bt_ls == s1_tp_ls):
            single_end = s0_bt_ls
            _, s0bbs = s0[len(s0)-1-len_end]
            _, s1bbs = s1[len_end]             
            _, s0e = s0[0]
            _, s1e = s1[len(s1)-1]
            if dist(s0e.position, s1bbs.position) < dist(s1e.position, s0bbs.position):# s0[0] is head, no flip
                assert 0 == 0
            else:# s1[48] is head, both need flip
                s1_flip_flag = True
                s0_flip_flag = True
                s0_end = True
        else:
            assert 0 == 1
        
        pairs_tp_tp_ls = []
        if s0_flip_flag:
            s0.reverse()    #may not use this if wish to preserve sequence location.
        if s1_flip_flag:
            s1.reverse()
        if s0_end:
            for i in range(len_arm):
                pairs_tp_tp_ls.append((s0[len_end+i],s1[i]))
        else:
            for i in range(len_arm):
                pairs_tp_tp_ls.append((s0[i],s1[len_end+i]))
        
        return pairs_tp_tp_ls, single_end, s0_end



        if dist(s0e.position, s1bbs.position) < dist(s1e.position, s0bbs.position):
            s0_leading_flag = True
        else:
            s0_leading_flag = False
#arm=3, end=2: skip 2
        pared_bases_tp_ls = []
        if reverse_flag: #==False
            if s0_leading_flag:
                for i in range(len_arm-1):
                    pared_bases_tp = (i,s0[len_arm-1-i],s1[len(s1)-len_end-len_arm+i])
                    pared_bases_tp_ls.append(pared_bases_tp)
            else:
                for i in range(len_arm-1)
                    pared_bases_tp = (i,s0[len(s1)-len_end-len_arm+i],s1[len_arm-1-i])
                    pared_bases_tp_ls.append(pared_bases_tp)
        else:
            if s0_leading_flag:
                for i in range(len_arm-1):
                    pared_bases_tp = (i,s0[len_arm-1-i],s1[len(s1)-len_end-len_arm+i])
                    pared_bases_tp_ls.append(pared_bases_tp)
            else:
                for i in range(len_arm-1)
                    pared_bases_tp = (i,s0[len(s1)-len_end-len_arm+i],s1[len_arm-1-i])
                    pared_bases_tp_ls.append(pared_bases_tp)

'''


def dist(t1, t2):
    return np.sqrt(np.sum(np.square(np.array(t1) - np.array(t2))))
