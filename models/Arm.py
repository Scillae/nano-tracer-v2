import numpy as np
from collections import OrderedDict


class Arm:
    def __init__(self, arm_id: int, ns_dims: list, leading_strid: int, strands: dict):
        """
        Arm is the dsDNA part of a nanostar
        :param arm_id: number id of the arm
        :param ns_dims: a list describing the dimensions of the nanostars, [len_arm, len_cen, len_end]
        :param leading_strid: indicating the leading strand of the arm (w/o ssDNA end)
        :param strands: list of strands
        """
        self.arm_id = arm_id
        assert len(strands) == 2
        self.strand_id_0, self.strand_id_1 = [strand_id for strand_id in strands]
        self.base_pairs, self.single_end = self.pair(ns_dims, leading_strid, strands)

    # pair bases
    def pair(self, ns_dims: list, leading_strid: int, strands: dict):
        """
        Pair bases in the given two strands to form an arm. Bps in arm count from center, starting from 1
        :param ns_dims: a list describing the dimensions of the nanostars, [len_arm, len_cen, len_end]
        :param leading_strid: indicating the leading strand of the arm (w/o ssDNA end)
        :param strands: list of strands
        """
        len_arm, len_cen, len_end = ns_dims
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


def dist(t1, t2):
    return np.sqrt(np.sum(np.square(np.array(t1) - np.array(t2))))
