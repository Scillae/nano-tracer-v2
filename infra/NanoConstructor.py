from models import Strand, Base, NanoStar, Arm
# from models import Strand, Base, NanoStar, NanoMesh, Arm
from infra import TimeSeries
from collections import OrderedDict
from utils import save_load as SL
import numpy as np
import copy
import pickle
import os.path


class NanoConstructor:
    def __init__(self, strands_series, ns_dims, arm_number):
        """
        init NanoStar TimeSeries Constructor
        :param strands_series: TimeSeries of strands used to construct nanostar series
        :param ns_dims: topological dimensions of the nanostar to be constructed, [len_arm, len_cen, len_end]
        :param arm_number: how many arms the nanostar has
        """
        self.time_series = None
        self.strands = strands_series
        self.dim = ns_dims
        self.arm_number = arm_number  # nanomesh not using this.
        # expansion site for nanomesh - nanostars in network - as a branch in self.construct()

    def construct(self, box_dim=None):  # box_dim hacking
        """
        construct TimeSeries of nanostars (can be branched to create other objects)
        :param box_dim: size of the simulation box, required if performing periodic boundary condition correction
        :return: nanostar series
        """
        self.time_series = TimeSeries(box_dim=self.strands.params['box_dim'], ns_dims=self.dim,
                                      arm_number=self.arm_number)

        for t_stamp in self.strands:
            ns = NanoStar(self.strands[t_stamp], self.dim, self.arm_number, box_dim=box_dim)
            self.time_series[t_stamp] = ns
        return self.time_series
