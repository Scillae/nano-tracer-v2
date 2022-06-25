from models import Strand, Base, NanoStar, NanoMesh, Arm, TimeMachine
from collections import OrderedDict
from utils import save_load as SL
import numpy as np
import copy
import pickle
import os.path


class NanoConstructor:
    def __init__(self, strands_tm, dims_ls, arm_num):
        """
        init NanoStar?
        :param strands_tm:
        :param dims_ls:
        :param arm_num:
        """
        self.time_machine = None
        self.dim = dims_ls
        self.strands = strands_tm
        self.arm_num = arm_num  # nanomesh not using this.

    def construct(self, box_dim=None):  # box_dim hacking
        """
        construct TimeMachine
        :param box_dim:
        :return: time machine
        """
        self.time_machine = TimeMachine()
        self.time_machine.box_dim = box_dim  # box_dim. A hacking solution.
        for t_stamp in self.strands.timeseries:
            ns = NanoStar(self.strands.time_capsule[t_stamp], self.dim, self.arm_num, box_dim=box_dim)
            self.time_machine.add_instance(t_stamp, ns)
        return self.time_machine
