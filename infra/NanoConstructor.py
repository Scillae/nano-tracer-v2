from models import Strand, Base, NanoStar, NanoMesh, Arm, TimeMachine
from collections import OrderedDict
from utils.tools import save_load as SL
import numpy as np
import copy
import pickle
import os.path

class NanoConstructor:
    def __init__(self, strands_series, dims_ls, arm_num):
        self.time_machine = None
        self.strands = strands_series
        self.dim = dims_ls
        self.arm_num = arm_num # nanomesh not using this.
        # expansion site for nanomesh - nanostars in network - as a branch in self.construct()

    def construct(self, box_dim = None): # box_dim hacking
        self.time_machine = TimeMachine()
        for t_stamp in self.strands.timeseries:
            ns = NanoStar(self.strands.time_capsule[t_stamp],self.dim, self.arm_num, box_dim=box_dim)
            self.time_machine.add_instance(t_stamp, ns)
        return self.time_machine


