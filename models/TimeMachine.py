import pickle
import os.path
from collections import OrderedDict
from models import Strand, Base, Arm

class TimeMachine:
    def __init__(self, path="data/cache.pkl"):
        self.timeseries = []
        self.path = path
        # TODO: sort
        self.cache_writer = None
        self.cache_reader = None
        self.mode = ''
        self.time_capsule = OrderedDict()
        self.box_dim = None # box_dim hacking
            
    # def save(self):
    #   if self.cache_writer is not None:
    #     self.cache_writer.close()
    #   self.mode = ''

    # def set_write_mode(self):
    #   if self.cache_reader is not None:
    #     self.cache_reader.close()
    #   if self.cache_writer is None:
    #     self.cache_writer = open(self.path, 'wb')
    #   self.mode = 'w'
    # def set_read_mode(self):
    #   if self.cache_writer is not None:
    #     self.cache_writer.close()
    #   if self.cache_reader is None:
    #     self.cache_reader = open(self.path, 'rb')
    #   self.mode = 'r'

    def add_instance(self, timestamp, instance):
        if timestamp == -1:
            raise Exception('Error reading timestamp, get value -1')
        if timestamp in self.timeseries:
            print('Duplicate timestamp detected')
        else:
            self.timeseries.append(timestamp)
        self.time_capsule[timestamp] = instance

    def get_instance_by_time(self, timestamp):
        return self.time_capsule.get(timestamp)

    # def get_strand_by_time(self, strand_id, timestamp):
    #     '''
    #     '''
    #     if timestamp not in self.timeseries:
    #         raise Exception('Timestamp not exists')
    #     # TODO: handle exception strand_id does not exist
    #     return self.time_capsule[timestamp].get(strand_id)
