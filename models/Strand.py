from collections import OrderedDict
class Strand:
    def __init__(self, strand_id, base_seq = None, timestamp = None):
        self.strand_id = strand_id
        self.timestamp = timestamp
        self.base_sequence = OrderedDict({}) if base_seq is None else base_seq

    def add_base(self, base_entity):
        self.base_sequence[base_entity.base_id] = base_entity

    @staticmethod
    def parse(s):
        pass