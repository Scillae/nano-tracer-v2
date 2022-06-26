from collections import OrderedDict

from models import TimeSeries


class Strand:
    def __init__(self, strand_id, base_seq=None, timestamp=None):
        """
        init Strand
        :param strand_id:
        :param base_seq:
        :param timestamp:
        """
        self.strand_id = strand_id
        self.timestamp = timestamp
        self.base_sequence = OrderedDict({}) if base_seq is None else base_seq

    def add_base(self, base_entity):
        """
        add base to strand
        :param base_entity:
        :return:
        """
        self.base_sequence[base_entity.base_id] = base_entity

    @staticmethod
    def parse(s):
        pass

