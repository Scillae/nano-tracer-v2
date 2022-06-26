from collections import OrderedDict


class Strand:
    def __init__(self, strand_id:int, base_seq:OrderedDict=None, timestamp:int=None):
        """
        Strand is a series of bases(nucleotides). Whether the bases are paired remain unknown.
        :param strand_id: id of this strand
        :param base_seq: the series of bases
        :param timestamp: at what time the strand instance appears (for debug use)
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
