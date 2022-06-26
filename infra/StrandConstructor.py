import os
import pickle
import os.path
from models import Strand, Base
from infra import TimeSeries
from utils import assignment_parser, nextline, formatter
from collections import OrderedDict

number_to_base = {0: 'A', 1: 'G', 2: 'C', 3: 'T'}

base_to_number = {'A': 0, 'a': 0, 'G': 1, 'g': 1,
                  'C': 2, 'c': 2, 'T': 3, 't': 3,
                  'U': 3, 'u': 3, 'D': 4}


class StrandConstructor:
    """
    Input topology and configuration file, return a time series of Strands
    This Constructor is currently INCOMPATIBLE with RNA!
    """

    def __init__(self, top_file: str, traj_file: str):
        """
        init reader with data file paths
        :param top_file: path to top file
        :param traj_file: path to traj file
        """
        self.top_cursor = open(top_file)
        self.traj_cursor = open(traj_file)
        self.format = {
            'timestamp': (int,),
            'box': tuple([float for i in range(3)]),
            'energy': tuple([float for i in range(3)]),
            'nucleotide_tr': tuple(
                [tuple([float for i in range(3)]) for i in range(5)]),
            'nucleotide_tp': tuple(
                [int, lambda x: base_to_number[x], int, int]),
            # TODO handle undefined base
            'counts': (int, int),
        }
        self.strands = OrderedDict()
        self.timestamp = -1
        self.time_series = None

    def read_single_strand(self) -> int:
        """
        read one strand to self.strands
        :return: status code, where -1 and -2 indicate line missing from traj and top file respectively, otherwise presents the base amount.
        """
        # initialize strand
        strand = None

        # nucleotide count and strand count
        line = nextline(self.top_cursor)
        if line is None:
            self.top_cursor.seek(0)
            return -2
        nucleotide_cnt, strand_cnt = formatter(self.format['counts'], line)

        line = nextline(self.traj_cursor)
        if line is None:
            return -1
        # timestamp
        self.timestamp = timestamp = formatter(self.format['timestamp'],
                                               assignment_parser(line)[-1])[0]
        # box size
        line = nextline(self.traj_cursor)
        box = formatter(self.format['box'], assignment_parser(line)[-1])

        # info
        line = nextline(self.traj_cursor)
        total, potential, kinetic = formatter(self.format['energy'],
                                              assignment_parser(line)[-1])
        # nucleotide
        base_incre = 0
        tr_line = nextline(self.traj_cursor)
        tp_line = nextline(self.top_cursor)
        line_counter = 0
        while tr_line and tp_line:
            tr_params = formatter(self.format['nucleotide_tr'], tr_line)
            tp_params = formatter(self.format['nucleotide_tp'], tp_line)
            params = tr_params + (base_incre,) + tp_params
            base_incre += 1

            # get strand from dict by id, create new strand if not exists
            strand_id = tp_params[0]

            if self.strands.get(strand_id) is None:
                strand = self.strands[strand_id] = Strand(
                    strand_id, timestamp=timestamp)
            if strand is None or strand.strand_id != strand_id:
                strand = self.strands[strand_id]

            strand.add_base(Base.parse_list(params))
            if base_incre < nucleotide_cnt:
                tr_line = nextline(self.traj_cursor)
                tp_line = nextline(self.top_cursor)
            else:
                break

        return base_incre

    def read_data(self) -> TimeSeries:
        """
        read all strands in specified files
        :return:
        """
        self.time_series = TimeSeries()
        res = self.read_single_strand()
        while res != -1:
            if res == -2:
                # reverse the read strands
                strands_r = OrderedDict()
                for strand_id, strand in self.strands.items():
                    base_seq_r_ls = list(strand.base_sequence.items())
                    base_seq_r_ls.reverse()
                    od = OrderedDict(base_seq_r_ls)
                    strands_r[strand_id] = Strand(strand_id, od)
                self.strands = strands_r
                self.time_series[self.timestamp] = self.strands
                self.strands = OrderedDict()
            else:
                print(f'  Strands: {len(self.strands)}, timestamp: {list(self.strands.values())[0].timestamp}')
                for i in self.strands.values():
                    print(
                        f'    id: {i.strand_id}, base_count: {len(i.base_sequence)}')
            res = self.read_single_strand()
        print('Reach end of file.')
        print(f'Total time series length: {len(self.time_series)}')
        return self.time_series
