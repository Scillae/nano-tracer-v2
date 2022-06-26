from collections import OrderedDict


class TimeSeries(OrderedDict):
    def __init__(self, keys: list = None, vals: list = None):
        super(TimeSeries, self).__init__()
        # check keys/vals
        if keys and vals and len(set(keys)) == len(vals) and len(keys) == len(vals):
            self.add_instances(keys, vals)

    def add_instances(self, keys: list, vals: list):
        # log duplicates
        if not set(self.keys()).isdisjoint(vals):
            print("Duplicated keys")
        self.update(zip(keys, vals))
