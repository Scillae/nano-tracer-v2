from collections import OrderedDict


class TimeSeries(OrderedDict):
    def __init__(self, keys:list = None, vals:list = None, **kwargs):
        '''
        TimeSeries is an customized OrderedDict in which parameters can be attatched.
        :param keys: list of keys
        :param vals: list of values
        :param kwargs: will be saved as parameters
        Note that the two lists should be matched!
        '''
        super(TimeSeries, self).__init__()
        # check keys/vals. if existing, try to fill in
        if keys and vals:
            self.add_instances(keys, vals)
        # keyword parameters saved in property: params (dictionary)
        self.params = kwargs
    
    def add_instances(keys: list, vals: list):
        '''
        Add instance(s) into the current TimeSeries. Duplicated Keys will cause overwritting.
        :param keys: list of keys
        :param vals: list of values
        Note that the two lists should be matched!
        '''
        if not set(keys).isdisjoint(self.keys):
            print("Duplicated Keys! Proceeding, but OVERWRITTING with new values.")
        if len(kwargs['keys']) == len(kwargs['vals']):
            self.update(zip(keys, vals))
        else:
            raise Exception("Keys are not matched with Vals! Cannot Update TimeSeries Dictionary.")
    
