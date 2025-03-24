# these are support functions for Matilda

import h5py
import numpy as np
import logging
import copy
import xarray as xr
## support stuff here
# Function to recursively read a group and store its datasets in a dictionary
def read_group_to_dict(group):
    data_dict = {}
    for key, item in group.items():
        if isinstance(item, h5py.Dataset):
            # Read the dataset
            data = item[()]
             # Check if the dataset is bytes
            if isinstance(data, bytes):
                # Decode bytes to string
                data = data.decode('utf-8')
            # Check if the dataset is an array with a single element
            elif hasattr(data, 'size') and data.size == 1:
                # Convert to a scalar (number or string)
                data = data.item()
            data_dict[key] = data
        elif isinstance(item, h5py.Group):
            # If the item is a group, recursively read its contents
            data_dict[key] = read_group_to_dict(item)
    return data_dict

def results_to_dataset(results):
    results = copy.deepcopy(results)
    ds = xr.Dataset()
    ds['USAXS_int'] = ('q',results['ReducedData']['UPD'])
    ds['q'] = results['ReducedData']['Q_array']
    del results['ReducedData']['UPD']
    del results['ReducedData']['Q_array']
    ds.update(results['ReducedData'])
    for our_name,raw_name in [('AR_angle','ARangles'),
                              ('TimePerPoint','TimePerPoint'),
                              ('Monitor','Monitor'),
                              ('UPD','UPD_array'),
                             ]:
        ds[our_name] = ('flyscan_bin',results['RawData'][raw_name])
        del results['RawData'][raw_name]
    for our_name,raw_name in [('AmpGain','AmpGain'),
                              ('AmpReqGain','AmpReqGain'),
                              ('amp_change_channel','Channel')
                             ]:
        ds[our_name] = ('amp_change_channel',results['RawData'][raw_name])
        del results['RawData'][raw_name]
                                      
    ds.attrs.update(results['RawData']['metadata'])
    del results['RawData']['metadata']
    ds.attrs['instrument'] = results['RawData']['instrument']
    del results['RawData']['instrument']
    ds.update(results['RawData'])

    return ds

