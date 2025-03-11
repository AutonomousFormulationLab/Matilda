# these are support functions for Matilda

import h5py
import numpy as np
import logging

## support stuff here
# Function to recursively read a group and store its datasets in a dictionary
def read_group_to_dict(group):
    data_dict = {}
    for key, item in group.items():
        if isinstance(item, h5py.Dataset):
            # Read the dataset
            data = item[()]
            # Check if the dataset is an array with a single element
            if data.size == 1:
                # Convert to a scalar (number or string)
                data = data.item()
            data_dict[key] = data
        elif isinstance(item, h5py.Group):
            # If the item is a group, recursively read its contents
            data_dict[key] = read_group_to_dict(item)
    return data_dict
