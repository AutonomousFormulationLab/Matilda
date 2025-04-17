import h5py
import numpy as np


def flatten_dict(d, parent_key='', sep='/'):
    """Flatten a nested dictionary."""
    items = [] 
    for k, v in d.items():
        new_key = f"{parent_key}{sep}{k}" if parent_key else k
        if isinstance(v, dict):
            items.extend(flatten_dict(v, new_key, sep=sep).items())
        else:
            items.append((new_key, v))
    return dict(items)

def unflatten_dict(d, sep='/'):
    """Unflatten a dictionary."""
    result_dict = {}
    for k, v in d.items():
        keys = k.split(sep)
        d = result_dict
        for key in keys[:-1]:
            d = d.setdefault(key, {})
        d[keys[-1]] = v
    return result_dict

def save_dict_to_hdf5(d, filename, group_name):
    """Save a dictionary to a specific group in an HDF5 file."""
    flat_dict = flatten_dict(d)
    with h5py.File(filename, 'a') as h5file:  # Use 'a' to append to the file
        group = h5file.require_group(group_name)
        for key, value in flat_dict.items():
            # Convert value to a numpy array if it's not already
            if isinstance(value, str):
                # Encode strings to bytes
                value = np.string_(value)
            elif not isinstance(value, np.ndarray):
                value = np.array(value)
            group.create_dataset(key, data=value)

def load_dict_from_hdf5(filename, group_name):
    """Load a dictionary from a specific group in an HDF5 file."""
    with h5py.File(filename, 'r') as h5file:
        group = h5file[group_name]
        flat_dict = {}
        for key in group.keys():
            value = group[key][()]
            # Decode bytes to strings if necessary
            if isinstance(value, bytes):
                value = value.decode('utf-8')
            flat_dict[key] = value
    return unflatten_dict(flat_dict)

# # Example usage
# nested_dict = {
#     'level1': {
#         'level2': {
#             'value1': 1,
#             'value2': [1, 2, 3]
#         },
#         'level2b': {
#             'value3': 3.14
#         }
#     },
#     'level1b': {
#         'value4': 'hello'
#     }
# }

# # Save the dictionary to a specific group
# save_dict_to_hdf5(nested_dict, 'nested_dict.h5', 'my_group')

# # Load the dictionary from the specific group
# loaded_dict = load_dict_from_hdf5('nested_dict.h5', 'my_group')
# print(loaded_dict)
