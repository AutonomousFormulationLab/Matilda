'''
    this contains needed hdf5 support for matilda
    created by Argo on 2025-04-05, not tested yet. 

'''
import h5py
import numpy as np

def save_dict_to_hdf5(dic, location, h5file):
    """
    Save a dictionary to an HDF5 file.

    Parameters:
    dic (dict): The dictionary to save.
    filename (str): The name of the HDF5 file.
    """
    def recursively_save_dict_contents_to_group(h5file, path, dic):
        for key, item in dic.items():
            if isinstance(item, dict):
                # Create a new group for nested dictionaries
                group = h5file.create_group(path + key)
                recursively_save_dict_contents_to_group(h5file, path + key + '/', item)
            else:
                # Save numpy arrays and other data types
                h5file[path + key] = item

    recursively_save_dict_contents_to_group(h5file, location, dic)

# # Example usage
# data_dict = {
#     'array': np.array([1, 2, 3]),
#     'value': 42,
#     'nested': {
#         'string': 'hello',
#         'array2': np.array([4, 5, 6])
#     }
# }

#save_dict_to_hdf5(data_dict, 'data.h5')

def load_dict_from_hdf5(hdf_file, location):
    """
    Load a dictionary from an HDF5 file.

    Parameters:
    filename (str): The name of the HDF5 file.

    Returns:
    dict: The loaded dictionary.

    location (str): The location in the HDF5 file to load from.
    'root:DisplayData/'
    """
    def recursively_load_dict_contents_from_group(h5file, path):
        ans = {}
        for key, item in h5file[path].items():
            if isinstance(item, h5py._hl.group.Group):
                ans[key] = recursively_load_dict_contents_from_group(h5file, path + key + '/')
            else:
                ans[key] = item[()]
        return ans

    return recursively_load_dict_contents_from_group(hdf_file,location)



