import h5py

class HDF5Data:
    def __init__(self):
        self.data = {}

    def load_from_hdf5(self, file_path, paths=None):
        """
        Load specific parts of an HDF5 file.

        :param file_path: Path to the HDF5 file.
        :param paths: List of paths to load. If None, load the entire file.
        """
        with h5py.File(file_path, 'r') as hdf_file:
            if paths is None:
                # Load the entire file
                self._load_group(hdf_file, self.data)
            else:
                # Load only specified paths
                for path in paths:
                    self._load_path(hdf_file, path)

    def _load_group(self, group, data_dict):
        for key, item in group.items():
            if isinstance(item, h5py.Dataset):
                data_dict[key] = item[()]
            elif isinstance(item, h5py.Group):
                data_dict[key] = {}
                self._load_group(item, data_dict[key])

    def _load_path(self, hdf_file, path):
        """
        Load a specific path from the HDF5 file.

        :param hdf_file: Open HDF5 file object.
        :param path: Path to load.
        """
        try:
            item = hdf_file[path]
            keys = path.strip('/').split('/')
            data_dict = self.data
            for key in keys[:-1]:
                data_dict = data_dict.setdefault(key, {})
            if isinstance(item, h5py.Dataset):
                data_dict[keys[-1]] = item[()]
            elif isinstance(item, h5py.Group):
                data_dict[keys[-1]] = {}
                self._load_group(item, data_dict[keys[-1]])
        except KeyError:
            print(f"Path {path} not found in the HDF5 file.")

    def __repr__(self):
        return repr(self.data)

# Usage
file_path = 'USAXS.h5'
hdf5_data = HDF5Data()

# Load specific paths
paths_to_load = ['/entry/metadata', '/entry/flyScan']
hdf5_data.load_from_hdf5(file_path, paths=paths_to_load)

# Now you can access the loaded data
print(hdf5_data)