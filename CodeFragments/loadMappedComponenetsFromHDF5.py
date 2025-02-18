
import h5py

class HDF5Data:
    def __init__(self):
        # Initialize attributes that will be populated from the HDF5 file
        self.mac1 = None            #array
        self.FS_ScanTime = None     #variable
        self.name = None            #string

    def load_from_hdf5(self, file_path, mapping):
        """
        Load specific parts of an HDF5 file into class attributes based on a mapping.

        :param file_path: Path to the HDF5 file.
        :param mapping: Dictionary mapping class attributes to HDF5 paths.
        """
        with h5py.File(file_path, 'r') as hdf_file:
            for attr, path in mapping.items():
                try:
                    item = hdf_file[path]
                    if isinstance(item, h5py.Dataset):
                        setattr(self, attr, item[()])
                    else:
                        print(f"Path {path} is not a dataset.")
                except KeyError:
                    print(f"Path {path} not found in the HDF5 file.")

    def __repr__(self):
        return f"HDF5Data(mca1={self.mca1}, FS_ScanTime={self.FS_ScanTime}, name={self.name})"

# Usage
file_path = 'USAXS.h5'
RAW = HDF5Data()

# Define a mapping from class attributes to HDF5 paths
mapping = {
    'mca1': '/entry/flyScan/mca1',
    'FS_ScanTime': '/entry/flyScan/FS_ScanTime',
    'name': '/entry/instrument/name'
}

# Load data based on the mapping
RAW.load_from_hdf5(file_path, mapping)

# Now you can access the loaded data as class attributes
print(RAW)