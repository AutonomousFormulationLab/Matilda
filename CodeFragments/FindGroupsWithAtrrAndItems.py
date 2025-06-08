import h5py

def find_matching_groups(hdf5_file_path, required_attributes, required_items):
    def check_group(name, obj):
        if isinstance(obj, h5py.Group):
            print(name)
            # Check attributes
            attributes_match = all(
                attr in obj.attrs and obj.attrs[attr] == value
                for attr, value in required_attributes.items()
            )
            
            # Check items
            items_match = True
            for item, expected_value in required_items.items():
                if item in obj:
                    actual_value = obj[item][()]
                    # Decode byte strings to regular strings if necessary
                    if isinstance(actual_value, bytes):
                        actual_value = actual_value.decode('utf-8')
                    if actual_value != expected_value:
                        items_match = False
                        break
                else:
                    items_match = False
                    break
            
            if attributes_match and items_match:
                matching_group_paths.append(name)

    matching_group_paths = []

    with h5py.File(hdf5_file_path, 'r') as hdf5_file:
        hdf5_file.visititems(check_group)

    return matching_group_paths


# Example usage:
hdf5_file_path = 'C:/Users/ilavsky/Desktop/TestNexus.hdf'
required_attributes = {'canSAS_class': 'SASentry', 'NX_class': 'NXsubentry'}
required_items = {'definition': 'NXcanSAS'}

matching_groups = find_matching_groups(hdf5_file_path, required_attributes, required_items)
print(matching_groups)
