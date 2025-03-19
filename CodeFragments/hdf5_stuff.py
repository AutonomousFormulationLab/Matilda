
# append Q ands Int to hdf5 
# Example NumPy arrays
Intensity = np.array([1, 2, 3, 4, 5])
Qvector = np.array([10, 20, 30, 40, 50])

# Open or create an HDF5 file
with h5py.File('data.h5', 'a') as hdf_file:
    # Create or access the DisplayData group
    if 'DisplayData' not in hdf_file:
        display_data_group = hdf_file.create_group('DisplayData')
    else:
        display_data_group = hdf_file['DisplayData']
    
    # Append Intensity array
    if 'Intensity' in display_data_group:
        # Resize the existing dataset to accommodate new data
        existing_data = display_data_group['Intensity']
        existing_data.resize((existing_data.shape[0] + Intensity.shape[0],))
        existing_data[-Intensity.shape[0]:] = Intensity
    else:
        # Create a new dataset for Intensity
        display_data_group.create_dataset('Intensity', data=Intensity, maxshape=(None,))
    
    # Append Qvector array
    if 'Qvector' in display_data_group:
        # Resize the existing dataset to accommodate new data
        existing_data = display_data_group['Qvector']
        existing_data.resize((existing_data.shape[0] + Qvector.shape[0],))
        existing_data[-Qvector.shape[0]:] = Qvector
    else:
        # Create a new dataset for Qvector
        display_data_group.create_dataset('Qvector', data=Qvector, maxshape=(None,))

