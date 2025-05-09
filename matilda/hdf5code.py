'''
    this contains needed hdf5 support for matilda
    used by saving blank BL_QRS data. 
    use code: 
    save_dict_to_hdf5(data_dict, 'data.h5') to save a dictionary to an hdf5 file.
    use code: load_dict_from_hdf5('data.h5') to load a dictionary from an hdf5 file.
    use code: load_dict_from_hdf5('data.h5','root:DisplayData/') to load a dictionary from an hdf5 file.

'''
import h5py
import numpy as np
import six  #what is this for???
import datetime

def saveNXcanSAS(Sample,path, filename):
    
    #read stuff from the data dictionary
    Intensity = Sample["CalibratedData"]["Intensity"]
    Q = Sample["CalibratedData"]["Q"]
    Error = Sample["CalibratedData"]["Error"]
    dQ = Sample["CalibratedData"]["dQ"]
    units = Sample["CalibratedData"]["units"]
    Kfactor = Sample["CalibratedData"]["Kfactor"]
    OmegaFactor = Sample["CalibratedData"]["OmegaFactor"]
    BlankName = Sample["CalibratedData"]["BlankName"]
    thickness = Sample["CalibratedData"]["thickness"]
    label = Sample["RawData"]["Filename"]
    timeStamp = Sample["RawData"]["metadata"]["timeStamp"]
    SampleName = Sample["RawData"]["sample"]["name"]

    # SMR_Int =Sample["CalibratedData"]["SMR_Int"]
    # SMR_Error =Sample["CalibratedData"]["SMR_Error"]
    # SMR_Qvec =Sample["CalibratedData"]["SMR_Qvec"]
    # SMR_dQ =Sample["CalibratedData"]["SMR_dQ"]

    # create the HDF5 NeXus file
    with h5py.File(path+'/'+filename, "w") as f:
        # point to the default data to be plotted
        f.attrs['default']          = SampleName
        # give the HDF5 root some more attributes
        f.attrs['file_name']        = filename
        f.attrs['file_time']        = timeStamp 
        f.attrs['instrument']       = '12IDE USAXS'
        f.attrs['creator']          = 'Matilda NeXus writer'
        f.attrs['NeXus_version']    = '4.3.0' #2025-5-9 4.3.0 is rc, it is current. 
        f.attrs['HDF5_version']     = six.u(h5py.version.hdf5_version)
        f.attrs['h5py_version']     = six.u(h5py.version.version)

        # create the NXentry group
        nxentry = f.create_group(SampleName)
        nxentry.attrs['NX_class'] = 'NXentry'
        nxentry.attrs['canSAS_class'] = 'SASentry'
        nxentry.attrs['default'] = 'sasdata'
        nxentry.attrs['title'] = SampleName
        
        # create the NXdata group for I(Q) for the avergaed data
        nxdata = nxentry.create_group('sasdata')
        nxdata.attrs['NX_class'] = 'NXdata'
        nxdata.attrs['canSAS_class'] = 'SASdata'
        nxdata.attrs['signal'] = 'I'      # Y axis of default plot
        nxdata.attrs['I_axes'] = 'Q'      # X axis of default plot
        #nxdata.attrs['Q_indices'] = [1]    # TODO not sure what this means

        # Y axis data
        ds = nxdata.create_dataset('I', data=Intensity)
        ds.attrs['units'] = '1/cm'
        ds.attrs['uncertainties'] = 'Idev'
        ds.attrs['long_name'] = 'cm2/cm3'    # suggested X axis plot label

        # X axis data
        ds = nxdata.create_dataset('Q', data=Q)
        ds.attrs['units'] = '1/angstrom'
        ds.attrs['long_name'] = 'Q (A^-1)'    # suggested Y axis plot label
        ds.attrs['resolutions'] = 'Qdev'
       
        # d X axis data
        ds = nxdata.create_dataset('Qdev', data=dQ)
        ds.attrs['units'] = '1/angstrom'
        ds.attrs['long_name'] = 'Q (A^-1)'   
        # dI axis data
        ds = nxdata.create_dataset('Idev', data=Error)
        ds.attrs['units'] = 'cm2/cm3'
        ds.attrs['long_name'] = 'Uncertainties'  

        # create the NXinstrument metadata group
        # nxinstr = nxentry.create_group('instrument')
        # nxinstr.attrs['NX_class'] = 'NXinstrument'
        # nxinstr.attrs['canSAS_class'] = 'SASinstrument'

        # nxprocess = nxinstr.create_group('simulation_engine')
        # nxprocess.attrs['NX_class'] = 'NXprocess'
        # nxprocess.attrs['canSAS_class'] = 'SASprocess'
        # nxprocess.attrs['name'] = '12IDE USAXS'
        # nxprocess.attrs['date'] = timestamp # @TODO: get timestamp from simulation run and embed here.
        # nxprocess.attrs['description'] = 'USAXS or SWAXS data'

        # sim_notes = nxprocess.create_group('NOTE')
        # sim_notes.attrs['NX_class'] = 'NXnote'

        # sim_notes.attrs['description'] = 'UUSAXS or SWAXS data'
        # sim_notes.attrs['author'] = '12IDE USAXS/SAXS/WAXS'
        # sim_notes.attrs['data'] = 'TBA' #@TODO

        # nxsample = nxentry.create_group('sample')
        # nxsample.attrs['NX_class'] = 'NXsample'
        # nxsample.attrs['canSAS_class'] = 'SASsample'
        
        # nxsample.attrs['name'] = SampleName
        # nxsample.attrs['description'] = 'SAMPLE DESCRIPTION GOES HERE'
        # #nxsample.attrs['type'] = 'simulated data'

        #comp.create_dataset
        
    print("wrote file:", filename)





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



