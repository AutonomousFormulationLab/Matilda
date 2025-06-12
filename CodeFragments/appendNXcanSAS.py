'''
this is place to add code for writing 1D NXcanSAS data
This is in Nika: NEXUS_WriteNx1DCanSASdata
https://github.com/usnistgov/PyHyperScattering/blob/main/src/PyHyperScattering/Nexus.py
or 
https://github.com/usnistgov/PyHyperScattering/blob/main/src/PyHyperScattering/FileIO.py


'''

# LoadNexus has code to load NXcanSAS data. Peter had also code to write it, I belive. 

import warnings
import xarray as xr
import numpy as np
import math
import pathlib
import h5py
import pathlib

'''
def saveNexus(self,fileName,compression=5):
    data = self._obj
    timestamp = datetime.datetime.now()
    # figure out if xr is a raw or integrated array
    
    axes = list(data.indexes.keys())
    array_to_save = data.variable.to_numpy()
    dims_of_array_to_save = data.variable.dims

    dim_to_index = {}
    index_to_dim = {}
    
    for n,dim in enumerate(dims_of_array_to_save):
        dim_to_index[dim] = n
        index_to_dim[n] = dim
    
    if 'pix_x' in axes:
        self.pyhyper_type='raw'
        nonspatial_coords = axes
        nonspatial_coords.remove('pix_x')
        nonspatial_coords.remove('pix_y')
    elif 'qx' in axes:
        self.pyhyper_type='qxy'
        nonspatial_coords = axes
        nonspatial_coords.remove('qx')
        nonspatial_coords.remove('qy')
    elif 'chi' in axes:
        self.pyhyper_type='red2d'
        nonspatial_coords = axes
        nonspatial_coords.remove('chi')
        nonspatial_coords.remove('q')
    elif 'q' in axes:
        self.pyhyper_type='red1d'
        nonspatial_coords = axes
        nonspatial_coords.remove('q')
    else:
        raise Exception(f'Invalid PyHyper_type {self.pyhyper_type}.  Cannot write Nexus.')
    raw_axes = list(data.indexes.keys())
        
        # create the HDF5 NeXus file
    with h5py.File(fileName, "w") as f:
        # point to the default data to be plotted
        f.attrs['default']          = 'entry'
        # give the HDF5 root some more attributes
        f.attrs['file_name']        = fileName
        f.attrs['file_time']        = str(timestamp)
        #f.attrs['instrument']       = 'CyRSoXS v'
        f.attrs['creator']          = 'PyHyperScattering NeXus writer'
        f.attrs['creator_version']  = phs_version
        f.attrs['NeXus_version']    = '4.3.0'
        f.attrs['HDF5_version']     = six.u(h5py.version.hdf5_version)
        f.attrs['h5py_version']     = six.u(h5py.version.version)

        # create the NXentry group
        nxentry = f.create_group('entry')
        nxentry.attrs['NX_class'] = 'NXentry'
        nxentry.attrs['canSAS_class'] = 'SASentry'
        nxentry.attrs['default'] = 'data'
        #nxentry.create_dataset('title', data='SIMULATION NAME GOES HERE') # do we have a sample name field?

        #setup general file stuff
        nxdata = nxentry.create_group('sasdata')
        nxdata.attrs['NX_class'] = 'NXdata'
        nxdata.attrs['canSAS_class'] = 'SASdata'
        nxdata.attrs['canSAS_version'] = '0.1' #required for Nika to read the file.
        nxdata.attrs['signal'] = 'I'      # Y axis of default plot
            

        
        '''if self.pyhyper_type == 'raw':
            nxdata.attrs['I_axes'] = 'pix_x,pix_y'         # X axis of default plot
            nxdata.attrs['Q_indices'] = f'[{dim_to_index["pix_x"]},{dim_to_index["pix_y"]}]'               
        else:
            if self.pyhyper_type == 'qxy':
                nxdata.attrs['I_axes'] = 'Qx,Qy'         # X axis of default plot
                nxdata.attrs['Q_indices'] = f'[{dim_to_index["Qx"]},{dim_to_index["Qy"]}]'  
            elif self.pyhyper_type == 'red2d':
                nxdata.attrs['I_axes'] = 'q,chi'         # X axis of default plot
                nxdata.attrs['Q_indices'] = f'[{dim_to_index["q"]},{dim_to_index["chi"]}]' 
            elif self.pyhyper_type == 'red1d':
                nxdata.attrs['I_axes'] = 'q'         # X axis of default plot
                nxdata.attrs['Q_indices'] = f'[{dim_to_index["q"]}]' 
            else:
                raise Exception(f'Invalid PyHyper_type {self.pyhyper_type}.  Cannot write Nexus.')
        '''
        
        ds = nxdata.create_dataset('I', data=array_to_save,compression=compression)
        ds.attrs['units'] = 'arbitrary'
        ds.attrs['long_name'] = 'Intensity (arbitrary units)'    # suggested X axis plot label
        # the following are to enable compatibility with Nika canSAS loading
        # ds.attrs['signal'] = 1
        #ds.attrs['axes'] = 'Qx,Qy'
        I_axes = '['
        Q_indices = '['
        for axis in raw_axes:
            I_axes += f'{axis},'
            Q_indices += f'{dim_to_index[axis]},'
            if type(data.indexes[axis]) == pandas.core.indexes.multi.MultiIndex:
                idx = data.indexes[axis]
                I_axes = I_axes[:-1]+'('
                lvls = idx.levels
                multiindex_arrays = defaultdict(list)
                for row in idx:
                    for n,level in enumerate(lvls):
                        multiindex_arrays[level.name].append(row[n])
                for level in lvls:
                    ds = nxdata.create_dataset(level.name, data=multiindex_arrays[level.name])
                    I_axes += f'{level.name};'
                    ds.attrs['PyHyper_origin'] = axis
                I_axes = I_axes[:-1]+'),'
            else:
                ds = nxdata.create_dataset(data.indexes[axis].name, data=data.indexes[axis].values)
                if 'q' in axis:
                    ds.attrs['units'] = '1/angstrom'
                elif 'chi' in axis:
                    ds.attrs['units'] = 'degree'
                #ds.attrs['long_name'] = 'Qx (A^-1)'    # suggested Y axis plot label
        I_axes = I_axes[:-1]+']'
        Q_indices = Q_indices[:-1]+']'
        nxdata.attrs['I_axes'] = I_axes
        nxdata.attrs['Q_indices'] = Q_indices
                    
        residual_attrs = nxentry.create_group('attrs')
        residual_attrs = self._serialize_attrs(residual_attrs,data.attrs.items())
        '''for k,v in data.attrs.items():
            print(f'Serializing {k}...')
            print(f'Data: type {type(v)}, data {v}')
            if type(v)==datetime.datetime:
                ds = residual_attrs.create_dataset(k,data=v.strftime('%Y-%m-%dT%H:%M:%SZ'))
                ds.attrs['phs_encoding'] = 'strftime-%Y-%m-%dT%H:%M:%SZ'
            elif type(v)==dict:
                ds = residual_attrs.create_group(k)
                
            else:
                try:
                    residual_attrs.create_dataset(k, data=v)
                except TypeError:
                    ds = residual_attrs.create_dataset(k, data=json.dumps(v))
                    ds.attrs['phs_encoding'] = 'json'
        '''
    print("wrote file:", fileName)

'''
    
'''
def save(xr,fileName):
    
    # figure out if xr is a raw or integrated array
        # create the HDF5 NeXus file
    with h5py.File(fileName, "w") as f:
        # point to the default data to be plotted
        f.attrs['default']          = 'entry'
        # give the HDF5 root some more attributes
        f.attrs['file_name']        = fileName
        f.attrs['file_time']        = timestamp
        f.attrs['instrument']       = 'CyRSoXS v'
        f.attrs['creator']          = 'PyHyperScattering NeXus writer'
        f.attrs['NeXus_version']    = '4.3.0'
        f.attrs['HDF5_version']     = six.u(h5py.version.hdf5_version)
        f.attrs['h5py_version']     = six.u(h5py.version.version)

        # create the NXentry group
        nxentry = f.create_group('entry')
        nxentry.attrs['NX_class'] = 'NXentry'
        nxentry.attrs['canSAS_class'] = 'SASentry'
        nxentry.attrs['default'] = 'data'
        nxentry.create_dataset('title', data='SIMULATION NAME GOES HERE') #@TODO

        # figure out if one image or more
        if 'system' in xr.dimensions:
            # writing a stack of images
            pass
        else: 
            # writing a single image
            
            if single_image_energy is not None:
                imgpos = coords['energy'][single_image_energy].index()
            else:
                imgpos = 0
            nxdata = nxentry.create_group('sasdata_singleimg')
            nxdata.attrs['NX_class'] = 'NXdata'
            nxdata.attrs['canSAS_class'] = 'SASdata'
            #nxdata.attrs['canSAS_version'] = '0.1' #required for Nika to read the file.
            nxdata.attrs['signal'] = 'I'      # Y axis of default plot
            nxdata.attrs['I_axes'] = 'Qx,Qy'         # X axis of default plot
            nxdata.attrs['Q_indices'] = '[0,1]'   # use "mr" as the first dimension of I00

            # X axis data
            ds = nxdata.create_dataset('I', data=data['img'][imgpos])
            ds.attrs['units'] = 'arbitrary'
            ds.attrs['long_name'] = 'Intensity (arbitrary units)'    # suggested X axis plot label
            # the following are to enable compatibility with Nika canSAS loading
           # ds.attrs['signal'] = 1
            #ds.attrs['axes'] = 'Qx,Qy'

            # Y axis data
            ds = nxdata.create_dataset('Qx', data=data['qx'][0])
            ds.attrs['units'] = '1/angstrom'
            ds.attrs['long_name'] = 'Qx (A^-1)'    # suggested Y axis plot label

            ds = nxdata.create_dataset('Qy', data=data['qy'][0])
            ds.attrs['units'] = '1/angstrom'
            ds.attrs['long_name'] = 'Qy (A^-1)'    # suggested Y axis plot label


        if write_stack_qxy:
            # create the NXdata group for I(Qx,Qy)
            nxdata = nxentry.create_group('sasdata_energyseries')
            nxdata.attrs['NX_class'] = 'NXdata'
            nxdata.attrs['canSAS_class'] = 'SASdata'
            nxdata.attrs['signal'] = 'I'      # Y axis of default plot
            nxdata.attrs['I_axes'] = 'Qx,Qy,E'         # X axis of default plot
            nxdata.attrs['Q_indices'] = '[0,1]'   # use "mr" as the first dimension of I00

            # X axis data
            ds = nxdata.create_dataset('I', data=np.swapaxes(np.swapaxes(data['img'],0,1),1,2))
            ds.attrs['units'] = 'arbitrary'
            ds.attrs['long_name'] = 'Simulated Intensity (arbitrary units)'    # suggested X axis plot label

            # Y axis data
            ds = nxdata.create_dataset('Qx', data=data['qx'][0])
            ds.attrs['units'] = '1/angstrom'
            ds.attrs['long_name'] = 'Qx (A^-1)'    # suggested Y axis plot label

            ds = nxdata.create_dataset('Qy', data=data['qy'][0])
            ds.attrs['units'] = '1/angstrom'
            ds.attrs['long_name'] = 'Qy (A^-1)'    # suggested Y axis plot label


            ds = nxdata.create_dataset('E', data=coords['energy'])
            ds.attrs['units'] = 'eV'
            ds.attrs['long_name'] = 'Simulation Energy (eV)'    # suggested Y axis plot label
        if write_stack_qphi:
            # create the NXdata group for I(Q,phi)
            nxdata = nxentry.create_group('sasdata_unwrap')
            nxdata.attrs['NX_class'] = 'NXdata'
            nxdata.attrs['canSAS_class'] = 'SASdata'
            nxdata.attrs['signal'] = 'I'      # Y axis of default plot
            nxdata.attrs['I_axes'] = 'E,chi,Q'         # X axis of default plot
            nxdata.attrs['Q_indices'] = [2]   # use "mr" as the first dimension of I00

            # X axis data
            ds = nxdata.create_dataset('I', data=data['img'])
            ds.attrs['units'] = 'arbitrary'
            ds.attrs['long_name'] = 'Simulated Intensity (arbitrary units)'    # suggested X axis plot label

            # Y axis data
            ds = nxdata.create_dataset('Q', data=coords['q'])
            ds.attrs['units'] = '1/angstrom'
            ds.attrs['long_name'] = 'Q (A^-1)'    # suggested Y axis plot label

            ds = nxdata.create_dataset('chi', data=coords['chi'])
            ds.attrs['units'] = 'degree'
            ds.attrs['long_name'] = 'azimuthal angle chi (deg)'    # suggested Y axis plot label


            ds = nxdata.create_dataset('E', data=coords['energy'])
            ds.attrs['units'] = 'eV'
            ds.attrs['long_name'] = 'Simulation Energy (eV)'    # suggested Y axis plot label
        if write_oned_traces:
            # create the NXdata group for I(Q,E) at two fixed orientations horizontal and vertical
            nxdata = nxentry.create_group('sasdata_horizontal')
            nxdata.attrs['NX_class'] = 'NXdata'
            nxdata.attrs['canSAS_class'] = 'SASdata'
            nxdata.attrs['signal'] = 'I'      # Y axis of default plot
            nxdata.attrs['I_axes'] = 'E,Q'         # X axis of default plot
            nxdata.attrs['Q_indices'] = [1]   # use "mr" as the first dimension of I00

            # X axis data
            ds = nxdata.create_dataset('I', data=data['Ihoriz'])
            ds.attrs['units'] = 'arbitrary'
            ds.attrs['long_name'] = 'Simulated Intensity (arbitrary units)'    # suggested X axis plot label

            # Y axis data
            ds = nxdata.create_dataset('Q', data=coords['q'])
            ds.attrs['units'] = '1/angstrom'
            ds.attrs['long_name'] = 'Q (A^-1)'    # suggested Y axis plot label

            ds = nxdata.create_dataset('E', data=coords['energy'])
            ds.attrs['units'] = 'eV'
            ds.attrs['long_name'] = 'Simulated Photon Energy (eV)'    # suggested Y axis plot label

             # create the NXdata group for I(Q,E) at two fixed orientations horizontal and vertical
            nxdata = nxentry.create_group('sasdata_vertical')
            nxdata.attrs['NX_class'] = 'NXdata'
            nxdata.attrs['canSAS_class'] = 'SASdata'
            nxdata.attrs['signal'] = 'I'      # Y axis of default plot
            nxdata.attrs['I_axes'] = 'E,Q'         # X axis of default plot
            nxdata.attrs['Q_indices'] = [1]   # use "mr" as the first dimension of I00

            # X axis data
            ds = nxdata.create_dataset('I', data=data['Ivert'])
            ds.attrs['units'] = 'arbitrary'
            ds.attrs['long_name'] = 'Simulated Intensity (arbitrary units)'    # suggested X axis plot label

            # Y axis data
            ds = nxdata.create_dataset('Q', data=coords['q'])
            ds.attrs['units'] = '1/angstrom'
            ds.attrs['long_name'] = 'Q (A^-1)'    # suggested Y axis plot label

            ds = nxdata.create_dataset('E', data=coords['energy'])
            ds.attrs['units'] = 'eV'
            ds.attrs['long_name'] = 'Simulated Photon Energy (eV)'    # suggested Y axis plot label

       
        # create the NXinstrument metadata group
        nxinstr = nxentry.create_group('instrument')
        nxinstr.attrs['NX_class'] = 'NXinstrument'
        nxinstr.attrs['canSAS_class'] = 'SASinstrument'

        nxprocess = nxinstr.create_group('simulation_engine')
        nxprocess.attrs['NX_class'] = 'NXprocess'
        nxprocess.attrs['canSAS_class'] = 'SASprocess'
        nxprocess.attrs['name'] = 'CyRSoXS Simulation Engine'
        nxprocess.attrs['date'] = timestamp # @TODO: get timestamp from simulation run and embed here.
        nxprocess.attrs['description'] = 'Simulation of RSoXS pattern from optical constants del/beta and morphology'

        sim_notes = nxprocess.create_group('NOTE')
        sim_notes.attrs['NX_class'] = 'NXnote'

        sim_notes.attrs['description'] = 'Simulation Engine Input Parameters/Run Data'
        sim_notes.attrs['author'] = 'CyRSoXS PostProcessor'
        sim_notes.attrs['data'] = 'Run metadata goes here' #@TODO

        for key in config:
            if 'Energy' in key:
                units = 'eV'
            elif 'Angle' in key:
                units = 'degree'
            elif 'PhysSize' in key:
                units = 'nm' #@TODO: is this correct?
            else:
                units = ''

            metads = sim_notes.create_dataset(key,(config[key],),dtype='f')
            metads.attrs['units'] = units
        nxsample = nxentry.create_group('sample')
        nxsample.attrs['NX_class'] = 'NXsample'
        nxsample.attrs['canSAS_class'] = 'SASsample'
        
        nxsample.attrs['name'] = 'SAMPLE NAME GOES HERE'
        nxsample.attrs['description'] = 'SAMPLE DESCRIPTION GOES HERE'
        nxsample.attrs['type'] = 'simulated data'

        #comp.create_dataset
        
    print("wrote file:", fileName)
    if 'pix_x' in xr.dimensions:
        pass
    elif 'q' in xr.dimensions:
        pass
    else:
        raise NotImplementedError(f'I do not support xarrays with dimensions of {xr.dimensions}')
    
def load(path):
    if type(path) is str:
        raise NotImplementedError

'''
            