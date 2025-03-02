import h5py
import xarray as xr     #this needs xarray to run. https://docs.xarray.dev/en/stable/index.html
import json
import warnings
import datetime
import pandas           #needs pandas for multiindex support



def loadNexus(filename):
    with h5py.File(filename, "r") as f:    
        ds = xr.DataArray(f['entry']['sasdata']['I'],
                  dims=_parse_Iaxes(f['entry']['sasdata'].attrs['I_axes']),
                 coords = _make_coords(f))

        
        loaded_attrs = _unserialize_attrs(f['entry']['attrs'],{})
        
        '''        for entry in f['entry']['attrs']:
            #print(f'Processing attribute entry {entry}')
            try:
                encoding = f['entry']['attrs'][entry].attrs['phs_encoding']
                #print(f'Found data with a labeled encoding: {encoding}')
                if encoding == 'json':
                    loaded_attrs[entry] = json.loads(f['entry']['attrs'][entry][()].decode())
                elif encoding == 'dict-expanded':
                    loaded_attrs[entry] = self.load_attrs(entry)
                elif 'strftime' in encoding:
                    loaded_attrs[entry] = datetime.datetime.strptime(str(f['entry']['attrs'][entry][()].decode()),
                                                           encoding.replace('strftime-',''))
                else:
                    warnings.warn(f'Unknown phs_encoding {encoding} while loading {entry}.  Possible version mismatch.  Loading as string.',stacklevel=2)
                    loaded_attrs[entry] = f['entry']['attrs'][entry][()]
            except KeyError:
                loaded_attrs[entry] = f['entry']['attrs'][entry][()]'''
        #print(f'Loaded: {loaded_attrs}')
        ds.attrs.update(loaded_attrs)

    return ds


def _unserialize_attrs(hdf,attrdict):
    for entry in hdf:
        #print(f'Processing attribute entry {entry}')
        try:
            encoding = hdf[entry].attrs['phs_encoding']
            #print(f'Found data with a labeled encoding: {encoding}')
            if encoding == 'json':
                attrdict[entry] = json.loads(hdf[entry][()].decode())
            elif encoding == 'dict-expanded':
                attrdict[entry] = _unserialize_attrs(hdf[entry],{})
            elif 'strftime' in encoding:
                attrdict[entry] = datetime.datetime.strptime(str(hdf[entry][()].decode()),
                                                       encoding.replace('strftime-',''))
            else:
                warnings.warn(f'Unknown phs_encoding {encoding} while loading {entry}.  Possible version mismatch.  Loading as string.',stacklevel=2)
                attrdict[entry] = hdf[entry][()]
        except KeyError:
            attrdict[entry] = hdf[entry][()]        
    return attrdict
def _parse_Iaxes(axes,suppress_multiindex=True):
    axes = axes.replace('[','').replace(']','')
    axes_parts = axes.split(',')
    axes = []
    if suppress_multiindex:
        for part in axes_parts:
            if '(' in part:
                #print(f'multiindex: {part}')
                part = part.split('(')[0]
                #print(f'set part to {part}')
            axes.append(part)
    else:
        axes = axes_parts
    return axes

def _parse_multiindex_Iaxes(axis):
    axis = axis.replace(')','')
    axis = axis.split('(')[1]
    return axis.split(';')

def _make_coords(f):
    axes = _parse_Iaxes(f['entry']['sasdata'].attrs['I_axes'],suppress_multiindex=True)
    axes_raw = _parse_Iaxes(f['entry']['sasdata'].attrs['I_axes'],suppress_multiindex=False)

    coords = {}
    for n,axis in enumerate(axes_raw):
        if '(' in axis:
            levels = _parse_multiindex_Iaxes(axis)
            vals = []
            names = []
            for level in levels:
                names.append(level)
                vals.append(f['entry']['sasdata'][level])
            #print(names)
            #print(vals)
            coords[axes[n]] = pandas.MultiIndex.from_arrays(vals,names=names)
        else:
            coords[axes[n]] = f['entry']['sasdata'][axis]

    return coords

