# this file will convert SAXS and WAXS area detector data from the HDF5 format to the 1Ddata
# we will use pyFai to do the conversion

import pyFAI
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from pyFAI.integrator.azimuthal import AzimuthalIntegrator
import h5py
from supportFunctions import read_group_to_dict
import pprint as pp

# Load your 2D image data
# For this example, let's assume you have a NumPy array `image_data` representing your detector image
#image_data = np.random.random((1000, 1000))  # Replace with your actual image data

## main code here
#def ImportFlyscan(path, filename):
# Open the HDF5 file and read its content, parse content in numpy arrays and dictionaries
#with h5py.File(path+"/"+filename, 'r') as file:
with h5py.File("WAXS.hdf", 'r') as file:
    #read various data sets
    dataset = file['/entry/data/data'] 
    my2DData = np.array(dataset)
    #metadata
    instrument_group = file['/entry/instrument']
    instrument_dict = read_group_to_dict(instrument_group)
    #metadata
    metadata_group = file['/entry/Metadata']
    metadata_dict = read_group_to_dict(metadata_group)    
    #pp.pprint(metadata_dict)
    # wavelength, convert to m
    wavelength = instrument_dict["monochromator"]["wavelength"]* 1e-10 
     # pixel_size, convert to m
    pixel_size1 = instrument_dict["detector"]["x_pixel_size"]* 1e-3 
    pixel_size2 = instrument_dict["detector"]["y_pixel_size"]* 1e-3 
    # detector_distance, convert to m
    detector_distance = instrument_dict["detector"]["distance"]* 1e-3 
    # poni1, point of intercept. 
    poni1 = instrument_dict["detector"]["beam_center_y"] 
    poni2 = instrument_dict["detector"]["beam_center_x"]
    poni1 = poni1 * pixel_size1 
    poni2 = poni2 * pixel_size2 
    if "pin_ccd_tilt_x" in metadata_dict:
        rot1 = metadata_dict["pin_ccd_tilt_x"]*np.pi/180
        rot2 = metadata_dict["pin_ccd_tilt_y"]*np.pi/180
    else:
        rot1 = metadata_dict["waxs_ccd_tilt_x"]*np.pi/180
        rot2 = metadata_dict["waxs_ccd_tilt_y"]*np.pi/180     
        
plt.imshow(my2DData, cmap='viridis',norm=LogNorm())  # You can choose different colormaps
plt.colorbar()  # Optional: Add a colorbar to show the scale
plt.title('2D Array Visualization')
plt.show()

# Define your detector geometry
# You need to specify parameters like the detector distance, pixel size, and wavelength
#detector_distance = 0.1  # in meters
#pixel_size = 0.0001  # in meters
#wavelength = 1.54e-10  # in meters (for example, Cu K-alpha)

# Create an AzimuthalIntegrator object
ai = AzimuthalIntegrator(dist=detector_distance, poni1=poni1, poni2=poni2, rot1=rot1, rot2=rot2,
                         pixel1=pixel_size1, pixel2=pixel_size2, 
                         wavelength=wavelength)

# Perform azimuthal integration
# You can specify the number of bins for the integration
npt = 500  # Number of bins
q, intensity = ai.integrate1d(my2DData, npt, correctSolidAngle=False, unit="q_A^-1")

# Plot the integrated intensity
plt.plot(q, intensity)
plt.xlabel("q (1/Å)")
plt.ylabel("Intensity")
plt.title("Azimuthal Integration")
plt.show()