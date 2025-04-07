# this file will convert SAXS and WAXS area detector data from the HDF5 format to the 1Ddata
# we will use pyFai to do the conversion

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from pyFAI.integrator.azimuthal import AzimuthalIntegrator
import h5py
from supportFunctions import read_group_to_dict, filter_nested_dict
import pprint as pp
import os
import tifffile as tiff
import logging
from convertNikaTopyFAI import convert_Nika_to_Fit2D

# Load your 2D image data
# For this example, let's assume you have a NumPy array `image_data` representing your detector image
#image_data = np.random.random((1000, 1000))  # Replace with your actual image data

## main code here
def ImportAndReduceAD(path, filename):
    # Open the HDF5 file and read its content, parse content in numpy arrays and dictionaries
    #with h5py.File(path+"/"+filename, 'r') as file:
    with h5py.File(path+'/'+filename, 'r') as file:
        #read various data sets
        logging.info(f"Read file :{filename}")
        dataset = file['/entry/data/data'] 
        my2DData = np.array(dataset)
        #metadata
        instrument_group = file['/entry/instrument']
        instrument_dict = read_group_to_dict(instrument_group)
        #metadata
        keys_to_keep = ['I000_cts', 'I00_cts', 'I00_gain', 'I0_cts', 'I0_gated',
                        'I0_gain', 'I_scaling', 'Pin_TrI0', 'Pin_TrI0gain', 'Pin_TrI0gain','Pin_TrPD','Pin_TrPDgain',
                        'PresetTime', 'monoE', 'pin_ccd_center_x_pixel','pin_ccd_center_y_pixel',
                        'pin_ccd_tilt_x', 'pin_ccd_tilt_y', 'wavelength', 'waxs_ccd_center_x', 'waxs_ccd_center_y',
                        'waxs_ccd_tilt_x', 'waxs_ccd_tilt_y', 'waxs_ccd_center_x_pixel', 'waxs_ccd_center_y_pixel',
                        'scaler_freq'                     
                    ]        
        metadata_group = file['/entry/Metadata']
        metadata_dict = read_group_to_dict(metadata_group)
        metadata_dict = filter_nested_dict(metadata_dict, keys_to_keep)

        #logging.info(f"Metadata : {metadata_group}")
        #pp.pprint(metadata_dict)
        # wavelength, keep in A for Fit2D
        wavelength = instrument_dict["monochromator"]["wavelength"]
        # pixel_size, convert to um
        pixel_size1 = instrument_dict["detector"]["x_pixel_size"]*1000 #in mm in NIka, in micron in Fit2D
        pixel_size2 = instrument_dict["detector"]["y_pixel_size"]*1000 #in mm in NIka, in micron in Fit2D
        # detector_distance, keep in mm for Fit2D
        detector_distance = instrument_dict["detector"]["distance"] #in mm in Nika, in mm in Fit2D
        #logging.info(f"Read metadata")
        # convert to Fit2D syntax and then use pyFAI function to convert to poni
        #TODO
        if "pin_ccd_tilt_x" in metadata_dict:
            usingWAXS=0
            BCY = instrument_dict["detector"]["beam_center_y"]      # they are swapped later in conversion routine
            BCX= instrument_dict["detector"]["beam_center_x"]       # they are swapped later in conversion routine
            HorTilt = metadata_dict["pin_ccd_tilt_x"]              #keep in degrees for Fit2D
            VertTilt = metadata_dict["pin_ccd_tilt_y"]              #keep in degrees for Fit2D
        else:
            usingWAXS=1
            BCY = instrument_dict["detector"]["beam_center_y"]      # they are swapped later in conversion routine
            BCX = instrument_dict["detector"]["beam_center_x"]      # they are swapped later in conversion routine
            HorTilt = metadata_dict["waxs_ccd_tilt_x"]             #keep in degrees for Fit2D
            VertTilt = metadata_dict["waxs_ccd_tilt_y"]             #keep in degrees for Fit2D    

        #logging.info(f"Finished reading metadata")
    poni = convert_Nika_to_Fit2D(detector_distance, pixel_size1, BCX, BCY, HorTilt, VertTilt, wavelength)
    # poni is geomtry file for pyFAI, created by converting first to Fit2D and then calling pyFAI conversion function.

    #create mask here. Duplicate the my2DData and set all values above 1e7 to NaN
    #this is for waxs ONLY, DIFFERENT FOR saxs
    if usingWAXS:
        mask = np.copy(my2DData)
        mask = 0*mask   # set all values to zero
        mask[my2DData > 1e7] = 1
    else:
        mask = np.copy(my2DData)
        mask = 0*mask   # set all values to zero
        mask[my2DData < 0] = 1
        # Set the first 4 rows to 1
        mask[:, :4] = 1
        # Set rows 192 to 195 to 1
        mask[:, 242:245] = 1

    #logging.info(f"Finished creating mask")

    # Create an AzimuthalIntegrator object
    # old method ai = AzimuthalIntegrator(dist=detector_distance, poni1=poni1, poni2=poni2, rot1=rot1, rot2=rot2,
    #                        pixel1=pixel_size1, pixel2=pixel_size2, 
    #                        wavelength=wavelength)
    ai = AzimuthalIntegrator(poni)

    # Perform azimuthal integration
    #   You can specify the number of bins for the integration
    #   set npt to larger of dimmension of my2DData  `
    npt = max(my2DData.shape)
        #npt = 1000  # Number of bins, if should be lower
    q, intensity = ai.integrate1d(my2DData, npt, mask=mask, correctSolidAngle=True, unit="q_A^-1")
    #logging.info(f"Finished 2d to 1D conversion")
    result = {"Intensity":np.ravel(intensity), "Q_array":np.ravel(q)}
    #pp.pprint(result)
    return result



# intensity=np.log(intensity)
# # Plot the integrated intensity
# plt.plot(q, intensity)
# plt.xlabel("q (1/Å)")
# plt.ylabel("Intensity")
# plt.title("Azimuthal Integration")
# plt.show()


def reduceADToQR(path, filename):
        tempFilename= os.path.splitext(filename)[0]
        tempSample = {"RawData":{"Filename": tempFilename}}
        # label = data_dict["RawData"]["Filename"]
        # Q_array = data_dict["ReducedData"]["Q_array"]
        # Intensity = data_dict["ReducedData"]["UPD"]
        tempSample["ReducedData"]=ImportAndReduceAD(path, filename)
        #pp.pprint(tempSample)
        #pp.pprint(tempSample["RawData"]["Filename"])
        return tempSample


def PlotResults(data_dict):
    # Find Peak center and create Q vector.
    Q_array = data_dict["ReducedData"]["Q_array"]
    Intensity = data_dict["ReducedData"]["Intensity"]
    
    # Plot ydata against xdata
    plt.figure(figsize=(6, 12))
    plt.plot(Q_array, Intensity, linestyle='-')  # You can customize the marker and linestyle
    plt.title('Plot of Intensity vs. Q')
    plt.xlabel('log(Q) [1/A]')
    plt.ylabel('Intensity')
    plt.xscale('linear')
    plt.yscale('linear')
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    Sample = dict()
    Sample=reduceADToQR("/home/parallels/Github/Matilda","SAXS.hdf")
    #Sample["ReducedData"]=test("/home/parallels/Github/Matilda","LaB6_45deg.tif")
    #pp.pprint(Sample)
    PlotResults(Sample)
