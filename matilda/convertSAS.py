'''
    Will convert SAXS and WAXS area detector data from the HDF5 format to the 1Ddata
    we will use pyFai to do the conversion
    we will convert Nika parameters to Fit2D format and then use pyFAI to convert to poni format
    Only some metadata are kept to keep all more reasnobale on size

'''

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
from hdf5code import save_dict_to_hdf5, load_dict_from_hdf5


## main code here
def ImportAndReduceAD(path, filename, deleteExisting=False):
    # Open the HDF5 file and read its content, parse content in numpy arrays and dictionaries
    location = 'entry/displayData/'
    with h5py.File(path+'/'+filename, 'r+') as hdf_file:
        # Check if the group 'displayData' exists
        if deleteExisting:
            # Delete the group
            del hdf_file[location]
            print("Deleted existing group 'entry/displayData'.")

        if location in hdf_file:
            # exists, so lets reuse the data from the file
            Sample = dict()
            Sample = load_dict_from_hdf5(hdf_file, location)
            print("Used existing data")
            q = Sample["ReducedData"]["Q_array"]
            intensity = Sample["ReducedData"]["Intensity"]
            result = {"Intensity":np.ravel(intensity), "Q_array":np.ravel(q)}  
            return result
        
        
        else:
            Sample = dict()
            #read various data sets
            logging.info(f"Read file :{filename}")
            dataset = hdf_file['/entry/data/data'] 
            my2DData = np.array(dataset)
            #metadata
            instrument_group = hdf_file['/entry/instrument']
            instrument_dict = read_group_to_dict(instrument_group)
            #metadata
            keys_to_keep = ['I000_cts', 'I00_cts', 'I00_gain', 'I0_cts', 'I0_gated',
                            'I0_gain', 'I_scaling', 'Pin_TrI0', 'Pin_TrI0gain', 'Pin_TrI0gain','Pin_TrPD','Pin_TrPDgain',
                            'PresetTime', 'monoE', 'pin_ccd_center_x_pixel','pin_ccd_center_y_pixel',
                            'pin_ccd_tilt_x', 'pin_ccd_tilt_y', 'wavelength', 'waxs_ccd_center_x', 'waxs_ccd_center_y',
                            'waxs_ccd_tilt_x', 'waxs_ccd_tilt_y', 'waxs_ccd_center_x_pixel', 'waxs_ccd_center_y_pixel',
                            'scaler_freq'                     
                        ]        
            metadata_group = hdf_file['/entry/Metadata']
            metadata_dict = read_group_to_dict(metadata_group)
            metadata_dict = filter_nested_dict(metadata_dict, keys_to_keep)
            # wavelength, keep in A for Fit2D
            wavelength = instrument_dict["monochromator"]["wavelength"]
            # pixel_size, keep in mm, converted in convert_Nika_to_Fit2D to micron for Fit2D and then to m for pyFAI... 
            pixel_size1 = instrument_dict["detector"]["x_pixel_size"] #in mm in NIka, will convert to micron for Fit2D later
            # assume pixels are square, therefore size2 is not needed. No idea how to fix this in Nika or pyFAI for that matter anyway. 
            #pixel_size2 = instrument_dict["detector"]["y_pixel_size"] #in mm in NIka, will convert to micron for Fit2D later
            # detector_distance, keep in mm for Fit2D
            detector_distance = instrument_dict["detector"]["distance"] #in mm in Nika, in mm in Fit2D
            #logging.info(f"Read metadata")
            if "pin_ccd_tilt_x" in metadata_dict:                       # this is SAXS
                usingWAXS=0
                BCX= instrument_dict["detector"]["beam_center_x"]       #  This will be swapped later in convert_Nika_to_Fit2D 
                BCY = instrument_dict["detector"]["beam_center_y"]      #  This will be swapped later in convert_Nika_to_Fit2D
                HorTilt = metadata_dict["pin_ccd_tilt_x"]               #   keep in degrees for Fit2D
                VertTilt = metadata_dict["pin_ccd_tilt_y"]              #   keep in degrees for Fit2D
            else:                                                       # and this is WAXS
                usingWAXS=1
                BCX = instrument_dict["detector"]["beam_center_x"]      #  This will be swapped later in convert_Nika_to_Fit2D
                BCY = instrument_dict["detector"]["beam_center_y"]      #  This will be swapped later in convert_Nika_to_Fit2D
                HorTilt = metadata_dict["waxs_ccd_tilt_x"]              #   keep in degrees for Fit2D
                VertTilt = metadata_dict["waxs_ccd_tilt_y"]             #   keep in degrees for Fit2D    

            #logging.info(f"Finished reading metadata")
        
        # poni is geometry file for pyFAI, created by converting first to Fit2D and then calling pyFAI conversion function.
        my_poni = convert_Nika_to_Fit2D(SSD=detector_distance, pix_size=pixel_size1, BCX=BCX, BCY=BCY, HorTilt=HorTilt, VertTilt=VertTilt, wavelength=wavelength)

        #create mask here. Duplicate the my2DData and set all values above 1e7 to NaN for WAXS or for SAXS mask all negative intensities
        # the differecne is due to Pilatus vs Eiger handing bad pixels differently. Dectris issue... 
        if usingWAXS:
            mask = np.copy(my2DData)
            mask = 0*mask   # set all values to zero
            mask[my2DData > 1e7] = 1
            mask[:, 511:516] = 1
            mask[:, 1026:1041] = 1
            mask[:, 1551:1556] = 1
        else:
            mask = np.copy(my2DData)
            mask = 0*mask   # set all values to zero
            mask[my2DData < 0] = 1
            # Set the first 4 rows to 1
            mask[:, :4] = 1
            # Set rows 192 to 195 to 1
            mask[:, 242:245] = 1

        #logging.info(f"Finished creating mask")
        
        #now define integrator... using pyFAI here. 
        # this does not work, they really do not have way to pass whole poni in? 
        #ai = AzimuthalIntegrator(poni=my_poni)
        # but this works fine
        ai = AzimuthalIntegrator(dist=my_poni.dist, poni1=my_poni.poni1, poni2=my_poni.poni2, rot1=my_poni.rot1, rot2=my_poni.rot2,
                            rot3=my_poni.rot3, pixel1=my_poni.detector.pixel1, pixel2=my_poni.detector.pixel2, 
                            wavelength=my_poni.wavelength)
        
        #   You can specify the number of bins for the integration
        #   set npt to larger of dimmension of my2DData
        if usingWAXS:
            npt = max(my2DData.shape)
        else:
            npt=200 
        #npt = 1000  # Number of bins, if should be lower
        # Perform azimuthal integration
        q, intensity = ai.integrate1d(my2DData, npt, mask=mask, correctSolidAngle=True, unit="q_A^-1")
        #logging.info(f"Finished 2d to 1D conversion")
        Sample["ReducedData"] = dict()
        Sample["ReducedData"]["Q_array"] = q
        Sample["ReducedData"]["Intensity"] = intensity
        save_dict_to_hdf5(Sample, location, hdf_file)
        print("Appended new data to 'entry/displayData'.")
        result = {"Intensity":np.ravel(intensity), "Q_array":np.ravel(q)}
        return result


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
    Sample=reduceADToQR("./TestData/TestTiltData","LaB6_tilt7v_0049.hdf")
    #Sample["ReducedData"]=test("/home/parallels/Github/Matilda/TestData","LaB6_45deg.tif")
    #pp.pprint(Sample)
    PlotResults(Sample)



## test for tilts using LaB6 45 deg tilted detector from GSAXS-II goes here
# to the best of my undestanding, the images loaded from tiff file are mirrored and the values here are just weird. 
# def test(path, filename):
#     # read data from tiff file and read the data 
#     # tiff files are actually loaded differently than HDF5 files. Looks like they are mirrored. 
#     my2DData = tiff.imread(path+'/'+filename)
#     wavelength = 0.10798 # in A
#     # pixel_size
#     pixel_size1 = 0.1 # x in Nika, in mm
#     #pixel_size2 = 0.1 # y in Nika, in mm
#     # detector_distance, in mm
#     detector_distance = 1004.91 # in Nika, in mm 
#     # Nika BCX and BCY in pixels
#     BCY = 886.7     # this is for hdf5 x in Nika
#     BCX = 1048.21   # this is for hdf5 y in Nika
#     # read Nika HorTilt and VertTilt 
#     VertTilt  = -44.7   # this is negative value for horizontal tilt in Nika
#     HorTilt = 0.02      # this is value for vertical tilt in Nika, not sure if this shoudl be negative. 
#     # poni is geometry file for pyFAI, created by converting first to Fit2D and then calling pyFAI conversion function.
#     my_poni = convert_Nika_to_Fit2D(detector_distance, pixel_size1, BCX, BCY, HorTilt, VertTilt, wavelength)
#     # setup integrator geometry
#     ai = AzimuthalIntegrator(dist=my_poni.dist, poni1=my_poni.poni1, poni2=my_poni.poni2, rot1=my_poni.rot1, rot2=my_poni.rot2,
#                        rot3=my_poni.rot3, pixel1=my_poni.detector.pixel1, pixel2=my_poni.detector.pixel2, 
#                        wavelength=my_poni.wavelength)
#     #create mask here. Duplicate the my2DData and set all values to be masked to NaN, not used here. 
#     mask = np.copy(my2DData)
#     mask = 0*mask           # set all values to zero
#     # Perform azimuthal integration
#     # You can specify the number of bins for the integration
#     #set npt to larger of dimmension of my2DData  `
#     npt = max(my2DData.shape)
#     q, intensity = ai.integrate1d(my2DData, npt, mask=mask, correctSolidAngle=True, unit="q_A^-1")
#     result = {"Intensity":np.ravel(intensity), "Q_array":np.ravel(q)}
#     return result
