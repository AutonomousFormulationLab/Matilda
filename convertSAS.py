# this file will convert SAXS and WAXS area detector data from the HDF5 format to the 1Ddata
# we will use pyFai to do the conversion

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from pyFAI.integrator.azimuthal import AzimuthalIntegrator
import h5py
from supportFunctions import read_group_to_dict
import pprint as pp
import os
import tifffile as tiff
import logging


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
        metadata_group = file['/entry/Metadata']
        metadata_dict = read_group_to_dict(metadata_group)    
        logging.info(f"Metadata : {metadata_group}")
        #pp.pprint(metadata_dict)
        # wavelength, convert to m
        wavelength = instrument_dict["monochromator"]["wavelength"]* 1e-10 
        # pixel_size, convert to m
        pixel_size1 = instrument_dict["detector"]["x_pixel_size"]* 1e-3 
        pixel_size2 = instrument_dict["detector"]["y_pixel_size"]* 1e-3 
        # detector_distance, convert to m
        detector_distance = instrument_dict["detector"]["distance"]* 1e-3 
        # poni1, point of intercept. 
        logging.info(f"Read metadata")
        if "pin_ccd_tilt_x" in metadata_dict:
            usingWAXS=0
            BCY = instrument_dict["detector"]["beam_center_y"] # based on Peter's code this shoudl be opposite
            BCX = instrument_dict["detector"]["beam_center_x"] # poni2 shoudl be x and poni1 should be y
            BCY = BCY * pixel_size2 
            BCX = BCX * pixel_size1 
            rotX = metadata_dict["pin_ccd_tilt_x"]*np.pi/180
            rotY = metadata_dict["pin_ccd_tilt_y"]*np.pi/180
        else:
            usingWAXS=1
            BCY = instrument_dict["detector"]["beam_center_y"] # based on Peter's code this shoudl be opposite
            BCX = instrument_dict["detector"]["beam_center_x"] # poni2 shoudl be x and poni1 should be y
            BCY = BCY * pixel_size2 
            BCX = BCX * pixel_size1 
            rotX = metadata_dict["waxs_ccd_tilt_x"]*np.pi/180
            rotY = metadata_dict["waxs_ccd_tilt_y"]*np.pi/180     

        logging.info(f"Finished reading metadata")

    # now we need to do correction on geometry. 
    # Nika rotates the detector around the beam center while pyFAI around sample
    # this corrections testing is based on pyFAI-calib2 testing for 45 deg tilt tiff file from Bob GSAS-II
    # first correct distacne:
    #detector_distance = detector_distance*abs(np.cos(rotX)*np.cos(rotY))
    # now calculate poni1 and poni2
    # testing shows, that poni1 is related to BCY and rotY
    # but rot1 is = rotX
    # poni1 = BCX + detector_distance*np.tan(rotX)
    # poni2 = BCY + detector_distance*np.tan(rotY)
    # rot1=rotX
    # rot2=rotY   
    #override for for now below:
    #TODO: need to figrue this out more on new test dat aset in hdf5 file with 
    #proper calibration in parameters...        
    rot1=0
    rot2=0
    poni1=BCX
    poni2 = BCY
    # print(f"wavelength: {wavelength}")
    # print(f"pixel_size1: {pixel_size1}")
    # print(f"pixel_size2: {pixel_size2}")
    # print(f"detector_distance: {detector_distance}")
    # print(f"poni1: {poni1}")
    # print(f"poni2: {poni2}")
    # print(f"rot1: {rot1}")
    # print(f"rot2: {rot2}")    
    # print(f"rotX: {rotX}")
    # print(f"rotY: {rotY}")
                
    # plt.imshow(my2DData, cmap='viridis',norm=LogNorm())  # You can choose different colormaps
    # plt.colorbar()  # Optional: Add a colorbar to show the scale
    # plt.title('2D Array Visualization')
    # plt.show()

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

    logging.info(f"Finished creating mask")

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
    #set npt to larger of dimmension of my2DData  `
    npt = max(my2DData.shape)
    #npt = 1000  # Number of bins
    q, intensity = ai.integrate1d(my2DData, npt, mask=mask, correctSolidAngle=True, unit="q_A^-1")
    logging.info(f"Finished 2d to 1D conversion")

    result = {"Intensity":np.ravel(intensity), "Q_array":np.ravel(q)}
    #pp.pprint(result)
    return result


## test fro tilts goes here
def test(path, filename):
    # read data from tiff file and read the data 
    my2DData = tiff.imread(filename)
    # plt.imshow(my2DData, cmap='viridis',norm=LogNorm())  # You can choose different colormaps
    # plt.colorbar()  # Optional: Add a colorbar to show the scale
    # plt.title('2D Array Visualization')
    # plt.show()    # wavelength, convert to m
    wavelength = 0.10798* 1e-10 
    # pixel_size, convert to m
    pixel_size1 = 0.1* 1e-3 # x in Nika, in m
    pixel_size2 = 0.1* 1e-3 # y in Nika, in m
    # detector_distance, convert to m
    detector_distance = 1004.91* 1e-3 
    # poni1, point of intercept in y in Nika, in m
    # poni2, point of intercept in x in Nika, in m
    # read Nika BCX and BCY and convert to m
    BCX = 886.7     # x in Nika
    BCY = 1048.21   # y in Nika, width in Python
    BCY = BCY * pixel_size2 #in m now
    BCX = BCX * pixel_size1 #in m now
    # read Nika rot1 and rot2 and convert to radians
    rotX = 44.7*np.pi/180   # x direction in Nika
    rotY = 0.02*np.pi/180   # y direction in Nika
    # now corrections based on pyfain-geom2 testing for 45 deg tilt
    # first correct distacne:
    detector_distance = detector_distance*abs(np.cos(rotX)*np.cos(rotY))
    # now calculate poni1 and poni2
    # confusingly, the poni1 is realted to BCY and poni2 is related to BCX
    # is this issue with reading tiff vs hdf5 image orientations? 
    poni2 = BCX + detector_distance*np.tan(rotX)
    poni1 = BCY + detector_distance*np.tan(rotY)
    rot1=rotX
    rot2=rotY              
    #print(f"wavelength: {wavelength}")
    #print(f"pixel_size1: {pixel_size1}")
    #print(f"pixel_size2: {pixel_size2}")
    # print(f"detector_distance: {detector_distance}")
    # print(f"poni1: {poni1}")
    # print(f"poni2: {poni2}")
    # print(f"rot1: {rot1}")
    # print(f"rot2: {rot2}")
                
    # plt.imshow(my2DData, cmap='viridis',norm=LogNorm())  # You can choose different colormaps
    # plt.colorbar()  # Optional: Add a colorbar to show the scale
    # plt.title('2D Array Visualization')
    # plt.show()

    #create mask here. Duplicate the my2DData and set all values above 1e7 to NaN
    #this is for waxs ONLY, DIFFERENT FOR saxs
    mask = np.copy(my2DData)
    mask = 0*mask   # set all values to zero
 

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
    #set npt to larger of dimmension of my2DData  `
    npt = max(my2DData.shape)
    #npt = 1000  # Number of bins
    q, intensity = ai.integrate1d(my2DData, npt, mask=mask, correctSolidAngle=True, unit="q_A^-1")
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