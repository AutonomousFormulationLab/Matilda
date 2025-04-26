'''
Here we develop new code which then moves to proper package
'''
import h5py
import numpy as np
#from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pprint as pp
from supportFunctions import subtract_data #read_group_to_dict, filter_nested_dict, check_arrays_same_length
import os
from convertUSAXS import rebinData
from hdf5code import save_dict_to_hdf5, load_dict_from_hdf5
from convertUSAXS import importFlyscan, calculatePD_Fly, beamCenterCorrection, smooth_r_data



#develop calibrated flyscan and step scan data routines
# this will check if NXcanSAS data exist and if not, it will create properly calibrated NXcanSAS
# data. It needs Blank and Sample data.
def reduceFlyscan(path, filename, deleteExisting=False):
    # Open the HDF5 file in read/write mode
    #location = 'entry/displayData/'
    with h5py.File(path+'/'+filename, 'r+') as hdf_file:
            # Check if the group 'location' exists, if yes, bail out as this is all needed. 
            # if deleteExisting:
            #     # Delete the group
            #     del hdf_file[location]
            #     print("Deleted existing group 'entry/displayData'.")

            # if location in hdf_file:
            #     # exists, so lets reuse the data from the file
            #     Sample = dict()
            #     Sample = load_dict_from_hdf5(hdf_file, location)
            #     print("Used existing data")
            #     return Sample
            # else:
                Sample = dict()
                Sample["RawData"]=importFlyscan(path, filename)                 #import data
                Sample["reducedData"]= calculatePD_Fly(Sample)                  # Creates PD_Intesnity with corrected gains and background subtraction
                Sample["reducedData"].update(calculatePDError(Sample))          # Calculate UPD error, mostly the same as in Igor                
                Sample["reducedData"].update(beamCenterCorrection(Sample,useGauss=0)) #Beam center correction
                Sample["reducedData"].update(smooth_r_data(Sample["reducedData"]["PD_intensity"],     #smooth data data
                                                        Sample["reducedData"]["Q_array"], 
                                                        Sample["reducedData"]["PD_range"], 
                                                        Sample["reducedData"]["PD_error"], 
                                                        Sample["RawData"]["TimePerPoint"],
                                                        replaceNans=True))                 
                Sample["BlankData"]=getBlankFlyscan()
                Sample["reducedData"].update(normalizeBlank(Sample))          # Normalize sample by dividing by transmission for subtraction
                Sample["CalibratedData"]=(calibrateAndSubtractFlyscan(Sample))
                #pp.pprint(Sample)
                # TODO: calibration
                # TODO: fix rebinning for 3 input waves returning 4 waves with dQ
                Sample["CalibratedData"].update(rebinData(Sample, num_points=500, isSMRData=True))         #Rebin data
                # TODO: desmearing here
                # Create the group and dataset for the new data inside the hdf5 file for future use. 
                #save_dict_to_hdf5(Sample, location, hdf_file)
                #print("Appended new data to 'entry/displayData'.")
                return Sample

def normalizeBlank(Sample):
    # This is a simple normalization of the blank data to the sample data. 
    # It will be used for background subtraction.
    PeakIntensitySample = Sample["reducedData"]["Maximum"]
    PeakIntensityBlank = Sample["BlankData"]["Maximum"]
    PeakToPeakTransmission = PeakIntensitySample/PeakIntensityBlank


    PD_intensity = Sample["reducedData"]["PD_intensity"]
    PD_intensity = PD_intensity / PeakToPeakTransmission
    PD_error = Sample["reducedData"]["PD_error"]
    PD_error = PD_error / PeakToPeakTransmission
    result = {"PD_intensity":PD_intensity,
            "PD_error":PD_error,
            "PeakToPeakTransmission":PeakToPeakTransmission
            }
    return result
    

def calibrateAndSubtractFlyscan(Sample):
    # This is a step wehre we subtract and calibrate the sample and Blank. 
    PD_intensity = Sample["reducedData"]["PD_intensity"]
    BL_PD_intensity = Sample["BlankData"]["PD_intensity"]
    PD_error = Sample["reducedData"]["PD_error"]
    BL_PD_error = Sample["BlankData"]["PD_error"]
    Q_array = Sample["reducedData"]["Q_array"]
    BL_Q_array = Sample["BlankData"]["Q_array"]

    SMR_Qvec, SMR_Int, SMR_Error = subtract_data(Q_array, PD_intensity,PD_error, BL_Q_array, BL_PD_intensity, BL_PD_error)
    # TODO: trim, calibrate, 
    # find Qmin as the first point where we get above 3% of the background avleu and larger than instrument resolution
    IntRatio = PD_intensity / BL_PD_intensity
    # find point where the IntRatio is larger than 1.03
    FWHMSample = Sample["reducedData"]["FWHM"]
    FWHMBlank = Sample["BlankData"]["FWHM"]
    wavelength =  Sample["reducedData"]["wavelength"]
    PeakToPeakTransmission =  Sample["reducedData"]["PeakToPeakTransmission"]
    SaTransCounts = Sample['RawData']['metadata']['trans_pin_counts']
    SaTransGain = Sample['RawData']['metadata']['trans_pin_gain']
    SaI0Counts = Sample['RawData']['metadata']['trans_I0_counts']
    SaI0Gain = Sample['RawData']['metadata']['trans_I0_gain'] 
    BlTransCounts = Sample['BlankData']['BlTransCounts']
    BlTransGain = Sample['BlankData']['BlTransGain']
    BlI0Counts = Sample['BlankData']['BlI0Counts']
    BlI0Gain = Sample['BlankData']['BlI0Gain']
    MeasuredTransmission = ((SaTransCounts/SaTransGain)/(SaI0Counts/SaI0Gain))/((BlTransCounts /BlTransGain )/(BlI0Counts/BlI0Gain))
    MSAXSCorrection = MeasuredTransmission / PeakToPeakTransmission
    QminSample = 4*np.pi*np.sin(np.radians(FWHMSample)/2)/wavelength
    QminBlank = 4*np.pi*np.sin(np.radians(FWHMBlank)/2)/wavelength
    indexSample = np.searchsorted(Q_array, QminSample)
    indexBlank = np.searchsorted(Q_array, QminBlank)
    indexRatio =  np.searchsorted(IntRatio, 1.03)
    largest_value = max(indexSample, indexBlank, indexRatio)
    SMR_Qvec = SMR_Qvec[largest_value-1 : ]    
    SMR_Int = SMR_Int[largest_value-1 : ]    
    SMR_Error = SMR_Error[largest_value-1 : ]
    # now calibration... 
    SDD = Sample["RawData"]["metadata"]['detector_distance']
    UPDSize =  Sample["RawData"]["metadata"]['UPDsize']
    thickness = Sample["RawData"]["sample"]['thickness']
    BLPeakMax = Sample["BlankData"]["Maximum"]
    BlankName = Sample["BlankData"]["BlankName"]
    OmegaFactor= (UPDSize/SDD)*np.radians(FWHMBlank)
    Kfactor=BLPeakMax*OmegaFactor*thickness * 0.1 
    #apply calibration
    SMR_Int =  SMR_Int / (Kfactor*MSAXSCorrection) 
    SMR_Error = SMR_Error/ (Kfactor*MSAXSCorrection) 

    return {"SMR_Qvec":SMR_Qvec,
            "SMR_Int":SMR_Int,
            "SMR_Error":SMR_Error,
            "Kfactor":Kfactor,
            "OmegaFactor":OmegaFactor,
            "BlankName":BlankName,
            "thickness":thickness,
            "units":"[cm2/cm3]"
            }

def getBlankFlyscan():
      # need to get proper info from tiled on last collected Blank
      # for now fake by uusing defalut from test data
      # Then either get data from file, if exist or reduce, append, and return. 
      # We need the BL_QRS and calibration data.
    BlankPath="C:/Users/ilavsky/Documents/GitHub/Matilda/TestData/TestSet/02_21_Megan_usaxs" 
    BlankFile="HeaterBlank_0060.h5"
    # Open the HDF5 file in read/write mode
    location = 'entry/blankData/'
    with h5py.File(BlankPath+'/'+BlankFile, 'r+') as hdf_file:
            # Check if the group 'location' exists, if yes, bail out as this is all needed. 
            # if deleteExisting:
            #     # Delete the group
            #     del hdf_file[location]
            #     print("Deleted existing group 'entry/blankData'.")

            if location in hdf_file:
                # exists, so lets reuse the data from the file
                Blank = dict()
                Blank = load_dict_from_hdf5(hdf_file, location)
                print("Used existing data")
                return Blank
            else:
                Blank = dict()
                Blank["RawData"]=importFlyscan(BlankPath, BlankFile)         #import data
                BlTransCounts = Blank['RawData']['metadata']['trans_pin_counts']
                BlTransGain = Blank['RawData']['metadata']['trans_pin_gain']
                BlI0Counts = Blank['RawData']['metadata']['trans_I0_counts']
                BlI0Gain = Blank['RawData']['metadata']['trans_I0_gain']
                Blank["BlankData"]= calculatePD_Fly(Blank)                  # Creates PD_Intensity with corrected gains and background subtraction
                Blank["BlankData"].update({"BlankName":BlankFile})          # add the name of the blank file
                Blank["BlankData"].update({"BlTransCounts":BlTransCounts})  # add the BlTransCounts
                Blank["BlankData"].update({"BlTransGain":BlTransGain})      # add the BlTransGain
                Blank["BlankData"].update({"BlI0Counts":BlI0Counts})        # add the BlI0Counts
                Blank["BlankData"].update({"BlI0Gain":BlI0Gain})            # add the BlTransGain
                Blank["BlankData"].update(calculatePDError(Blank, isBlank=True))          # Calculate UPD error, mostly the same as in Igor                
                Blank["BlankData"].update(beamCenterCorrection(Blank,useGauss=0, isBlank=True)) #Beam center correction
                Blank["BlankData"].update(smooth_r_data(Blank["BlankData"]["PD_intensity"],     #smooth data data
                                                        Blank["BlankData"]["Q_array"], 
                                                        Blank["BlankData"]["PD_range"], 
                                                        Blank["BlankData"]["PD_error"], 
                                                        Blank["RawData"]["TimePerPoint"],
                                                        replaceNans=True )) 
                # we need to return just the BlandData part 
                BlankData=dict()
                BlankData=Blank["BlankData"]
                # Create the group and dataset for the new data inside the hdf5 file for future use. 
                #save_dict_to_hdf5(Blank, location, hdf_file)
                #print("Appended new data to 'entry/blankData'.")
                return BlankData



def calculatePDError(Sample, isBlank=False):
    #OK, another incarnation of the error calculations...
    UPD_array = Sample["RawData"]["UPD_array"]
    # USAXS_PD = Sample["reducedData"]["PD_intensity"]
    MeasTimeCts = Sample["RawData"]["TimePerPoint"]
    Frequency=1e6   #this is frequency of clock fed into mca1
    MeasTime = MeasTimeCts/Frequency    #measurement time in seconds per point
    if isBlank:
        UPD_gains=Sample["BlankData"]["UPD_gains"]
        UPD_bkgErr = Sample["BlankData"]["UPD_bkgErr"]    
    else:
        UPD_gains=Sample["reducedData"]["UPD_gains"]
        UPD_bkgErr = Sample["reducedData"]["UPD_bkgErr"]    

    Monitor = Sample["RawData"]["Monitor"]
    I0AmpGain=Sample["RawData"]["metadata"]["I0AmpGain"]
    VToFFactor = Sample["RawData"]["VToFFactor"]/10      #this is mca1 frequency, HDF5 writer 1.3 and above needs /10 
    SigmaUSAXSPD=np.sqrt(UPD_array*(1+0.0001*UPD_array))		#this is our USAXS_PD error estimate, Poisson error + 1% of value
    SigmaPDwDC=np.sqrt(SigmaUSAXSPD**2+(MeasTime*UPD_bkgErr)**2) #This should include now measured error for background
    SigmaPDwDC=SigmaUSAXSPD/(Frequency*UPD_gains)
    A=(UPD_array)/(VToFFactor[0]*UPD_gains)		#without dark current subtraction
    SigmaMonitor= np.sqrt(Monitor)		            #these calculations were done for 10^6 
    ScaledMonitor = Monitor
    SigmaRwave=np.sqrt((A**2 * SigmaMonitor**4)+(SigmaPDwDC**2 * ScaledMonitor**4)+((A**2 + SigmaPDwDC**2) * ScaledMonitor**2 * SigmaMonitor**2))
    SigmaRwave=SigmaRwave/(ScaledMonitor*(ScaledMonitor**2-SigmaMonitor**2))
    SigmaRwave=SigmaRwave * I0AmpGain			#fix for use of I0 gain here, the numbers were too low due to scaling of PD by I0AmpGain
    PD_error=SigmaRwave/4		#2025-04 these values are simply too large on new APS-U USAXS instrument
    result = {"PD_error":PD_error}
    return result



def test_matildaLocal():

    Sample = dict()
    #does the file exists?
    # e = os.path.isfile("C:/Users/ilavsky/Documents/GitHub/Matilda/TestData/USAXS.h5")
    # if not e:
    #     print("File not found")
    #     return
    # else:
    #     print("File found")
    #open the file
    Sample = reduceFlyscan("C:/Users/ilavsky/Documents/GitHub/Matilda/TestData/TestSet/02_21_Megan_usaxs","PPOH_25C_2_0068.h5",deleteExisting=True)    
    Q_array = Sample["reducedData"]["Q_array"]
    UPD = Sample["reducedData"]["PD_intensity"]
    Error = Sample["reducedData"]["PD_error"]
 

if __name__ == "__main__":
    #test_matilda()
    test_matildaLocal()