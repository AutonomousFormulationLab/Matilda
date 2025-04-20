'''
Here we develop new code which then moves to proper package
'''
import h5py
import numpy as np
#from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pprint as pp
#from supportFunctions import read_group_to_dict, filter_nested_dict, check_arrays_same_length
import os
from convertUSAXS import rebinData
#from hdf5code import save_dict_to_hdf5, load_dict_from_hdf5
from convertUSAXS import importFlyscan, calculatePD_Fly, beamCenterCorrection



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
                Sample["RawData"]=importFlyscan(path, filename)         #import data
                Sample["ReducedData"]= calculatePD_Fly(Sample)       # Correct gains
                Sample["ReducedData"].update(calculatePDError(Sample))         # Calculate UPD error                
                Sample["ReducedData"].update(beamCenterCorrection(Sample,useGauss=0)) #Beam center correction
                # TODO: need to create errors wave here
                # TODO: Blank/background subtraction
                # TODO:     figure out proper blank and if needed, reduce to BL_QRS
                # TODO: calibration
                # TODO: fix rebinning for 3 input waves returning 4 waves with dQ
                Sample["ReducedData"].update(rebinData(Sample))         #Rebin data
                # TODO: desmearing here
                # Create the group and dataset for the new data inside the hdf5 file for future use. 
                #save_dict_to_hdf5(Sample, location, hdf_file)
                #print("Appended new data to 'entry/displayData'.")
                return Sample




def calculatePDError(Sample):
    #OK, another incarnation of the error calculations...
    UPD_array = Sample["RawData"]["UPD_array"]
    # USAXS_PD = Sample["ReducedData"]["PD_intensity"]
    #MeasTime = Sample["RawData"]["TimePerPoint"]
    UPD_gains=Sample["ReducedData"]["UPD_gains"]
    Frequency=1e6   #this is frequency of clock fed into mca1
    Monitor = Sample["RawData"]["Monitor"]
    I0AmpGain=Sample["RawData"]["metadata"]["I0AmpGain"]
    VToFFactor = Sample["RawData"]["VToFFactor"]/10      #this is mca1 frequency, HDF5 writer 1.3 and above needs /10 
    SigmaUSAXSPD=np.sqrt(UPD_array*(1+0.0001*UPD_array))		#this is our USAXS_PD error estimate, Poisson error + 1% of value
	#SigmaPDwDC=np.sqrt(SigmaUSAXSPD^2+(MeasTime*ErrorParameters[UPD_gains-1])^2) #This should be measured error for background
    SigmaPDwDC=SigmaUSAXSPD/(Frequency*UPD_gains)
    A=(UPD_array)/(VToFFactor[0]*UPD_gains)		#without dark current subtraction
    SigmaMonitor= np.sqrt(Monitor)		            #these calculations were done for 10^6 
    ScaledMonitor = Monitor
    SigmaRwave=np.sqrt((A**2 * SigmaMonitor**4)+(SigmaPDwDC**2 * ScaledMonitor**4)+((A**2 + SigmaPDwDC**2) * ScaledMonitor**2 * SigmaMonitor**2))
    SigmaRwave=SigmaRwave/(ScaledMonitor*(ScaledMonitor**2-SigmaMonitor**2))
    SigmaRwave=SigmaRwave * I0AmpGain			#fix for use of I0 gain here, the numbers were too low due to scaling of PD by I0AmpGain
    PD_error=SigmaRwave/5		#2025-04 these values are simply too large on new APS-U USAXS instrument
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
    Sample = reduceFlyscan("C:/Users/ilavsky/Documents/GitHub/Matilda/TestData","USAXS.h5",deleteExisting=True)    
    Q_array = Sample["ReducedData"]["Q_array"]
    UPD = Sample["ReducedData"]["PD_intensity"]
    Error = Sample["ReducedData"]["PD_error"]
 

if __name__ == "__main__":
    #test_matilda()
    test_matildaLocal()