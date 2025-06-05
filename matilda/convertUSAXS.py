''' 

    Import data from flyscan or step scan hdf5 file
    and convert to QRdata

    First check if group root:DisplayData exists, if not, load data, process and save for future use. 
    If it exists, load data from it.
    This is to save time for future use.

    Main routines: 
    reduceStepScanToQR
    reduceFlyscanToQR
'''

import h5py
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import pprint as pp
from supportFunctions import read_group_to_dict, filter_nested_dict, check_arrays_same_length
from supportFunctions import gaussian, modifiedGauss
import os
from rebinData import rebin_QRSdata
from hdf5code import save_dict_to_hdf5, load_dict_from_hdf5



## Flyscan main code here

## importFlyscan loads data from flyscan NX file. It should be same for QR pass as well as for calibrated data processing. 
def importFlyscan(path, filename):
    # Open the HDF5 file and read its content, parse content in numpy arrays and dictionaries
    with h5py.File(path+"/"+filename, 'r') as file:
        #read various data sets
        #figure out how many points are in AR angles, this has 1 more point that mca data, usually
        dataset = file['/entry/flyScan/AR_PulsePositions'] 
        ARangles =  np.ravel(np.array(dataset))
        num_elements = ARangles.size - 1 
        ARangles= ARangles[-num_elements:]
        #time per point
        dataset = file['/entry/flyScan/mca1'] 
        TimePerPoint = np.ravel(np.array(dataset))  [-num_elements:]
        #I0 - Monitor
        dataset = file['/entry/flyScan/mca2'] 
        Monitor = np.ravel(np.array(dataset))   [-num_elements:]
        #UPD
        dataset = file['/entry/flyScan/mca3'] 
        UPD_array = np.ravel(np.array(dataset)) [-num_elements:]
        #Arrays for UPD gain changes
        dataset = file['/entry/flyScan/changes_DDPCA300_ampGain'] 
        AmpGain = np.ravel(np.array(dataset))
        dataset = file['/entry/flyScan/changes_DDPCA300_ampReqGain'] 
        AmpReqGain = np.ravel(np.array(dataset))
        dataset = file['/entry/flyScan/changes_DDPCA300_mcsChan'] 
        Channel = np.ravel(np.array(dataset))            
        dataset = file['/entry/flyScan/mca_clock_frequency'] 
        vTof = np.ravel(np.array(dataset))    
        #metadata
        keys_to_keep = ['AR_center', 'ARenc_0', 'DCM_energy', 'DCM_theta', 'I0AmpGain','detector_distance',
                        'timeStamp',
                        'trans_pin_counts','trans_pin_gain','trans_pin_time','trans_I0_counts','trans_I0_gain',
                        'UPDsize', 'trans_I0_counts', 'trans_I0_gain', 'upd_bkg0', 'upd_bkg1','upd_bkg2','upd_bkg3',
                        'upd_bkgErr0','upd_bkgErr1','upd_bkgErr2','upd_bkgErr3','upd_bkgErr4','upd_bkg_err0',
                        'upd_bkg4','DDPCA300_gain0','DDPCA300_gain1','DDPCA300_gain2','DDPCA300_gain3','DDPCA300_gain4',
                        'upd_amp_change_mask_time0','upd_amp_change_mask_time1','upd_amp_change_mask_time2','upd_amp_change_mask_time3','upd_amp_change_mask_time4',
                    ]
        metadata_group = file['/entry/metadata']
        metadata_dict = read_group_to_dict(metadata_group)
        metadata_dict = filter_nested_dict(metadata_dict, keys_to_keep)
        #Instrument
        keys_to_keep = ['monochromator', 'energy', 'wavelength']
        instrument_group = file['/entry/instrument']
        instrument_dict = read_group_to_dict(instrument_group)
        instrument_dict = filter_nested_dict(instrument_dict, keys_to_keep)
        # sample
        sample_group = file['/entry/sample']
        sample_group = read_group_to_dict(sample_group)

    # Call the function with your arrays
    check_arrays_same_length(ARangles, TimePerPoint, Monitor, UPD_array)
    #Package these results into dictionary
    data_dict = {"Filename": os.path.splitext(filename)[0],
                "ARangles":ARangles, 
                "TimePerPoint": TimePerPoint, 
                "Monitor":Monitor, 
                "UPD_array": UPD_array,
                "AmpGain": AmpGain,
                "Channel": Channel,
                "VToFFactor": vTof,
                "AmpReqGain": AmpReqGain,
                "metadata": metadata_dict,
                "instrument": instrument_dict,
                "sample": sample_group,
                }
    
    return data_dict

# 
def calculatePD_Fly(data_dict):
        # create the gains array and corrects UPD for it.
        # Masks deadtimes and range changes
        # get the needed data from dictionary
    ARangles = data_dict["RawData"]["ARangles"]
    AmpGain = data_dict["RawData"]["AmpGain"]
    AmpReqGain = data_dict["RawData"]["AmpReqGain"]
    Channel = data_dict["RawData"]["Channel"]
    metadata_dict = data_dict["RawData"]["metadata"]
    instrument_dict = data_dict["RawData"]["instrument"]
    UPD_array = data_dict["RawData"]["UPD_array"]
    TimePerPoint = data_dict["RawData"]["TimePerPoint"]
    Monitor = data_dict["RawData"]["Monitor"]
    VToFFactor = data_dict["RawData"]["VToFFactor"]

    
        # Create Gains arrays - one for requested and one for real
    I0AmpGain = metadata_dict["I0AmpGain"]
    num_elements = UPD_array.size 
    AmpGain_array = np.full(num_elements, AmpGain[len(AmpGain)-1])
    AmpGainReq_array = np.full(num_elements,AmpGain[len(AmpReqGain)-1])

        # Iterate over the Channel array to get index pairs
    for i in range(0, len(Channel)-2, 1):
        start_index = int(Channel[i])
        end_index = int(Channel[i + 1])
            
        # Ensure indices are within bounds
        if start_index < 0 or end_index > len(AmpGain_array) or start_index > end_index:
            raise ValueError("Invalid index range in Channel array.")
        
        # Fill the new array with AmpGain and AmpGainReq values between the two indices
        if(start_index!=end_index):
            AmpGain_array[start_index:end_index] = AmpGain[i]    
            AmpGainReq_array[start_index:end_index] = AmpReqGain[i]
        else:
            AmpGain_array[start_index] = AmpGain[i]    
            AmpGainReq_array[start_index] = AmpReqGain[i]

        # Create a new array res with the same shape as AmpGain_array, initialized with NaN
    GainsIndx = np.full(AmpGain_array.shape, np.nan)
    Gains = np.full(AmpGain_array.shape, np.nan)
    updBkg = np.full(AmpGain_array.shape, np.nan)
    updBkgErr = np.full(AmpGain_array.shape, np.nan)

        # Use a boolean mask to find where two arrays agree
    mask = AmpGain_array == AmpGainReq_array

        # Set the values in res where two agree
    GainsIndx[mask] = AmpGain_array[mask]
        #set to Nan also points in channel array that are not in mask
    for i in range(0, len(Channel)-2, 1):
        s = int(Channel[i])
        GainsIndx[s] = np.nan

        #next, replace the values in Gains array with values looked up from metadata dictionary
        #for now, lets look only for DDPCA300_gain+"gainNumber"
        # also need to create gain measureds background 'upd_bkg0', 'upd_bkg1','upd_bkg2','upd_bkg3', 'upd_bkg4'
    for i in range(0, len(GainsIndx)-1, 1):
        if np.isnan(GainsIndx[i]):
            continue
        gainName = 'DDPCA300_gain'+str(int(GainsIndx[i]))
        Gains[i] = metadata_dict[gainName]
        updBkgName = 'upd_bkg'+str(int(GainsIndx[i]))
        updBkg[i] =  metadata_dict[updBkgName]
        try:
            updBkgErrName = 'upd_bkgErr'+str(int(GainsIndx[i]))
            updBkgErr[i] =  metadata_dict[updBkgErrName]
        except:
            updBkgErrName = 'upd_bkg_err'+str(int(GainsIndx[i]))     #typo in Flyscan schema below 1.3 (before June 2025) 
            updBkgErr[i] =  metadata_dict[updBkgErrName]

        #mask amplifier dead times. This is done by comparing table fo deadtimes from metadata with times after range change. 
    Frequency=VToFFactor[0]/10   #this is frequency of clock fed ito mca1/10 for HDF5 writer 1.3 and higher
    TimeInSec = TimePerPoint/Frequency
    #print("Exp. time :", sum(TimeInSec))
    for i in range(0, len(Channel)-1, 1):
        startPnt=Channel[i]
        deadtimeName = 'upd_amp_change_mask_time'+str(int(AmpReqGain[i]))
        deadtime = metadata_dict[deadtimeName]
        elapsed = 0
        indx = int(startPnt)
        while elapsed < deadtime:
            elapsed+=TimeInSec[indx]
            Gains[indx]=np.nan
            #print("Index is:", indx)
            #print("elapsed time is:",elapsed) 
            indx += 1

        #Correct UPD for gains and monitor counts and amplifier gain. 
        # Frequency=1e6  #this is to keep in sync with Igor code. 
    PD_Intensity = ((UPD_array-TimeInSec*updBkg)/(Frequency*Gains)) / (Monitor/I0AmpGain)  
    PD_error = 0.01*PD_Intensity   #this is fake error for QR conversion
    result = {"Intensity":PD_Intensity,
              "Error":PD_error,
              "PD_range":GainsIndx,
              "UPD_gains":Gains,
              "UPD_bkgErr":updBkgErr}
    return result

def rebinData(data_dict,num_points=200, isSMRData=False):
    # Rebin data to 200+peak area points.
    if isSMRData:
        SMR_Int = data_dict["CalibratedData"]["SMR_Int"]
        SMR_Qvec = data_dict["CalibratedData"]["SMR_Qvec"]
        SMR_Error = data_dict["CalibratedData"]["SMR_Error"]
        SMR_QvecNew, SMR_IntNew, SMR_ErrorNew, SMR_dQ = rebin_QRSdata(SMR_Qvec, SMR_Int,SMR_Error, num_points)
        results = {"SMR_Qvec":SMR_QvecNew,
                "SMR_Int":SMR_IntNew,
                "SMR_Error":SMR_ErrorNew,
                "SMR_dQ":SMR_dQ,
                }    
    else:
        Q_array = data_dict["reducedData"]["Q"]
        R_array = data_dict["reducedData"]["Intensity"]
        S_array = data_dict["reducedData"]["Error"]
        Q_arrayNew, R_arrayNew, S_arrayNew = rebin_QRSdata(Q_array, R_array,S_array, num_points)
        results = {"Q":Q_arrayNew,
                "Intensity":R_arrayNew,
                "Error":S_arrayNew  
                }
    return results


## Stepscan main code here
def ImportStepScan(path, filename):
    # Open the HDF5 file and read its content, parse content in numpy arrays and dictionaries
    with h5py.File(path+"/"+filename, 'r') as file:
        #read various data sets
        #AR angle
        dataset = file['/entry/data/a_stage_r'] 
        ARangles = np.ravel(np.array(dataset))         
        #time per point
        dataset = file['/entry/data/seconds'] 
        TimePerPoint = np.ravel(np.array(dataset))         
        # I0 gain
        dataset = file['/entry/data/I0_autorange_controls_gain'] 
        I0gain = np.ravel(np.array(dataset)) 
        #I0 - Monitor
        dataset = file['/entry/data/I0'] 
        Monitor = np.ravel(np.array(dataset))  
        #UPD
        dataset = file['/entry/data/UPD'] 
        UPD_array = np.ravel(np.array(dataset))
        #Arrays for gain changes
        dataset = file['/entry/data/upd_autorange_controls_gain'] 
        AmpGain = np.ravel(np.array(dataset))
            #dataset = file['/entry/data/upd_autorange_controls_reqrange'] 
            #AmpReqGain = np.ravel(np.array(dataset))       #this contains only 0 values, useless... 
        #metadata
        keys_to_keep = ['SAD_mm', 'SDD_mm', 'thickness', 'title', 'useSBUSAXS',
                        'intervals', 
                    ]
        metadata_group = file['/entry/instrument/bluesky/metadata']
        metadata_dict = read_group_to_dict(metadata_group)     
        metadata_dict = filter_nested_dict(metadata_dict, keys_to_keep)
        #Instrument
        instrument_group = file['/entry/instrument/monochromator']
        instrument_dict = read_group_to_dict(instrument_group)        
        #Sample
        sample_group = file['/entry/sample']
        sample_dict = read_group_to_dict(sample_group)


    # Call the function with your arrays
    check_arrays_same_length(ARangles, TimePerPoint, Monitor, UPD_array)
    #Package these results into dictionary
    data_dict = {"Filename": os.path.splitext(filename)[0],
                "ARangles":ARangles, 
                "TimePerPoint": TimePerPoint, 
                "Monitor":Monitor, 
                "UPD_array": UPD_array,
                "AmpGain": AmpGain,
                "I0gain": I0gain,
                "sample": sample_dict,
                "metadata": metadata_dict,
                "instrument": instrument_dict,
                }
    
    return data_dict
    
def CorrectUPDGainsStep(data_dict):
        # here we will multiply UPD by gain and divide by monitor corrected for its gain.
        # get the needed data from dictionary
    AmpGain = data_dict["RawData"]["AmpGain"]
    UPD_array = data_dict["RawData"]["UPD_array"]
    Monitor = data_dict["RawData"]["Monitor"]
    I0gain = data_dict["RawData"]["I0gain"]
        # for some  reason, the AmpGain is shifted by one value so we need to duplicate the first value and remove end value. 
    first_value = AmpGain[0]
    AmpGain = np.insert(AmpGain, 0, first_value)
    AmpGain = AmpGain[:-1]
                # change gain masking may not any be necessary... 
                # # need to remove points where gain changes
                # # Find indices where the change occurs
                # change_indices = np.where(np.diff(AmpGain) != 0)[0]
                # change_indices = change_indices +1
                # # fix range changes
                # #Correct UPD for gains so we can find max value loaction
                # UPD_temp = (UPD_array*I0gain)/(AmpGain*Monitor)
                # #remove renage chanegs on thsi array
                # UPD_temp[change_indices] = np.nan
                # # now locate location of max value in UPD_array
                # max_index = np.nanargmax(UPD_temp)
                # # we need to limit change_indices to values less than the location of maximum (before peak) = max_index
                # # this removes the range changes only to before the peak location, does nto seem to work, really
                # #change_indices = change_indices[change_indices < max_index]
                # # Create a copy of the array to avoid modifying the original
                # AmpGain_new = AmpGain.astype(float)                 # Ensure the array can hold NaN values
                # # Set the point before each range change to NaN
                # if len(change_indices) > 0:
                #     AmpGain_new[change_indices] = np.nan
                #Correct UPD for gains with points we  want removed set to Nan

        #Correct UPD for gains and monitor
    UPD_corrected = (UPD_array*I0gain)/(AmpGain*Monitor)
    result = {"Intensity":UPD_corrected}
    return result

## Common steps go here
def beamCenterCorrection(data_dict, useGauss=1, isBlank=False):
    # Find Peak center and create Q vector.
        #RawData=data_dict["rawData"]
        #reducedData = data_dict["reducedData"]
    ARangles = data_dict["RawData"]["ARangles"]
    instrument_dict = data_dict["RawData"]["instrument"]
    if isBlank:
        UPD_array = data_dict["BlankData"]["Intensity"]
    else:
        UPD_array = data_dict["reducedData"]["Intensity"]

        #plt.figure(figsize=(6, 12))
        #plt.plot(ARangles, UPD_array, marker='o', linestyle='-')  # You can customize the marker and linestyle
        # Remove NaN values from both xdata and ydata
    nan_mask = ~np.isnan(ARangles) & ~np.isnan(UPD_array)
    xdata_clean = ARangles[nan_mask]
    ydata_clean = UPD_array[nan_mask]
        #plt.plot(xdata_clean, ydata_clean, marker='o', linestyle='-') 
        # Find the threshold for the top ~40% of UPD_array
    threshold = np.max(ydata_clean)/2.3
        #pp.pprint(threshold)
        #print(f"Threshold for the top 40% of UPD_array: {threshold}")
        # Filter the data to only include the top 40%
    mask = ydata_clean >= threshold
    xdata_filtered = xdata_clean[mask]
    ydata_filtered = ydata_clean[mask]
        #plt.plot(xdata_filtered, ydata_filtered, marker='o', linestyle='-')

    
    if useGauss:
        # Initial guess for the parameters: amplitude, mean, and standard deviation, 2 fgor d parameter
        
        initial_guess = [np.max(ydata_filtered), xdata_filtered[np.argmax(ydata_filtered)], 0.0001]
        # Fit the Gaussian function to the filtered data
        popt, _ = curve_fit(gaussian, xdata_filtered, ydata_filtered, p0=initial_guess)
        # Extract the fitted parameters
        amplitude, x0, sigma = popt
        # Calculate the FWHM
        fwhm = 2 * np.abs(np.sqrt(2 * np.log(2)) * sigma)        
        # Calculate the predicted y values using the fitted parameters
        y_pred = gaussian(xdata_filtered, *popt)       
    else:
        initial_guess = [np.max(ydata_filtered), xdata_filtered[np.argmax(ydata_filtered)], 0.0001, 1.98]
        #print(initial_guess)
        popt, _ = curve_fit(modifiedGauss, xdata_filtered, ydata_filtered, p0=initial_guess)

        # Extract the fitted parameters
        amplitude, x0, sigma, dparam = popt
      
        # Calculate the FWHM
        # Calculate the half maximum
        half_max = amplitude / 2

        # but next calculation needs to be done over larger q range
        threshold = amplitude/3
        mask = ydata_clean >= threshold
        xdata_calc = xdata_clean[mask]
        ydata_calc = ydata_clean[mask]
     
        # Calculate the predicted y values using the fitted parameters
        ydata_calc = modifiedGauss(xdata_calc, *popt)

        # Find where the array crosses the half maximum
        crossings = np.where((ydata_calc[:-1] < half_max) & (ydata_calc[1:] >= half_max) |
                     (ydata_calc[:-1] >= half_max) & (ydata_calc[1:] < half_max))[0]

        # Calculate fractional crossing indices using linear interpolation
        indices = []
        for i in crossings:
            y1, y2 = ydata_calc[i], ydata_calc[i + 1]
            if y2 != y1:
                fractional_index = i + (half_max - y1) / (y2 - y1)
                indices.append(fractional_index)
                
        # Calculate the FWHM
        if len(indices) >= 2:
            i = int(np.floor(indices[0]))
            y1 = xdata_calc[i]
            y2 = xdata_calc[i+1]
            xdata_calcl = y1 + (indices[0] - i) * (y2 - y1)
            i = int(np.floor(indices[-1]))
            y1 = xdata_calc[i]
            y2 = xdata_calc[i+1]
            xdata_calch = y1 + (indices[-1] - i) * (y2 - y1)           
            fwhm = np.abs(xdata_calch - xdata_calcl)
        else:
            fwhm = np.nan  # FWHM cannot be determined

        # Calculate the residuals
        y_pred = modifiedGauss(xdata_filtered, *popt)

    #use above calculated y_pred to get residuals
    residuals = ydata_filtered - y_pred

        # Calculate the chi-square
        # If you have measurement errors, replace 1 with the variance of ydata_filtered
    chi_square = np.sum((residuals**2) / 1)

        #Make wave vector
    Q_array = np.full(UPD_array.shape, 0)
        #AR_center = metadata_dict["AR_center"]
    try:
        # Try to get the value using the first key
        wavelength = instrument_dict["monochromator"]["wavelength"]
    except KeyError:
        #print(instrument_dict)
        # If the first key doesn't exist, try the second key
        wavelength = instrument_dict["wavelength"]

    Q_array = -1*(4*np.pi*np.sin(np.radians(ARangles-x0)/2)/wavelength)
        #Q_array = (4*np.pi*np.sin(np.radians(ARangles-x0)/2)/wavelength)
        #Q_array_log = np.sign(Q_array)*np.log(np.abs(Q_array))
        #pp.pprint(Q_array)

    results = {"Q":Q_array,
            "Chi-Square":chi_square,
            "Center":x0,
            "Maximum":amplitude,
            "FWHM":fwhm,
            "wavelength":wavelength,    
            }
    return results


def smooth_r_data(intensity, qvector, pd_range, r_error, meas_time, replaceNans=False):
    # Smoothing times for different ranges
    rwave_smooth_times = [0, 0, 0.01, 0.03, 0.06]   # these are [in sec] values for USAXS on 4/20/2025

    # Logarithm of intensity
    temp_int_log = np.log(intensity)
    # for some cases replace nans with interpolated values
    if replaceNans:
        # Create a mask for NaNs
        nan_mask = np.isnan(temp_int_log)
        # Indices of non-NaN values
        x_non_nan = qvector[~nan_mask]
        y_non_nan = temp_int_log[~nan_mask]
        # Create an interpolation function
        interp_func = interp1d(x_non_nan, y_non_nan, kind='linear', fill_value='extrapolate')
        # Replace NaNs with interpolated values
        temp_int_log[nan_mask] = interp_func(qvector[nan_mask])


    smooth_intensity = np.copy(temp_int_log)
    meas_time_sec = meas_time/1e6       # meas_time is still frequency, need time in seconds. 

    def linear_fit(x, a, b):
        return a + b * x

    for i in range(40, len(intensity)):
        if pd_range[i] == 1:
            tmp_time = rwave_smooth_times[0]
        elif pd_range[i] == 2:
            tmp_time = rwave_smooth_times[1]
        elif pd_range[i] == 3:
            tmp_time = rwave_smooth_times[2]
        elif pd_range[i] == 4:
            tmp_time = rwave_smooth_times[3]
        else:
            tmp_time = rwave_smooth_times[4]

        if meas_time_sec[i] > tmp_time:
            smooth_intensity[i] = temp_int_log[i]
        else:
            start_points = int(np.ceil(tmp_time / meas_time_sec[i])) + 1
            end_points = start_points

            if (i - start_points) < 0:
                raise ValueError("Bad data, cannot fix this. Likely Flyscan parameters were wrong")

            if i + end_points > len(intensity) - 1:
                end_points = len(intensity) - 1 - i

            if (pd_range[i - start_points] != pd_range[i]) or (pd_range[i + end_points] != pd_range[i]):
                temp_r = temp_int_log[i - start_points:i + end_points]
                temp_q = qvector[i - start_points:i + end_points]

                if len(temp_r) > np.isnan(temp_r).sum() + 5:
                    popt, _ = curve_fit(linear_fit, temp_q, temp_r)
                    smooth_intensity[i] = linear_fit(qvector[i], *popt)
                    r_error[i] /= 3
                else:
                    smooth_intensity[i] = temp_int_log[i]
                    r_error[i] = r_error[i]
            else:
                temp_r = temp_int_log[i - start_points:i + end_points + 1]
                temp_q = qvector[i - start_points:i + end_points + 1]
                start_x = temp_q[0]
                end_x = temp_q[-1]
                area = np.trapezoid(temp_r, temp_q)
                smooth_intensity[i] = area / (end_x - start_x)
                r_error[i] = r_error[i]

    intensity = np.exp(smooth_intensity)
    return  {"Intensity":intensity,
              "Error":r_error} 

def PlotResults(data_dict):
        # Plot UPD vs Q.
    Q = data_dict["reducedData"]["Q"]
    Intensity = data_dict["reducedData"]["Intensity"]
    
        # Plot ydata against xdata
    plt.figure(figsize=(6, 12))
    plt.plot(Q, Intensity, marker='o', linestyle='-')  # You can customize the marker and linestyle
    plt.title('Plot of UPD vs. Q')
    plt.xlabel('log(Q) [1/A]')
    plt.ylabel('UPD')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(True)
    plt.show()

# TODO: remove deleteExisting=True for operations
def reduceFlyscanToQR(path, filename, deleteExisting=True):
    # Open the HDF5 file in read/write mode
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
                return Sample
            else:
                Sample = dict()
                Sample["RawData"]=importFlyscan(path, filename)         #import data
                Sample["reducedData"]= calculatePD_Fly(Sample)       # Correct gains
                Sample["reducedData"].update(beamCenterCorrection(Sample,useGauss=0)) #Beam center correction
                Sample["reducedData"].update(rebinData(Sample))         #Rebin data
                # Create the group and dataset for the new data inside the hdf5 file for future use. 
                # these are not fully reduced data, this is for web plot purpose. 
                save_dict_to_hdf5(Sample, location, hdf_file)
                print("Appended new data to 'entry/displayData'.")
                return Sample
# TODO: remove deleteExisting=True for operations
def reduceStepScanToQR(path, filename, deleteExisting=True):
  # Open the HDF5 file in read/write mode
    location = 'entry/displayData/'
    with h5py.File(path+'/'+filename, 'r+') as hdf_file:
        if deleteExisting:
            # Delete the group
            del hdf_file[location]
            print("Deleted existing group 'entry/displayData'.")    
        
        if location in hdf_file:
                # # exists, reuse existing data
                Sample = dict()
                Sample = load_dict_from_hdf5(hdf_file, location)
                return Sample
        else:
                Sample = dict()
                Sample["RawData"]=ImportStepScan(path, filename)
                Sample["reducedData"]= CorrectUPDGainsStep(Sample)
                Sample["reducedData"].update(beamCenterCorrection(Sample,useGauss=1))
                # Create the group and dataset for the new data inside the hdf5 file for future use.
                # these are not fully reduced data, this is for web plot purpose.
                save_dict_to_hdf5(Sample, location, hdf_file)
                return Sample



if __name__ == "__main__":
    #Sample = dict()
    #Sample = reduceStepScanToQR("/home/parallels/Github/Matilda/TestData","USAXS_step.h5")
    #Sample["RawData"]=ImportStepScan("/home/parallels/Github/Matilda","USAXS_step.h5")
        #pp.pprint(Sample)
        #Sample["reducedData"]= CorrectUPDGainsStep(Sample)
        #Sample["reducedData"].update(BeamCenterCorrection(Sample))
        #pp.pprint(Sample["reducedData"])
    #PlotResults(Sample)
    #flyscan
    Sample = dict()
    Sample = reduceFlyscanToQR("./TestData","USAXS.h5",deleteExisting=True)
    # Sample["RawData"]=ImportFlyscan("/home/parallels/Github/Matilda","USAXS.h5")
    # #pp.pprint(Sample)
    # Sample["reducedData"]= CorrectUPDGainsFly(Sample)
    # Sample["reducedData"].update(BeamCenterCorrection(Sample))
    # #pp.pprint(Sample["reducedData"])
    PlotResults(Sample)

  
    
    
