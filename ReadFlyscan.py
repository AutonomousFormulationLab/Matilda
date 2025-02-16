#this file will import data from numpy flyscan file


import h5py
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pprint as pp


## support stuff here
# Function to recursively read a group and store its datasets in a dictionary
def read_group_to_dict(group):
    data_dict = {}
    for key, item in group.items():
        if isinstance(item, h5py.Dataset):
            # Read the dataset
            data = item[()]
            # Check if the dataset is an array with a single element
            if data.size == 1:
                # Convert to a scalar (number or string)
                data = data.item()
            data_dict[key] = data
        elif isinstance(item, h5py.Group):
            # If the item is a group, recursively read its contents
            data_dict[key] = read_group_to_dict(item)
    return data_dict

# Function to check if all arrays have the same length
def check_arrays_same_length(*arrays):
    lengths = [arr.size for arr in arrays]  # Get the size of each array
    if len(set(lengths)) != 1:  # Check if all lengths are the same
        raise ValueError("Not all arrays have the same length.")
    #else:
        #print("All arrays have the same length.")


# Gaussian function
def gaussian(x, a, x0, sigma):
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))


## main code here
def ImportFlyscan(filename):
    # Open the HDF5 file and read its content, parse content in numpy arrays and dictionaries
    with h5py.File(filename, 'r') as file:
        #read various data sets
            #figure out how many points are in AR angles, this has 1 more point that mca data, usually
        dataset = file['/entry/flyScan/AR_PulsePositions'] 
        ARangles =  np.ravel(np.array(dataset))
        num_elements = ARangles.size - 1 
        ARangles= ARangles[-num_elements:]
        #time per point
        dataset = file['/entry/flyScan/mca1'] 
        TimePerPoint = np.ravel(np.array(dataset))  #[-num_elements:]
        #I0 - Monitor
        dataset = file['/entry/flyScan/mca2'] 
        Monitor = np.ravel(np.array(dataset))   #[-num_elements:]
        #UPD
        dataset = file['/entry/flyScan/mca3'] 
        UPD_array = np.ravel(np.array(dataset)) [-num_elements:]
        #Arrays for gain changes
        dataset = file['/entry/flyScan/changes_DDPCA300_ampGain'] 
        AmpGain = np.ravel(np.array(dataset))
        dataset = file['/entry/flyScan/changes_DDPCA300_ampReqGain'] 
        AmpReqGain = np.ravel(np.array(dataset))
        dataset = file['/entry/flyScan/changes_DDPCA300_mcsChan'] 
        Channel = np.ravel(np.array(dataset))    
        #metadata
        metadata_group = file['/entry/metadata']
        metadata_dict = read_group_to_dict(metadata_group)
        #Instrument
        instrument_group = file['/entry/instrument']
        instrument_dict = read_group_to_dict(instrument_group)

    #print(instrument_dict["monochromator"]["wavelength"])


    # Call the function with your arrays
    check_arrays_same_length(ARangles, TimePerPoint, Monitor, UPD_array)
    #Package these results into distionary
    data_dict = {"ARangles":ARangles, 
                "TimePerPoint": TimePerPoint, 
                "Monitor":Monitor, 
                "UPD_array": UPD_array,
                "AmpGain": AmpGain,
                "Channel": Channel,
                "AmpReqGain": AmpReqGain,
                "metadata": metadata_dict,
                "instrument": instrument_dict,
                }
    
    return data_dict
    
def CorrectUPDGains(data_dict):
    # create the gains array and corrects UPD for it.
    # Masks deadtimes and raneg changes
    # get the needed data from dictionary
    ARangles = data_dict["RawData"]["ARangles"]
    AmpGain = data_dict["RawData"]["AmpGain"]
    AmpReqGain = data_dict["RawData"]["AmpReqGain"]
    Channel = data_dict["RawData"]["Channel"]
    metadata_dict = data_dict["RawData"]["metadata"]
    instrument_dict = data_dict["RawData"]["instrument"]
    UPD_array = data_dict["RawData"]["UPD_array"]
    TimePerPoint = data_dict["RawData"]["TimePerPoint"]
    
    # Create Gains arrays - one for requested and one for real
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
    for i in range(0, len(GainsIndx)-1, 1):
        if np.isnan(GainsIndx[i]):
            continue
        gainName = 'DDPCA300_gain'+str(int(GainsIndx[i]))
        Gains[i] = metadata_dict[gainName]

    #mask amplifier dead times. This is done by comparing table fo deadtimes from metadata with times after range change. 
    Frequency=1e6
    TimeInSec = TimePerPoint/Frequency
    print("Exp. time :", sum(TimeInSec))
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

    #Correct UPD for gains
    UPD_corrected = (UPD_array)/(Gains) 
    
#    # Plot ydata against xdata
#     plt.figure(figsize=(10, 6))
#     plt.plot(UPD_corrected)
#     plt.title('Plot of Data from /entry/flyScan/AmpGain_array')
#     plt.xlabel('Index')
#     plt.ylabel('Value')
#     plt.grid(True)
#     plt.show()
    
    #UPD_array_log=np.log(UPD_array)
    result = {"UPD":UPD_corrected}
    return result





def BeamCenterCorrection(data_dict):
    # Find Peak center and create Q vector.
    #RawData=data_dict["rawData"]
    #ReducedData = data_dict["ReducedData"]
    ARangles = data_dict["RawData"]["ARangles"]
    metadata_dict = data_dict["RawData"]["metadata"]
    instrument_dict = data_dict["RawData"]["instrument"]
    UPD_array = data_dict["ReducedData"]["UPD"]

    # Remove NaN values from both xdata and ydata
    nan_mask = ~np.isnan(ARangles) & ~np.isnan(UPD_array)
    xdata_clean = ARangles[nan_mask]
    ydata_clean = UPD_array[nan_mask]

    # Find the threshold for the top 50% of UPD_array
    threshold = np.percentile(ydata_clean, 50)

    print(f"Threshold for the top 50% of UPD_array: {threshold}")

    # Filter the data to only include the top 50%
    mask = ydata_clean >= threshold
    xdata_filtered = xdata_clean[mask]
    ydata_filtered = ydata_clean[mask]

    # Initial guess for the parameters: amplitude, mean, and standard deviation
    initial_guess = [np.max(ydata_filtered), xdata_filtered[np.argmax(ydata_filtered)], 1]

    # Fit the Gaussian function to the filtered data
    popt, _ = curve_fit(gaussian, xdata_filtered, ydata_filtered, p0=initial_guess)

    # Extract the fitted parameters
    amplitude, x0, sigma = popt

    # Calculate the FWHM
    fwhm = 2 * np.abs(np.sqrt(2 * np.log(2)) * sigma)

    # Calculate the predicted y values using the fitted parameters
    y_pred = gaussian(xdata_filtered, *popt)

    # Calculate the residuals
    residuals = ydata_filtered - y_pred

    # Calculate the chi-square
    # If you have measurement errors, replace 1 with the variance of ydata_filtered
    chi_square = np.sum((residuals**2) / 1)

    #Make wave vector
    Q_array = np.full(UPD_array.shape, 0)
    #AR_center = metadata_dict["AR_center"]
    wavelength = instrument_dict["monochromator"]["wavelength"]
    Q_array = -1*(4*np.pi*np.sin(np.radians(ARangles-x0)/2)/wavelength)
    #Q_array_log = np.sign(Q_array)*np.log(np.abs(Q_array))

    results = {"Q_array":Q_array,
            "Chi-Square":chi_square,
            "Center":x0,
            "Maximum":amplitude,
            "FWHM":fwhm    
            }
    return results


    # Print the metadata
    # Now, metadata_dict contains the data from the /entry/metadata group
    #pp.pprint(instrument_dict)

    # # Plot the data
    # plt.figure(figsize=(10, 6))
    # plt.plot(UPD)
    # plt.title('Plot of Data from /entry/flyScan/AmpGain_array')
    # plt.xlabel('Index')
    # plt.ylabel('Value')
    # plt.grid(True)
    # plt.show()

def PlotResults(data_dict):
    # Find Peak center and create Q vector.
    Q_array = data_dict["ReducedData"]["Q_array"]
    UPD = data_dict["ReducedData"]["UPD"]
    
    # Plot ydata against xdata
    plt.figure(figsize=(6, 12))
    plt.plot(Q_array, UPD, marker='o', linestyle='-')  # You can customize the marker and linestyle
    plt.title('Plot of UPD vs. Q')
    plt.xlabel('log(Q) [1/A]')
    plt.ylabel('UPD')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(True)
    plt.show()



#if __file__ == "__main__":
Sample = dict()
Sample["RawData"]=ImportFlyscan("USAXS.h5")
#pp.pprint(Sample)
Sample["ReducedData"]= CorrectUPDGains(Sample)
Sample["ReducedData"].update(BeamCenterCorrection(Sample))
#pp.pprint(Sample["ReducedData"])
PlotResults(Sample)

  
    
    