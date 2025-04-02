'''
    This package contains the code to subtract two SAS data sets
    It needs to operate on QRSdQ data Q1,R1,S1,dQ1 and Q2,R2,S2,dQ2
    Needs to use interpolation for Q1 values of the second data set (Q2, R2, S2)
    Needs to do linear interpolation of intensity using log-intensity values
    Needs to propagate errors through in meaningful and correct way
    dQ is simply kept for Q1
'''

import numpy as np
from scipy.interpolate import interp1d

def subtract_data(X1, Y1, E1, X2, Y2, E2):
    """
    Interpolates and subtracts the input data sets to return X1, Ydiff, and Esub.

    Parameters:
    X1 (array-like): X1 data points.
    Y1 (array-like): Y1 data points.
    E1 (array-like): Uncertainty in Y1.
    X2 (array-like): X2 data points.
    Y2 (array-like): Y2 data points.
    E2 (array-like): Uncertainty in Y2.

    Returns:
    tuple: (X1, Ydiff, Esub) where
        - X1 is the input X1 data points.
        - Ydiff is the difference between Y1 and interpolated Y2.
        - Esub is the propagated uncertainty.
    """
    # Step 1: Interpolate log(Y2) vs X2 to values of X1
    logY2 = np.log(Y2)
    logY2_interp_func = interp1d(X2, logY2, kind='linear', fill_value='extrapolate')
    logY2_interp = logY2_interp_func(X1)

    # Convert interpolated logY2 back to Y2interp
    Y2_interp = np.exp(logY2_interp)

    # Step 2: Subtract Y1 - Y2interp to obtain Ydiff
    Ydiff = Y1 - Y2_interp

    # Step 3: Linearly interpolate E2 vs X2 to have E2interp vs X1
    E2_interp_func = interp1d(X2, E2, kind='linear', fill_value='extrapolate')
    E2_interp = E2_interp_func(X1)

    # Step 4: Propagate uncertainties for subtraction to obtain Esub
    Esub = np.sqrt(E1**2 + E2_interp**2)

    # Return the three data sets: X1, Ydiff, and Esub
    return X1, Ydiff, Esub

# Example usage:
# X1, Y1, E1 = np.array([...]), np.array([...]), np.array([...])
# X2, Y2, E2 = np.array([...]), np.array([...]), np.array([...])
# X1, Ydiff, Esub = process_data(X1, Y1, E1, X2, Y2, E2)
# print("X1:", X1)
# print("Ydiff:", Ydiff)
# print("Esub:", Esub)


