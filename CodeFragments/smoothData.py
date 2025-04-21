'''
This is copy of Igor smooth data for our USAXS data
This is complciated procedure which handles log-log data and manages range changes
'''


import numpy as np
from scipy.optimize import curve_fit

def smooth_r_data(intensity, qvector, pd_range, r_error, meas_time):
    # Smoothing times for different ranges
    rwave_smooth_times = [0, 0, 0.01, 0.03, 0.06]   # these are [in sec] values for USAXS on 4/20/2025

    # Logarithm of intensity
    temp_int_log = np.log(intensity)
    smooth_intensity = np.copy(temp_int_log)

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

        if meas_time[i] > tmp_time:
            smooth_intensity[i] = temp_int_log[i]
        else:
            start_points = int(np.ceil(tmp_time / meas_time[i])) + 1
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
    return  {"PD_intensity":intensity,
              "PD_error":r_error} 

# Example usage:
# intensity, qvector, pd_range, r_error, meas_time = load_your_data()
# smoothed_intensity, updated_r_error = smooth_r_data(intensity, qvector, pd_range, r_error, meas_time)
