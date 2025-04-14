# Convert Nika data to pyFAI using Fit2D 
'''
    Convert Nika SAS geometry parameters to pyFAI using Fit2D format.
    Needs pyFAI library. 
    Convert Nika SDD, pix size, BCX, BCY, HorTilt, verTiilt 
    First convert into Fit2D format
    Then use pyFAI to convert Fit2D to poni format.
'''

from pyFAI.geometry.fit2d import Fit2dGeometry, convert_from_Fit2d
from pprint import pprint as pp
import numpy as np

# Fit2D parameter definitions/example
# direct_dist = 200.0  # Distance from sample to detector in mm
# center_x = 1000.0    # Beam center X in pixels
# center_y = 1000.0    # Beam center Y in pixels
# tilt = 0.0           # Tilt angle in degrees
# tilt_plan_rotation = 0.0  # Tilt plan rotation in degrees
# pixel_x = 75         # Pixel size in um (X direction)
# pixel_y = 75         # Pixel size in um (Y direction)
# spline_file = None   # Path to spline file if any
# detector = None      # Detector object if any
# wavelength = 0.59    # Wavelength in Angstroms (this is 21keV)


# the main callable function 
def convert_Nika_to_Fit2D(SSD, pix_size, BCX, BCY, HorTilt, VertTilt, wavelengtgh):
    '''
    Convert Nika data to Fit2D format.
    Needs pyFAI
    Convert Nika SDD, pix size, BCX, BCY, HorTilt, verTiilt
    foirst into Fit2D format
    Then use pyFAI to convert the data to Fit2D format.
    units: 
        SDD in mm, 
        pix_size in mm (convert below to um), 
        BCX, BCY in pixels, 
        HorTilt, VertTilt in degrees, 
        wavelength in A
    '''
    Fit2D_TiltDir, Fit2D_TiltAngle = convert_Nika_to_Fit2D_angles(HorTilt, VertTilt)     

    # Create a Fit2dGeometry object
    fit2d_geometry = Fit2dGeometry(
        directDist=SSD,
        centerX=BCY,
        centerY=BCX,
        tilt=Fit2D_TiltAngle,
        tiltPlanRotation=Fit2D_TiltDir,
        pixelX=pix_size*1000,
        pixelY=pix_size*1000,
        splineFile=None,
        detector=None,
        wavelength=wavelengtgh
    )
    #pp(fit2d_geometry)
    # Convert to PONI object
    poni_object = convert_from_Fit2d(fit2d_geometry)

    # Print the PONI object
    #pp(poni_object)
    return poni_object


def convert_Nika_to_Fit2D_angles(HorTilt, VertTilt):
    '''
    Convert angles from Nika format to Fit2D format.
    The conversion is based on the assumption that the angles are given in degrees.
    The conversion is done by negating the angles and converting them to radians.
    The function returns the converted angles in degrees.
    '''
    # Convert angles from degrees to radians
    Ax_rad = np.radians(HorTilt)
    Ay_rad = np.radians(VertTilt)
    
    # in any of the two tilts is negative, the tilt direction will be negative.
    sign = 1
    if Ax_rad < 0:
        sign = -1
    if Ay_rad < 0:
        sign = -1

    # Calculate the rotation vector components
    # Assuming small angles, the rotation vector can be approximated as:
    # Rx = [1, 0, 0] rotation around x-axis
    # Ry = [0, 1, 0] rotation around y-axis
    # Combined effect in the xy-plane
    Rx = np.array([0, np.sin(Ax_rad), np.cos(Ax_rad)])
    Ry = np.array([np.sin(Ay_rad), 0, np.cos(Ay_rad)])

    # Combined rotation vector
    R_combined = Rx + Ry

    # Calculate tilt direction angle (theta) in the xy-plane
    theta = sign*np.degrees(np.arctan2(R_combined[1], R_combined[0]))

     # Calculate tilt angle (phi) as the magnitude of the rotation vector
    phi = np.degrees(np.arcsin(np.linalg.norm([np.sin(Ax_rad), np.sin(Ay_rad)])))
    return theta, phi


# Test from measured data on 2025-04-05 WAXS
# # using Eigter 2M-W detector
# # normal orientation
# HorTilt = 0.037
# VertTilt = 0.088
# print(f"pyFAI/Fit2D Tilt Direction: 105.25 degrees, Tilt Angle: 0.56 degrees")
# theta, phi = convert_Nika_to_Fit2D_angles(HorTilt, VertTilt)
# print(f"Converted Tilt Direction: {theta:.2f} degrees, Tilt Angle: {phi:.2f} degrees")

# # LaB6_tilt1_0043
# HorTilt = 0.03
# VertTilt = 9.29
# print(f"pyFAI/Fit2D Tilt Direction: 0.995 degrees, Tilt Angle: 9.16 degrees")
# theta, phi = convert_Nika_to_Fit2D_angles(HorTilt, VertTilt)
# print(f"Converted Tilt Direction: {theta:.2f} degrees, Tilt Angle: {phi:.2f} degrees")

# # LaB6_tilt3h_0045
# HorTilt = -8.485
# VertTilt = -0.526
# print(f"pyFAI/Fit2D Tilt Direction: 90 degrees, Tilt Angle: 8.9744 degrees")
# theta, phi = convert_Nika_to_Fit2D_angles(HorTilt, VertTilt)
# print(f"Converted Tilt Direction: {theta:.2f} degrees, Tilt Angle: {phi:.2f} degrees")

# # LaB6_tilt5comb_0047
# HorTilt = -8.97
# VertTilt = 8.974
# print(f"pyFAI/Fit2D Tilt Direction: 45.7 degrees, Tilt Angle: 12.78 degrees")
# theta, phi = convert_Nika_to_Fit2D_angles(HorTilt, VertTilt)
# print(f"Converted Tilt Direction: {theta:.2f} degrees, Tilt Angle: {phi:.2f} degrees")

# # LaB6_tilt6comb_0048
# HorTilt = -6.713
# VertTilt = 5.441
# print(f"pyFAI/Fit2D Tilt Direction: 51.2 degrees, Tilt Angle: 8.52 degrees")
# theta, phi = convert_Nika_to_Fit2D_angles(HorTilt, VertTilt)
# print(f"Converted Tilt Direction: {theta:.2f} degrees, Tilt Angle: {phi:.2f} degrees")

# # LaB6_tilt7v
# HorTilt = 0.316
# VertTilt = -9.5
# print(f"pyFAI/Fit2D Tilt Direction: -178.9 degrees, Tilt Angle: 9.45 degrees")
# theta, phi = convert_Nika_to_Fit2D_angles(HorTilt, VertTilt)
# print(f"Converted Tilt Direction: {theta:.2f} degrees, Tilt Angle: {phi:.2f} degrees")

#poni = convert_Nika_to_Fit2D(402.53, 75, 325.03, 2946.74, 0.316, -9.489, 12.4/21)
#print(poni)