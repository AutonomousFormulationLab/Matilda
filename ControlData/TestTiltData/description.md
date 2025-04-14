Tilt test data
These are the tilt test data collected on USAXS instrument in APril 2025.
These data have various levesl of tils and are used to test the tilt correction algorithm.
They were calibrated using Nika package 

Here are some notes from testing:

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

Using data collected April 5, 2025

normal LaB6
LaB6_Normal_0042.hdf
Nika:
SDD= 417.13 mm
BCX = 271.59 pix
BCY = 2833.76 pix
HorTilt = 0.037 deg
VerTilt = 0.088 deg

Fit2D:
SDD = 416.23mm
BCX = 2832.945 px
BCY = 268.159 px
Tilt = 0.5600 deg
Tilt Plan rot = 105.2495 deg. 

 
LaB6_tilt1_0043.hdf
Nika:
SDD = 430.00 mm
BCX = 287.65
BCY = 2763.93
Hor tilt = 0.03
Vert tilt = 9.293

Fit2D:
SDD= 429.682
BCX = 2763.966
BCY = 287.85
Tilt = 9.159 deg
Tilt plan rotation = 0.955 deg. 





LaB6_tilt3h_0045.hdf (tilt in X)
Nika:
SDD = 406.04
BCX = 206.50
BCY = 2834.63
Tilt hor = -8.485
Tilt Vert = -0.526

Fit2D:
SDD = 408.477 mm
BCX = 2840.07
BCY = 203.876
Tilt = 8.9744
Tilt plan rot: 90 deg





LaB6_tilt5_comb_0047 combined tilt:
Nika
SDD = 420.66
BCX = 230.95
BCY = 2759.61
Tilt Hor = -8.968
Tilt Vert = 8.974

Fit2D:
SDD =. 420.53
BCX = 2759.91
BCY = 230.144
Tilt = 12.7822
Tilt plan rot: 45.704






LaB6_tilt6_comb_0048
Nika:
LaB6 tilt 6 comb

SDD = 417.69
BCX = 234.03
BCY = 2781.47
Tilt hor = -6.713
Tlt vert = 5.441

Fit2D:
SDD = 417.3
BCX = 2781.19
BCY = 234.4
Tilt = 8.52
Tilt plane = 51.2 deg







LaB6 tilt 7 vert
Nika
SDD = 402.53
BCX = 325.03
BCY = 2946.74
Tilt hor = 0.316
Tilt vert = -9.489

Fit2D:
SDD = 402.9
BCX = 2948.14
BCY = 324.2
Tilt = 9.45
Tilt plane = -178.942
	- eq. to tilt -9.45 deg and plane 1deg
