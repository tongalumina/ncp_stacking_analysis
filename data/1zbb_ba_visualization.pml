load 1zbb_ba_stack.pdb, main_obj
bg_color white
hide everything
show cartoon
util.cbc
set cartoon_transparency, 0.2
select ncp1, (chain q and resi 44-135) or (chain r and resi 24-102) or (chain s and resi 21-116) or (chain t and resi 34-121) or (chain u and resi 44-135) or (chain v and resi 24-102) or (chain w and resi 21-116) or (chain x and resi 34-121) or (chain Y and resi 179-325) or (chain Z and resi 23-169)
color palecyan, ncp1
select h3_pair_1, (chain q and resi 44-135) or (chain u and resi 44-135)
select ncp2, (chain A and resi 44-135) or (chain B and resi 24-102) or (chain C and resi 21-116) or (chain D and resi 34-121) or (chain E and resi 44-135) or (chain F and resi 24-102) or (chain G and resi 21-116) or (chain H and resi 34-121) or (chain I and resi 12-158) or (chain J and resi 190-336)
color lightorange, ncp2
select h3_pair_2, (chain A and resi 44-135) or (chain E and resi 44-135)
select ncp3, (chain Q and resi 44-135) or (chain R and resi 24-102) or (chain S and resi 21-116) or (chain T and resi 34-121) or (chain U and resi 44-135) or (chain V and resi 24-102) or (chain W and resi 21-116) or (chain X and resi 34-121) or (chain Y and resi 12-158) or (chain Z and resi 190-336)
color lightmagenta, ncp3
select h3_pair_3, (chain Q and resi 44-135) or (chain U and resi 44-135)
select ncp4, (chain a and resi 44-135) or (chain b and resi 24-102) or (chain c and resi 21-116) or (chain d and resi 34-121) or (chain e and resi 44-135) or (chain f and resi 24-102) or (chain g and resi 21-116) or (chain h and resi 34-121) or (chain I and resi 179-325) or (chain J and resi 23-169)
color palegreen, ncp4
select h3_pair_4, (chain a and resi 44-135) or (chain e and resi 44-135)
pseudoatom com1, pos=[81.296, 24.530, 191.869], color=palecyan, sphere_scale=1.0
cgo_arrow [81.296, 24.530, 191.869], [84.072, 19.776, 172.641], radius=0.3, color=red, name=axis1
cgo_arrow [81.296, 24.530, 191.869], [66.169, 11.493, 192.976], radius=0.3, color=blue, name=normal1
load_cgo([ BEGIN, TRIANGLE_FAN, COLOR, palecyan, ALPHA, 0.6, VERTEX, 81.296, 24.530, 191.869, VERTEX, 106.889, -4.246, 202.680, VERTEX, 105.721, -2.391, 208.561, VERTEX, 103.951, 0.127, 214.032, VERTEX, 101.623, 3.245, 218.956, VERTEX, 98.795, 6.888, 223.214, VERTEX, 95.536, 10.966, 226.700, VERTEX, 91.926, 15.377, 229.328, VERTEX, 88.055, 20.013, 231.034, VERTEX, 84.017, 24.761, 231.775, VERTEX, 79.912, 29.504, 231.534, VERTEX, 75.841, 34.123, 230.316, VERTEX, 71.904, 38.507, 228.152, VERTEX, 68.199, 42.546, 225.094, VERTEX, 64.816, 46.142, 221.218, VERTEX, 61.839, 49.206, 216.619, VERTEX, 59.341, 51.662, 211.411, VERTEX, 57.384, 53.450, 205.721, VERTEX, 56.015, 54.526, 199.691, VERTEX, 55.269, 54.863, 193.468, VERTEX, 55.164, 54.453, 187.205, VERTEX, 55.702, 53.307, 181.058, VERTEX, 56.870, 51.452, 175.176, VERTEX, 58.640, 48.934, 169.706, VERTEX, 60.968, 45.815, 164.781, VERTEX, 63.796, 42.173, 160.524, VERTEX, 67.055, 38.095, 157.038, VERTEX, 70.665, 33.684, 154.409, VERTEX, 74.536, 29.047, 152.704, VERTEX, 78.574, 24.299, 151.962, VERTEX, 82.679, 19.557, 152.203, VERTEX, 86.750, 14.937, 153.421, VERTEX, 90.687, 10.554, 155.586, VERTEX, 94.392, 6.515, 158.644, VERTEX, 97.775, 2.919, 162.520, VERTEX, 100.752, -0.145, 167.119, VERTEX, 103.250, -2.601, 172.327, VERTEX, 105.207, -4.389, 178.016, VERTEX, 106.576, -5.465, 184.047, VERTEX, 107.322, -5.802, 190.270, VERTEX, 107.427, -5.393, 196.532, VERTEX, 106.889, -4.246, 202.680, END ], plane1, 1)
pseudoatom com2, pos=[48.394, -23.320, 191.041], color=lightorange, sphere_scale=1.0
cgo_arrow [48.394, -23.320, 191.041], [48.232, -21.218, 171.153], radius=0.3, color=red, name=axis2
cgo_arrow [48.394, -23.320, 191.041], [63.587, -10.408, 192.614], radius=0.3, color=blue, name=normal2
load_cgo([ BEGIN, TRIANGLE_FAN, COLOR, lightorange, ALPHA, 0.6, VERTEX, 48.394, -23.320, 191.041, VERTEX, 22.380, 6.874, 194.445, VERTEX, 22.673, 5.778, 200.618, VERTEX, 23.598, 3.965, 206.556, VERTEX, 25.135, 1.481, 212.111, VERTEX, 27.244, -1.615, 217.148, VERTEX, 29.874, -5.244, 221.542, VERTEX, 32.960, -9.319, 225.185, VERTEX, 36.426, -13.738, 227.987, VERTEX, 40.187, -18.394, 229.879, VERTEX, 44.150, -23.170, 230.815, VERTEX, 48.217, -27.951, 230.772, VERTEX, 52.289, -32.617, 229.750, VERTEX, 56.265, -37.055, 227.776, VERTEX, 60.047, -41.154, 224.896, VERTEX, 63.542, -44.814, 221.184, VERTEX, 66.664, -47.945, 216.729, VERTEX, 69.337, -50.469, 211.641, VERTEX, 71.493, -52.325, 206.046, VERTEX, 73.081, -53.467, 200.082, VERTEX, 74.061, -53.866, 193.895, VERTEX, 74.409, -53.513, 187.638, VERTEX, 74.116, -52.417, 181.465, VERTEX, 73.191, -50.605, 175.527, VERTEX, 71.654, -48.120, 169.972, VERTEX, 69.545, -45.025, 164.935, VERTEX, 66.915, -41.395, 160.541, VERTEX, 63.829, -37.321, 156.898, VERTEX, 60.363, -32.901, 154.096, VERTEX, 56.602, -28.246, 152.204, VERTEX, 52.639, -23.469, 151.267, VERTEX, 48.572, -18.689, 151.311, VERTEX, 44.500, -14.022, 152.332, VERTEX, 40.524, -9.585, 154.307, VERTEX, 36.742, -5.486, 157.186, VERTEX, 33.247, -1.826, 160.899, VERTEX, 30.125, 1.305, 165.354, VERTEX, 27.452, 3.830, 170.442, VERTEX, 25.296, 5.686, 176.036, VERTEX, 23.708, 6.827, 182.001, VERTEX, 22.728, 7.227, 188.188, VERTEX, 22.380, 6.874, 194.445, END ], plane2, 1)
pseudoatom com3, pos=[79.280, -23.320, 46.085], color=lightmagenta, sphere_scale=1.0
cgo_arrow [79.280, -23.320, 46.085], [79.443, -21.218, 65.973], radius=0.3, color=red, name=axis3
cgo_arrow [79.280, -23.320, 46.085], [64.088, -10.408, 44.512], radius=0.3, color=blue, name=normal3
load_cgo([ BEGIN, TRIANGLE_FAN, COLOR, lightmagenta, ALPHA, 0.6, VERTEX, 79.280, -23.320, 46.085, VERTEX, 105.295, 6.874, 42.681, VERTEX, 105.002, 5.778, 36.508, VERTEX, 104.076, 3.965, 30.571, VERTEX, 102.540, 1.481, 25.015, VERTEX, 100.431, -1.615, 19.978, VERTEX, 97.801, -5.244, 15.585, VERTEX, 94.715, -9.319, 11.942, VERTEX, 91.249, -13.738, 9.140, VERTEX, 87.488, -18.394, 7.247, VERTEX, 83.525, -23.170, 6.311, VERTEX, 79.458, -27.951, 6.354, VERTEX, 75.386, -32.617, 7.376, VERTEX, 71.410, -37.054, 9.351, VERTEX, 67.628, -41.154, 12.230, VERTEX, 64.133, -44.814, 15.943, VERTEX, 61.011, -47.945, 20.398, VERTEX, 58.338, -50.469, 25.485, VERTEX, 56.182, -52.325, 31.080, VERTEX, 54.594, -53.467, 37.044, VERTEX, 53.614, -53.866, 43.231, VERTEX, 53.266, -53.513, 49.488, VERTEX, 53.558, -52.417, 55.662, VERTEX, 54.484, -50.605, 61.599, VERTEX, 56.021, -48.120, 67.155, VERTEX, 58.130, -45.025, 72.191, VERTEX, 60.760, -41.395, 76.585, VERTEX, 63.846, -37.321, 80.228, VERTEX, 67.312, -32.901, 83.030, VERTEX, 71.073, -28.246, 84.923, VERTEX, 75.036, -23.469, 85.859, VERTEX, 79.103, -18.689, 85.815, VERTEX, 83.175, -14.022, 84.794, VERTEX, 87.151, -9.585, 82.819, VERTEX, 90.933, -5.486, 79.940, VERTEX, 94.428, -1.826, 76.227, VERTEX, 97.550, 1.305, 71.772, VERTEX, 100.222, 3.830, 66.685, VERTEX, 102.379, 5.686, 61.090, VERTEX, 103.967, 6.827, 55.126, VERTEX, 104.947, 7.227, 48.939, VERTEX, 105.295, 6.874, 42.681, END ], plane3, 1)
pseudoatom com4, pos=[46.379, 24.530, 45.257], color=palegreen, sphere_scale=1.0
cgo_arrow [46.379, 24.530, 45.257], [43.603, 19.776, 64.484], radius=0.3, color=red, name=axis4
cgo_arrow [46.379, 24.530, 45.257], [61.506, 11.493, 44.150], radius=0.3, color=blue, name=normal4
load_cgo([ BEGIN, TRIANGLE_FAN, COLOR, palegreen, ALPHA, 0.6, VERTEX, 46.379, 24.530, 45.257, VERTEX, 20.786, -4.246, 34.446, VERTEX, 21.954, -2.391, 28.565, VERTEX, 23.724, 0.127, 23.094, VERTEX, 26.052, 3.245, 18.170, VERTEX, 28.880, 6.888, 13.912, VERTEX, 32.139, 10.965, 10.426, VERTEX, 35.749, 15.377, 7.798, VERTEX, 39.620, 20.013, 6.092, VERTEX, 43.658, 24.761, 5.351, VERTEX, 47.763, 29.503, 5.592, VERTEX, 51.834, 34.123, 6.809, VERTEX, 55.771, 38.507, 8.974, VERTEX, 59.476, 42.546, 12.032, VERTEX, 62.859, 46.142, 15.908, VERTEX, 65.836, 49.206, 20.507, VERTEX, 68.334, 51.662, 25.715, VERTEX, 70.291, 53.450, 31.404, VERTEX, 71.660, 54.526, 37.435, VERTEX, 72.406, 54.863, 43.658, VERTEX, 72.511, 54.453, 49.920, VERTEX, 71.973, 53.307, 56.068, VERTEX, 70.805, 51.452, 61.949, VERTEX, 69.035, 48.934, 67.420, VERTEX, 66.707, 45.815, 72.345, VERTEX, 63.879, 42.173, 76.602, VERTEX, 60.620, 38.095, 80.088, VERTEX, 57.010, 33.684, 82.716, VERTEX, 53.138, 29.047, 84.422, VERTEX, 49.101, 24.300, 85.164, VERTEX, 44.996, 19.557, 84.923, VERTEX, 40.925, 14.938, 83.705, VERTEX, 36.988, 10.554, 81.540, VERTEX, 33.283, 6.515, 78.482, VERTEX, 29.900, 2.919, 74.606, VERTEX, 26.923, -0.145, 70.008, VERTEX, 24.425, -2.601, 64.799, VERTEX, 22.468, -4.389, 59.110, VERTEX, 21.099, -5.465, 53.080, VERTEX, 20.353, -5.802, 46.856, VERTEX, 20.248, -5.393, 40.594, VERTEX, 20.786, -4.246, 34.446, END ], plane4, 1)
cgo_arrow [81.296, 24.530, 191.869], [48.394, -23.320, 191.041], radius=0.3, color=black, name=dist_vec_1_2
pseudoatom label_title_1_2, pos=[104.845, 0.605, 191.455]
label label_title_1_2, "NCP 1-NCP 2 Stack" 
pseudoatom label_1_2_Distance, pos=[104.845, -7.395, 191.455]
label label_1_2_Distance, "Distance: 58.1 A"
pseudoatom label_1_2_Rise, pos=[104.845, -13.395, 191.455]
label label_1_2_Rise, "Rise: 56.0 A"
pseudoatom label_1_2_Shift, pos=[104.845, -19.395, 191.455]
label label_1_2_Shift, "Shift: 15.3 A"
pseudoatom label_1_2_Shift_Orientation, pos=[104.845, -25.395, 191.455]
label label_1_2_Shift_Orientation, "Shift Orientation: 59.4 deg"
pseudoatom label_1_2_Symmetry_Axes_Orientation, pos=[104.845, -31.395, 191.455]
label label_1_2_Symmetry_Axes_Orientation, "Symmetry Axes Orientation: -20.6 deg"
pseudoatom label_1_2_Tilt, pos=[104.845, -37.395, 191.455]
label label_1_2_Tilt, "Tilt: 7.7 deg"
pseudoatom label_1_2_Tilt_Direction, pos=[104.845, -43.395, 191.455]
label label_1_2_Tilt_Direction, "Tilt Direction: -18.7 deg"
cgo_arrow [48.394, -23.320, 191.041], [79.280, -23.320, 46.085], radius=0.3, color=black, name=dist_vec_2_3
pseudoatom label_title_2_3, pos=[63.837, -13.320, 118.563]
label label_title_2_3, "NCP 2-NCP 3 Stack" 
pseudoatom label_2_3_Distance, pos=[63.837, -21.320, 118.563]
label label_2_3_Distance, "Distance: 148.2 A"
pseudoatom label_2_3_Rise, pos=[63.837, -27.320, 118.563]
label label_2_3_Rise, "Rise: 12.1 A"
pseudoatom label_2_3_Shift, pos=[63.837, -33.320, 118.563]
label label_2_3_Shift, "Shift: 147.7 A"
pseudoatom label_2_3_Shift_Orientation, pos=[63.837, -39.320, 118.563]
label label_2_3_Shift_Orientation, "Shift Orientation: -12.7 deg"
pseudoatom label_2_3_Symmetry_Axes_Orientation, pos=[63.837, -45.320, 118.563]
label label_2_3_Symmetry_Axes_Orientation, "Symmetry Axes Orientation: 170.8 deg"
pseudoatom label_2_3_Tilt, pos=[63.837, -51.320, 118.563]
label label_2_3_Tilt, "Tilt: 80.4 deg"
pseudoatom label_2_3_Tilt_Direction, pos=[63.837, -57.320, 118.563]
label label_2_3_Tilt_Direction, "Tilt Direction: -98.7 deg"
cgo_arrow [79.280, -23.320, 46.085], [46.379, 24.530, 45.257], radius=0.3, color=black, name=dist_vec_3_4
pseudoatom label_title_3_4, pos=[102.830, -59.395, 45.671]
label label_title_3_4, "NCP 3-NCP 4 Stack" 
pseudoatom label_3_4_Distance, pos=[102.830, -67.395, 45.671]
label label_3_4_Distance, "Distance: 58.1 A"
pseudoatom label_3_4_Rise, pos=[102.830, -73.395, 45.671]
label label_3_4_Rise, "Rise: 55.9 A"
pseudoatom label_3_4_Shift, pos=[102.830, -79.395, 45.671]
label label_3_4_Shift, "Shift: 15.6 A"
pseudoatom label_3_4_Shift_Orientation, pos=[102.830, -85.395, 45.671]
label label_3_4_Shift_Orientation, "Shift Orientation: 71.8 deg"
pseudoatom label_3_4_Symmetry_Axes_Orientation, pos=[102.830, -91.395, 45.671]
label label_3_4_Symmetry_Axes_Orientation, "Symmetry Axes Orientation: -20.7 deg"
pseudoatom label_3_4_Tilt, pos=[102.830, -97.395, 45.671]
label label_3_4_Tilt, "Tilt: 7.7 deg"
pseudoatom label_3_4_Tilt_Direction, pos=[102.830, -103.395, 45.671]
label label_3_4_Tilt_Direction, "Tilt Direction: -1.9 deg"
set label_color, black, label_*
set label_size, 16
set label_font_id, 7
hide labels, com*
group basis_vectors, com* axis* normal* plane*
group h3_pairs, h3_pair_*
group labels, label_*
group stacking_vectors, dist_vec_*
zoom visible