# THIS CODE GIVES THE MAP OF THE SET UP IN A .TXT FILE THAT THEN IS READ BY THE ids_hist_corr.C file 
# 10 clovers set up

import numpy as np
from collections import defaultdict, Counter
from itertools import combinations
import matplotlib.pyplot as plt


a = 8.2 # Side of the squared front face of the clover in cm
b = 7.0 # Depth in cm of the HpGe crystals
d = 16.24 + b/2 # Distance in cm from implantation point to middle of front face of the clover + depth of HpGe crystals

# clover angular distribution in spherical coordinates
clovers = {
    'g1_A': (np.radians(45),  np.radians(35)),
    'g1_B': (np.radians(135), np.radians(35)),
    'g2_A': (np.radians(45),  np.radians(90)),
    'g2_B': (np.radians(135),  np.radians(90)),
    'g3_A': (np.radians(45),  np.radians(180)),
    'g3_B': (np.radians(135),  np.radians(180)),
    'g4_A': (np.radians(45),  np.radians(270)),
    'g4_B': (np.radians(135),  np.radians(270)),
    'g5_A': (np.radians(45),  np.radians(325)),
    'g5_B': (np.radians(135), np.radians(325))
}

# Corrections to center of crytals in cartesian coordinates
crystal_offsets = {
    'red':   np.array([0,  a/4,  a/4]),
    'green': np.array([0,  a/4, -a/4]),
    'white': np.array([0, -a/4,  a/4]),
    'blue':  np.array([0, -a/4, -a/4])
}

# Function to change to cartesian coordinates
def spherical_to_cartesian(theta, phi, r=d):
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return np.array([x, y, z])

# Function to create local reference system for each clover
def build_local_basis(r_vector):
    r = r_vector / np.linalg.norm(r_vector)
    # Auxiliar up vector
    up = np.array([0, 0, 1])
    if np.allclose(np.abs(np.dot(r, up)), 1.0):
        up = np.array([0, 1, 0])
    # Calculation real up axis
    z_local = np.cross(r, up)
    z_local /= np.linalg.norm(z_local)
    # The third axis
    y_local = np.cross(z_local, r)
    M = np.column_stack((r, y_local, z_local)) # Local basis matrix
    return M

# ================================
# GENERATING POSITIONS OF CRYSTALS
# ================================

original_colors = ['red', 'green', 'white', 'blue']
crystals_positions = {}
# Iteration over all the clovers
for clover_label, (theta, phi) in clovers.items():
    clover_center = spherical_to_cartesian(theta, phi)
    r_vector = clover_center / np.linalg.norm(clover_center)
    M = build_local_basis(r_vector)

    # Iteration over all the crystals
    for color_idx, color in enumerate(original_colors):
        offset_local = crystal_offsets[original_colors[color_idx]]
        offset_global = M @ offset_local # Local offset to global coordinates
        crystal_pos = clover_center + offset_global # Final position of each clover
        crystals_positions[(clover_label, color)] = crystal_pos

# =======
# ANGLES
# =======

def angle_between(v1, v2):
    v1_u = v1 / np.linalg.norm(v1)
    v2_u = v2 / np.linalg.norm(v2)
    cos_angle = np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)
    angle = np.degrees(np.arccos(cos_angle))
    return angle

pairs = []
for (key1, pos1), (key2, pos2) in combinations(crystals_positions.items(), 2):
    angle = angle_between(pos1, pos2)
    pairs.append((key1, key2, angle))

# =================
# EXPORTING RESULTS 
# =================

with open('test_10clovers.txt', 'w', encoding='utf-8') as f:
    #f.write('Crystal pairs, Angle (º)\n')
    for (clover1, color1), (clover2, color2), angle in pairs:
        f.write(f"('{clover1}', '{color1}') - ('{clover2}', '{color2}'), {angle:.4f}\n")

