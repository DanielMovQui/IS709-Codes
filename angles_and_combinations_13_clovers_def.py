# THIS CODE GIVES THE MAP OF THE SET UP IN A .TXT FILE NOT VALID TO THE ids_histo_corr BUT MAKES A FANCY PLOT IN ORDER TO CHECK CRYSTAL PAIRING
# 13 clovers set up

import numpy as np
from collections import defaultdict, Counter
from itertools import combinations
#import matplotlib
#matplotlib.use('TkAgg')  
import matplotlib.pyplot as plt

# =====================
# DETECTOR GEOMETRY 
# =====================

a = 8.2 # Side of the squared front face of the clover in cm
b = 7.0 # Depth in cm of the HpGe crystals
d = 12.5 + b/2 # Distance in cm from implantation point to middle of front face of the clover + depth of HpGe crystals

# clover angular distribution in spherical coordinates
clovers = {
    'g1_A': (np.radians(45),  np.radians(36)),
    'g1_B': (np.radians(135), np.radians(36)),
    'g2_A': (np.radians(45),  np.radians(90)),
    'g2_B': (np.radians(90),  np.radians(90)),
    'g2_C': (np.radians(135), np.radians(90)),
    'g3_A': (np.radians(45),  np.radians(180)),
    'g3_B': (np.radians(90),  np.radians(180)),
    'g3_C': (np.radians(135), np.radians(180)),
    'g4_A': (np.radians(45),  np.radians(270)),
    'g4_B': (np.radians(90),  np.radians(270)),
    'g4_C': (np.radians(135), np.radians(270)),
    'g5_A': (np.radians(45),  np.radians(324)),
    'g5_B': (np.radians(135), np.radians(324))
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
"""
# ======
# ANGLES
# ======

def angle_between(v1, v2):
    v1_u = v1 / np.linalg.norm(v1)
    v2_u = v2 / np.linalg.norm(v2)
    cos_angle = np.clip(np.dot(v1_u, v2_u), -1.0, 1.0) # scalar product
    angle = np.degrees(np.arccos(cos_angle))
    return angle, cos_angle


angle_groups = defaultdict(list)
angles_list = []

for (key1, pos1), (key2, pos2) in combinations(crystals_positions.items(), 2):
    angle, _ = angle_between(pos1, pos2)
    angle_rounded = round(angle, 3)
    angle_groups[angle_rounded].append((key1, key2))
    angles_list.append(angle_rounded)

angle_counts = Counter(angles_list)

# =================
# EXPORTING RESULTS 
# =================

with open('test_13clovers_2.txt', 'w') as f:
    f.write("Angles between all the crystal pairs and combinations\n")
    f.write("---------------------------------------------------------\n")
    for angle in sorted(angle_counts):
        count = angle_counts[angle]
        f.write(f"Angle: {angle}°, Combinations: {count}\n")
        combos = angle_groups[angle]
        for (clover1, color1), (clover2, color2) in combos:
            f.write(f"  ('{clover1}', {color1}) - ('{clover2}', {color2})\n")
        f.write("\n")

print("Data exported to 'test_13clovers_2.txt'.")
"""

def angle_between(v1, v2):
    v1_u = v1 / np.linalg.norm(v1)
    v2_u = v2 / np.linalg.norm(v2)
    cos_angle = np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)
    angle = np.degrees(np.arccos(cos_angle))
    return angle, cos_angle

pairs = []
for (key1, pos1), (key2, pos2) in combinations(crystals_positions.items(), 2):
    angle = angle_between(pos1, pos2)
    pairs.append((key1, key2, angle))

# =================
# EXPORTING RESULTS 
# =================

with open('angles_and_combinations_13_clovers.txt', 'w', encoding='utf-8') as f:
    #f.write('Crystal pairs, Angle (º)\n')
    for (clover1, color1), (clover2, color2), angle in pairs:
        f.write(f"('{clover1}', '{color1}') - ('{clover2}', '{color2}'), {angle:.4f}\n")

# ====
# PLOT
# ====

cos_angles = []

for (key1, pos1), (key2, pos2) in combinations(crystals_positions.items(), 2):
    _, cos_angle = angle_between(pos1, pos2)
    cos_angles.append(cos_angle)


bins = np.linspace(-1, 1, 22)  # 22 bins between -1 and 1

# Plot
plt.figure(figsize=(10, 6))
plt.hist(cos_angles, bins=bins, edgecolor='black')
plt.xlim(-1,1)
plt.xlabel("cos(θ)")
plt.ylabel('Number of pairs')
plt.title('Pairs per cos(θ) with 13 clovers')
plt.grid(True)
plt.tight_layout()
plt.show()

        
