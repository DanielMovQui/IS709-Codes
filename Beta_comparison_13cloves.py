import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict, Counter
from itertools import combinations

a = 8.2 # Side of the squared front face of the clover in cm
b = 7.0 # Depth in cm of the HpGe crystals
d = 12.5 + b/2 # Distance in cm from implantation point to middle of front face of the clover + depth of HpGe crystals

# Set up 320
clovers_320 = {
    'g1_A': (np.radians(45),  np.radians(35)),
    'g1_B': (np.radians(135), np.radians(35)),
    'g2_A': (np.radians(45),  np.radians(90)),
    'g2_B': (np.radians(90),  np.radians(90)),
    'g2_C': (np.radians(135), np.radians(90)),
    'g3_A': (np.radians(45),  np.radians(180)),
    'g3_B': (np.radians(90),  np.radians(180)),
    'g3_C': (np.radians(135), np.radians(180)),
    'g4_A': (np.radians(45),  np.radians(270)),
    'g4_B': (np.radians(90),  np.radians(270)),
    'g4_C': (np.radians(135), np.radians(270)),
    'g5_A': (np.radians(45),  np.radians(320)),
    'g5_B': (np.radians(135), np.radians(320))
}

# Set up 325
clovers_325 = clovers_320.copy()
clovers_325['g5_A'] = (np.radians(45), np.radians(325))
clovers_325['g5_B'] = (np.radians(135), np.radians(325))

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

def angle_between(v1, v2):
    v1_u = v1 / np.linalg.norm(v1)
    v2_u = v2 / np.linalg.norm(v2)
    cos_angle = np.clip(np.dot(v1_u, v2_u), -1.0, 1.0) # scalar product
    angle = np.degrees(np.arccos(cos_angle))
    return angle, cos_angle

def calcular_cos_angles(clovers_dict):
    crystals_positions = {}
    original_colors = ['red', 'green', 'white', 'blue']

    for clover_label, (theta, phi) in clovers_dict.items():
        clover_center = spherical_to_cartesian(theta, phi)
        r_vector = clover_center / np.linalg.norm(clover_center)
        M = build_local_basis(r_vector)

        for color in original_colors:
            offset_local = crystal_offsets[color]
            offset_global = M @ offset_local
            crystal_pos = clover_center + offset_global
            crystals_positions[(clover_label, color)] = crystal_pos

    cos_angles = []
    for (key1, pos1), (key2, pos2) in combinations(crystals_positions.items(), 2):
        _, cos_angle = angle_between(pos1, pos2)
        cos_angles.append(cos_angle)

    return cos_angles

cos_angles_320 = calcular_cos_angles(clovers_320)
cos_angles_325 = calcular_cos_angles(clovers_325)

bins = np.linspace(-1, 1, 1001)  # 1000 bins entre -1 y 1
counts_320, bin_edges = np.histogram(cos_angles_320, bins=bins)
counts_325, _ = np.histogram(cos_angles_325, bins=bins)

plt.figure(figsize=(12, 6))
plt.hist(cos_angles_320, bins=1000, range=(-1, 1), color='blue', alpha=0.5, label='g5 at 320°')
plt.hist(cos_angles_325, bins=1000, range=(-1, 1), color='red', alpha=0.5, label='g5 at 325°')
plt.xlabel('cos(θ)')
plt.ylabel('Number of pairs')
plt.title('Pairs per cos(θ) with 13 clovers')
plt.legend(prop={'size': 16})
plt.grid(True)
plt.tight_layout()
plt.show()
